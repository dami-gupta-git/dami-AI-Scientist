"""
Nonlinear (MLP) stability probe on S1724 — companion to megascale_stability.py.

Runs MLP regression (1280→256→64→1) under the same three CV schemes
(random / protein-holdout / cluster-holdout) and 5 seeds as the Ridge probe,
then compares to check whether nonlinearity adds signal beyond Ridge.

Same question as result_3/5/7 for mechanism: does the MLP lift survive
family-holdout, or does it evaporate (leakage)?

Usage:
  cd esm2_mechanism
  python scripts/megascale_mlp.py

Outputs:
  results/megascale_stability/mlp_summary.json
  results/megascale_stability/mlp_per_protein_spearman.json
"""

import json
import os
import sys
import numpy as np
from scipy.stats import spearmanr, pearsonr
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from megascale_stability import (
    load_s1724_variants, assign_protein_clusters,
    random_split_cv, protein_split_cv, cluster_split_cv,
    per_protein_spearman,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
EMB  = os.path.join(DATA, "embeddings")
OUT  = os.path.join(ROOT, "results", "megascale_stability")

WT_MEAN_EMB  = os.path.join(EMB, "megascale_wt_mean.npy")
MUT_MEAN_EMB = os.path.join(EMB, "megascale_mut_mean.npy")

N_SEEDS = 5
N_FOLDS = 5

os.makedirs(OUT, exist_ok=True)


# ---------------------------------------------------------------------------
# MLP regression probe
# ---------------------------------------------------------------------------

def run_mlp_regression(X, y, splits, seed=42, hidden=(256, 64),
                        lr=1e-3, max_epochs=200, patience=15, batch_size=64):
    import torch
    import torch.nn as nn
    from torch.utils.data import DataLoader, TensorDataset

    rhos, rs, aurocs = [], [], []

    for fold_i, (tr, te) in enumerate(splits):
        X_tr = X[tr].astype(np.float32)
        X_te = X[te].astype(np.float32)
        y_tr = y[tr].astype(np.float32)
        y_te = y[te].astype(np.float32)

        # hold out 15% of train for early stopping
        rng = np.random.RandomState(seed + fold_i)
        idx = np.arange(len(X_tr))
        rng.shuffle(idx)
        n_val = max(1, int(0.15 * len(idx)))
        val_idx, fit_idx = idx[:n_val], idx[n_val:]

        mu = X_tr[fit_idx].mean(0); std = X_tr[fit_idx].std(0) + 1e-8
        X_fit = (X_tr[fit_idx] - mu) / std
        X_val = (X_tr[val_idx] - mu) / std
        X_te_n = (X_te - mu) / std

        layers = []
        prev = X_fit.shape[1]
        for h in hidden:
            layers += [nn.Linear(prev, h), nn.ReLU(), nn.Dropout(0.2)]
            prev = h
        layers.append(nn.Linear(prev, 1))
        model = nn.Sequential(*layers)

        opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-4)
        crit = nn.MSELoss()

        ds = TensorDataset(torch.tensor(X_fit), torch.tensor(y_tr[fit_idx]).unsqueeze(1))
        loader = DataLoader(ds, batch_size=batch_size, shuffle=True)

        best_val, patience_cnt, best_state = float("inf"), 0, None
        for epoch in range(max_epochs):
            model.train()
            for xb, yb in loader:
                opt.zero_grad()
                crit(model(xb), yb).backward()
                opt.step()
            model.eval()
            with torch.no_grad():
                vl = crit(model(torch.tensor(X_val)),
                          torch.tensor(y_tr[val_idx]).unsqueeze(1)).item()
            if vl < best_val - 1e-4:
                best_val, patience_cnt = vl, 0
                best_state = {k: v.clone() for k, v in model.state_dict().items()}
            else:
                patience_cnt += 1
                if patience_cnt >= patience:
                    break
        if best_state:
            model.load_state_dict(best_state)

        model.eval()
        with torch.no_grad():
            pred = model(torch.tensor(X_te_n)).squeeze(1).numpy()

        rho, _ = spearmanr(y_te, pred)
        r, _   = pearsonr(y_te, pred)
        med = np.median(y_te)
        binary = (y_te >= med).astype(int)
        au = float(roc_auc_score(binary, pred)) if binary.sum() > 0 and (1-binary).sum() > 0 else float("nan")

        rhos.append(float(rho))
        rs.append(float(r))
        aurocs.append(au)

    if not rhos:
        return {}
    return {
        "spearman_mean": float(np.mean(rhos)),
        "spearman_std":  float(np.std(rhos)),
        "pearson_mean":  float(np.mean(rs)),
        "pearson_std":   float(np.std(rs)),
        "auroc_mean":    float(np.nanmean(aurocs)),
        "auroc_std":     float(np.nanstd(aurocs)),
        "n_folds": len(rhos),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    variants = load_s1724_variants()
    proteins  = np.array([v["protein"] for v in variants])
    ddg       = np.array([v["ddg"]     for v in variants])
    cluster_map = assign_protein_clusters(variants)

    print(f"Loaded {len(variants)} variants across {len(set(proteins))} proteins")

    wt_mean  = np.load(WT_MEAN_EMB)
    mut_mean = np.load(MUT_MEAN_EMB)
    delta_mean = mut_mean - wt_mean
    print(f"Embeddings: {delta_mean.shape}")

    results_by_seed = []
    for seed in range(N_SEEDS):
        print(f"\n── Seed {seed} ──")
        splits_random  = random_split_cv(len(variants), N_FOLDS, seed)
        splits_protein = protein_split_cv(proteins, N_FOLDS, seed)
        splits_cluster = cluster_split_cv(proteins, cluster_map, N_FOLDS, seed)

        seed_result = {"seed": seed}
        for split_name, splits in [
            ("random",  splits_random),
            ("protein", splits_protein),
            ("cluster", splits_cluster),
        ]:
            res = run_mlp_regression(delta_mean, ddg, splits, seed=seed)
            seed_result[split_name] = res
            if res:
                print(f"  {split_name}: ρ={res['spearman_mean']:.3f}±{res['spearman_std']:.3f}  "
                      f"AUROC={res['auroc_mean']:.3f}")
        results_by_seed.append(seed_result)

    # per-protein leave-one-out
    print("\nPer-protein Spearman (leave-one-out, MLP)...")
    per_prot_mlp = {}
    for prot in sorted(set(proteins)):
        mask = (proteins == prot)
        if mask.sum() < 5:
            continue
        tr = np.where(~mask)[0]
        te = np.where(mask)[0]
        if len(tr) < 10:
            continue
        splits_loo = [(tr, te)]
        res = run_mlp_regression(delta_mean, ddg, splits_loo, seed=0)
        if res:
            per_prot_mlp[prot] = {"spearman": res["spearman_mean"], "n_variants": int(mask.sum())}

    prot_rhos = [v["spearman"] for v in per_prot_mlp.values()]
    print(f"  mean={np.mean(prot_rhos):.3f}  std={np.std(prot_rhos):.3f}  "
          f"min={min(prot_rhos):.3f}  max={max(prot_rhos):.3f}  n={len(prot_rhos)}")

    # aggregate
    summary = {}
    for split_name in ["random", "protein", "cluster"]:
        vals_rho   = [sr[split_name]["spearman_mean"] for sr in results_by_seed if sr.get(split_name)]
        vals_auroc = [sr[split_name]["auroc_mean"]    for sr in results_by_seed if sr.get(split_name)]
        summary[f"mlp_{split_name}"] = {
            "spearman_mean": float(np.mean(vals_rho)),
            "spearman_std":  float(np.std(vals_rho)),
            "auroc_mean":    float(np.nanmean(vals_auroc)),
            "auroc_std":     float(np.nanstd(vals_auroc)),
        }

    summary["mlp_per_protein"] = {
        "spearman_mean": float(np.mean(prot_rhos)),
        "spearman_std":  float(np.std(prot_rhos)),
        "n_proteins": len(prot_rhos),
    }

    # comparison table
    print(f"\n{'='*60}")
    print("MLP vs Ridge comparison (delta_mean):")
    print(f"  {'':20s}  {'random ρ':>10}  {'protein ρ':>10}  {'Δ':>8}  {'random AUROC':>13}  {'protein AUROC':>14}")
    for label, rnd_rho, prot_rho, rnd_au, prot_au in [
        ("Ridge (result_21)", 0.546, 0.280, 0.764, 0.642),
        ("MLP",
         summary["mlp_random"]["spearman_mean"],
         summary["mlp_protein"]["spearman_mean"],
         summary["mlp_random"]["auroc_mean"],
         summary["mlp_protein"]["auroc_mean"]),
    ]:
        print(f"  {label:20s}  {rnd_rho:>10.3f}  {prot_rho:>10.3f}  {rnd_rho-prot_rho:>8.3f}  {rnd_au:>13.3f}  {prot_au:>14.3f}")
    print(f"{'='*60}")

    with open(os.path.join(OUT, "mlp_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    with open(os.path.join(OUT, "mlp_per_protein_spearman.json"), "w") as f:
        json.dump(per_prot_mlp, f, indent=2)
    print(f"\nResults written to {OUT}/")


if __name__ == "__main__":
    main()
