"""
Megascale stability as a second ESM-2 positive control.

Tests whether ESM-2 delta embeddings predict ΔΔG (thermodynamic stability)
under random / protein / family-split CV, using the S1724 benchmark
(1,422 single-point missense mutations across 33 proteins).

Pre-registered hypotheses:
  H1: Spearman ρ ≥ 0.5 under random split (stability encoded)
  H2: ρ drops ≤ 0.05 under protein/family-split (family-robust)
  H3: Per-protein Spearman std ≤ 0.10 (tight distribution)

Decision table (plan_megascale_stability.md):
  ROBUST: random ρ ≥ 0.5, protein-split Δ ≤ 0.05, per-protein std ≤ 0.10
  WEAK:   random ρ 0.3–0.5
  HETEROGENEOUS: ρ ≥ 0.5, Δ ≤ 0.05, std ≥ 0.15
  LEAKY:  ρ ≥ 0.5, protein-split Δ ≥ 0.10
  NULL:   ρ < 0.3

Usage (GPU required for embedding extraction):
  cd esm2_mechanism
  python scripts/megascale_stability.py

Outputs:
  data/megascale_variants.json
  data/embeddings/megascale_wt_mean.npy
  data/embeddings/megascale_mut_mean.npy
  data/embeddings/megascale_wt_pos.npy
  data/embeddings/megascale_mut_pos.npy
  results/megascale_stability/summary.json
  results/megascale_stability/per_protein_spearman.json
"""

import json
import os
import sys
import zipfile
import numpy as np
from scipy.stats import spearmanr, pearsonr
from sklearn.linear_model import Ridge
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from experiment import get_esm2_embeddings_for_pairs, window_sequence, apply_missense, ESM2_MODEL_650M

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
EMB  = os.path.join(DATA, "embeddings")
BM_ZIP = os.path.join(DATA, "megascale", "benchmarks.zip")
PFAM_JSON = os.path.join(DATA, "pfam_families.json")
OUT = os.path.join(ROOT, "results", "megascale_stability")

VARIANTS_CACHE = os.path.join(DATA, "megascale_variants.json")
WT_MEAN_EMB  = os.path.join(EMB, "megascale_wt_mean.npy")
MUT_MEAN_EMB = os.path.join(EMB, "megascale_mut_mean.npy")
WT_POS_EMB   = os.path.join(EMB, "megascale_wt_pos.npy")
MUT_POS_EMB  = os.path.join(EMB, "megascale_mut_pos.npy")

N_SEEDS = 5
N_FOLDS = 5

os.makedirs(OUT, exist_ok=True)
os.makedirs(EMB, exist_ok=True)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_s1724_variants():
    """
    Parse S1724 benchmark (ThermoMutDB curated, single-point missense).

    WT sequence: S1724_PDB_sequences.csv gabetrimmed column (keyed by PDB_chain e.g. '1AJ3_A').
    Mutation position: paper_seq_muts column, 1-indexed within gabetrimmed.
    Returns list of dicts: {protein, mutation_code, wt_seq, mut_seq, var_pos, ddg}
    """
    if os.path.exists(VARIANTS_CACHE):
        print("Loading cached S1724 variants...")
        with open(VARIANTS_CACHE) as f:
            return json.load(f)

    import re
    import pandas as pd
    mut_pat = re.compile(r"^([A-Z])(\d+)([A-Z])$")

    with zipfile.ZipFile(BM_ZIP) as z:
        with z.open("benchmarks/S1724_thermomutdb_cleaned_withseq.csv") as f:
            df = pd.read_csv(f)

    # Single-point missense with ΔΔG and parseable mutation
    # paper_seq stores the MUTANT sequence; paper_seq_muts gives the 1-indexed position
    # within paper_seq where the mutation was applied (paper_seq[var_pos-1] == aa_mut).
    # WT is reconstructed by reversing the mutation: replace aa_mut -> aa_wt at var_pos.
    df = df[
        (df["mutation_type"] == "Single") &
        df["ddg"].notna() &
        df["paper_seq"].notna() &
        df["paper_seq_muts"].notna()
    ].copy()

    variants = []
    skipped = {"bad_mut_pat": 0, "pos_out_of_range": 0, "mut_not_in_paper_seq": 0}
    for _, row in df.iterrows():
        protein = str(row["PDB_wild"]).strip()
        mut_str = str(row["paper_seq_muts"]).strip()
        m = mut_pat.match(mut_str)
        if not m:
            skipped["bad_mut_pat"] += 1
            continue
        aa_wt, pos_str, aa_mut = m.groups()
        var_pos = int(pos_str)  # 1-indexed in paper_seq

        paper_seq = str(row["paper_seq"]).strip()
        if var_pos < 1 or var_pos > len(paper_seq):
            skipped["pos_out_of_range"] += 1
            continue
        if paper_seq[var_pos - 1] != aa_mut:
            # paper_seq should have the mutant residue at var_pos
            skipped["mut_not_in_paper_seq"] += 1
            continue

        # Reconstruct WT by reversing the mutation
        wt_seq  = paper_seq[:var_pos - 1] + aa_wt + paper_seq[var_pos:]
        mut_seq = paper_seq  # already the mutant

        variants.append({
            "protein": protein,
            "mutation_code": mut_str,
            "wt_seq": wt_seq,
            "mut_seq": mut_seq,
            "var_pos": var_pos,
            "ddg": float(row["ddg"]),
        })

    print(f"S1724: {len(variants)} single-point variants, skipped={skipped}")

    with open(VARIANTS_CACHE, "w") as f:
        json.dump(variants, f)
    print(f"Cached to {VARIANTS_CACHE}")
    return variants


# ---------------------------------------------------------------------------
# Protein-level clustering for family-split (MMseqs2 or sequence-identity)
# ---------------------------------------------------------------------------

def assign_protein_clusters(variants):
    """
    Group proteins by sequence identity for family-split CV.
    Uses exact-match on WT sequence as cluster ID (each unique WT = one protein domain).
    With only 33 proteins this is already a meaningful holdout.
    Returns: dict protein_id -> cluster_id (here: just protein_id itself, 1 cluster per protein)
    """
    # Each PDB_wild is one protein — cluster = protein identity.
    # For the family-split analogue we use MMseqs2 if available, else identity clustering.
    mmseqs_cache = os.path.join(DATA, "megascale_protein_clusters.json")
    if os.path.exists(mmseqs_cache):
        with open(mmseqs_cache) as f:
            return json.load(f)

    # Build unique WT sequences
    proteins = {}
    for v in variants:
        pid = v["protein"]
        if pid not in proteins:
            proteins[pid] = v["wt_seq"]

    # Try MMseqs2
    try:
        cluster_map = _run_mmseqs2(proteins)
    except Exception as e:
        print(f"MMseqs2 failed ({e}), falling back to identity clustering")
        cluster_map = {pid: pid for pid in proteins}

    with open(mmseqs_cache, "w") as f:
        json.dump(cluster_map, f)
    return cluster_map


def _run_mmseqs2(proteins, min_seq_id=0.20, coverage=0.20):
    """
    Run MMseqs2 easy-cluster on the 33 S1724 proteins.
    Returns dict: protein_id -> cluster_representative_id
    """
    import subprocess, tempfile, shutil

    if not shutil.which("mmseqs"):
        raise RuntimeError("mmseqs not found in PATH")

    with tempfile.TemporaryDirectory() as tmp:
        fasta = os.path.join(tmp, "seqs.fasta")
        with open(fasta, "w") as f:
            for pid, seq in proteins.items():
                f.write(f">{pid}\n{seq}\n")

        subprocess.run(
            ["mmseqs", "easy-cluster", fasta, os.path.join(tmp, "clust"),
             os.path.join(tmp, "mmseqs_tmp"),
             f"--min-seq-id", str(min_seq_id),
             "-c", str(coverage),
             "--cov-mode", "0", "-v", "0"],
            check=True, capture_output=True
        )

        tsv = os.path.join(tmp, "clust_cluster.tsv")
        cluster_map = {}
        with open(tsv) as f:
            for line in f:
                rep, member = line.strip().split("\t")
                cluster_map[member] = rep

    return cluster_map


# ---------------------------------------------------------------------------
# Embedding extraction
# ---------------------------------------------------------------------------

def extract_embeddings(variants):
    import torch
    device = "cuda" if torch.cuda.is_available() else "cpu"

    wt_seqs, mut_seqs, positions = [], [], []
    for v in variants:
        wt_win, new_pos = window_sequence(v["wt_seq"], v["var_pos"])
        # mut_seq is already built; re-derive windowed version for consistency
        mut_win = wt_win[:new_pos - 1] + v["mut_seq"][v["var_pos"] - 1] + wt_win[new_pos:]
        wt_seqs.append(wt_win)
        mut_seqs.append(mut_win)
        positions.append(new_pos)

    print(f"Extracting ESM-2 embeddings for {len(wt_seqs)} S1724 pairs on {device}...")
    wt_mean, mut_mean, wt_pos, mut_pos = get_esm2_embeddings_for_pairs(
        wt_seqs, mut_seqs, positions,
        model_name=ESM2_MODEL_650M, device=device, batch_size=64
    )
    np.save(WT_MEAN_EMB, wt_mean)
    np.save(MUT_MEAN_EMB, mut_mean)
    np.save(WT_POS_EMB, wt_pos)
    np.save(MUT_POS_EMB, mut_pos)
    print(f"Saved embeddings: shape {wt_mean.shape}")
    return wt_mean, mut_mean, wt_pos, mut_pos


# ---------------------------------------------------------------------------
# CV splits
# ---------------------------------------------------------------------------

def random_split_cv(n, n_folds=5, seed=42):
    idx = np.arange(n)
    np.random.RandomState(seed).shuffle(idx)
    splits = []
    for fold in np.array_split(idx, n_folds):
        te = fold
        tr = np.setdiff1d(idx, fold)
        splits.append((tr, te))
    return splits


def protein_split_cv(proteins, n_folds=5, seed=42):
    """Hold out whole proteins (analogous to gene-split)."""
    unique = np.array(sorted(set(proteins)))
    np.random.RandomState(seed).shuffle(unique)
    splits = []
    for fold_proteins in np.array_split(unique, n_folds):
        mask_te = np.isin(proteins, fold_proteins)
        tr = np.where(~mask_te)[0]
        te = np.where(mask_te)[0]
        if len(tr) >= 10 and len(te) >= 5:
            splits.append((tr, te))
    return splits


def cluster_split_cv(proteins, cluster_map, n_folds=5, seed=42):
    """Hold out whole MMseqs2/identity clusters (analogous to family-split)."""
    prot_clusters = np.array([cluster_map.get(p, p) for p in proteins])
    unique_clusters = np.array(sorted(set(prot_clusters)))
    np.random.RandomState(seed).shuffle(unique_clusters)
    splits = []
    for fold_clusters in np.array_split(unique_clusters, n_folds):
        mask_te = np.isin(prot_clusters, fold_clusters)
        tr = np.where(~mask_te)[0]
        te = np.where(mask_te)[0]
        if len(tr) >= 10 and len(te) >= 5:
            splits.append((tr, te))
    return splits


# ---------------------------------------------------------------------------
# Ridge regression probe
# ---------------------------------------------------------------------------

def run_ridge(X, y, splits):
    rhos, rs = [], []
    for tr, te in splits:
        sc = StandardScaler()
        Xtr = sc.fit_transform(X[tr])
        Xte = sc.transform(X[te])
        clf = Ridge(alpha=1.0)
        clf.fit(Xtr, y[tr])
        pred = clf.predict(Xte)
        rho, _ = spearmanr(y[te], pred)
        r, _   = pearsonr(y[te], pred)
        rhos.append(float(rho))
        rs.append(float(r))
    if not rhos:
        return {}
    return {
        "spearman_mean": float(np.mean(rhos)),
        "spearman_std":  float(np.std(rhos)),
        "pearson_mean":  float(np.mean(rs)),
        "pearson_std":   float(np.std(rs)),
        "n_folds": len(rhos),
    }


def auroc_at_median(y_true, y_pred):
    """Binary AUROC: above-median = positive."""
    med = np.median(y_true)
    binary = (y_true >= med).astype(int)
    if binary.sum() == 0 or (1 - binary).sum() == 0:
        return float("nan")
    return float(roc_auc_score(binary, y_pred))


def run_ridge_with_auroc(X, y, splits):
    rhos, rs, aurocs = [], [], []
    for tr, te in splits:
        sc = StandardScaler()
        Xtr = sc.fit_transform(X[tr])
        Xte = sc.transform(X[te])
        clf = Ridge(alpha=1.0)
        clf.fit(Xtr, y[tr])
        pred = clf.predict(Xte)
        rho, _ = spearmanr(y[te], pred)
        r, _   = pearsonr(y[te], pred)
        au     = auroc_at_median(y[te], pred)
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
# Per-protein Spearman distribution
# ---------------------------------------------------------------------------

def per_protein_spearman(X, y, proteins, use_delta_mean=True):
    """
    For each protein with ≥5 variants, fit Ridge on all others, predict on that protein.
    Analogous to result_17/18 per-stratum AUROC distributions.
    """
    unique = sorted(set(proteins))
    results = {}
    for prot in unique:
        mask = (proteins == prot)
        if mask.sum() < 5:
            continue
        tr = np.where(~mask)[0]
        te = np.where(mask)[0]
        if len(tr) < 10:
            continue
        sc = StandardScaler()
        Xtr = sc.fit_transform(X[tr])
        Xte = sc.transform(X[te])
        clf = Ridge(alpha=1.0)
        clf.fit(Xtr, y[tr])
        pred = clf.predict(Xte)
        rho, pval = spearmanr(y[te], pred)
        results[prot] = {
            "spearman": float(rho),
            "p_value": float(pval),
            "n_variants": int(mask.sum()),
        }
    return results


# ---------------------------------------------------------------------------
# Decision rule
# ---------------------------------------------------------------------------

def apply_decision_rule(random_rho, protein_rho, per_prot_std):
    if random_rho < 0.3:
        return "NULL"
    delta = random_rho - protein_rho
    if delta >= 0.10:
        return "LEAKY"
    if random_rho >= 0.5 and delta <= 0.05 and per_prot_std >= 0.15:
        return "HETEROGENEOUS"
    if random_rho >= 0.5 and delta <= 0.05 and per_prot_std <= 0.10:
        return "ROBUST"
    if 0.3 <= random_rho < 0.5:
        return "WEAK"
    return f"INTERMEDIATE (rho={random_rho:.3f}, delta={delta:.3f}, std={per_prot_std:.3f})"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    # ── 1. Load variants ──────────────────────────────────────────────────────
    variants = load_s1724_variants()
    print(f"Loaded {len(variants)} S1724 variants across "
          f"{len(set(v['protein'] for v in variants))} proteins")

    proteins = np.array([v["protein"] for v in variants])
    ddg      = np.array([v["ddg"]     for v in variants])

    # ── 2. Cluster assignment for family-split analogue ───────────────────────
    cluster_map = assign_protein_clusters(variants)
    n_clusters = len(set(cluster_map.values()))
    print(f"Protein clusters: {len(set(proteins))} proteins → {n_clusters} clusters")

    # ── 3. Embeddings ─────────────────────────────────────────────────────────
    if all(os.path.exists(p) for p in [WT_MEAN_EMB, MUT_MEAN_EMB, WT_POS_EMB, MUT_POS_EMB]):
        print("Loading cached embeddings...")
        wt_mean  = np.load(WT_MEAN_EMB)
        mut_mean = np.load(MUT_MEAN_EMB)
        wt_pos   = np.load(WT_POS_EMB)
        mut_pos  = np.load(MUT_POS_EMB)
    else:
        wt_mean, mut_mean, wt_pos, mut_pos = extract_embeddings(variants)

    delta_mean = mut_mean - wt_mean
    delta_pos  = mut_pos  - wt_pos

    print(f"Embeddings: delta_mean {delta_mean.shape}, delta_pos {delta_pos.shape}")

    # ── 4. Multi-seed CV ──────────────────────────────────────────────────────
    results_by_seed = []
    for seed in range(N_SEEDS):
        print(f"\n── Seed {seed} ──")

        splits_random  = random_split_cv(len(variants), N_FOLDS, seed)
        splits_protein = protein_split_cv(proteins, N_FOLDS, seed)
        splits_cluster = cluster_split_cv(proteins, cluster_map, N_FOLDS, seed)

        seed_result = {"seed": seed}

        for feat_name, X in [("delta_mean", delta_mean), ("delta_pos", delta_pos)]:
            for split_name, splits in [
                ("random",  splits_random),
                ("protein", splits_protein),
                ("cluster", splits_cluster),
            ]:
                key = f"{feat_name}_{split_name}"
                res = run_ridge_with_auroc(X, ddg, splits)
                seed_result[key] = res
                if res:
                    print(f"  {key}: ρ={res['spearman_mean']:.3f}±{res['spearman_std']:.3f}  "
                          f"AUROC={res['auroc_mean']:.3f}")

        results_by_seed.append(seed_result)

    # ── 5. Per-protein Spearman distribution ──────────────────────────────────
    print("\nPer-protein Spearman (leave-one-protein-out)...")
    per_prot = per_protein_spearman(delta_mean, ddg, proteins)
    prot_rhos = [v["spearman"] for v in per_prot.values()]
    per_prot_std = float(np.std(prot_rhos)) if prot_rhos else float("nan")
    per_prot_mean = float(np.mean(prot_rhos)) if prot_rhos else float("nan")

    print(f"  Per-protein ρ: mean={per_prot_mean:.3f}  std={per_prot_std:.3f}  "
          f"min={min(prot_rhos):.3f}  max={max(prot_rhos):.3f}  "
          f"n={len(prot_rhos)}")

    with open(os.path.join(OUT, "per_protein_spearman.json"), "w") as f:
        json.dump(per_prot, f, indent=2)

    # ── 6. Aggregate across seeds ─────────────────────────────────────────────
    summary = {}
    all_keys = set()
    for sr in results_by_seed:
        for k, v in sr.items():
            if isinstance(v, dict):
                all_keys.add(k)

    for key in sorted(all_keys):
        vals_rho   = [sr[key]["spearman_mean"] for sr in results_by_seed if key in sr and sr[key]]
        vals_auroc = [sr[key]["auroc_mean"]    for sr in results_by_seed if key in sr and sr[key]]
        if not vals_rho:
            continue
        summary[key] = {
            "spearman_mean": float(np.mean(vals_rho)),
            "spearman_std":  float(np.std(vals_rho)),
            "auroc_mean":    float(np.nanmean(vals_auroc)),
            "auroc_std":     float(np.nanstd(vals_auroc)),
        }

    summary["per_protein"] = {
        "spearman_mean": per_prot_mean,
        "spearman_std":  per_prot_std,
        "n_proteins":    len(prot_rhos),
    }

    # ── 7. Decision rule ──────────────────────────────────────────────────────
    dm_random  = summary.get("delta_mean_random",  {}).get("spearman_mean", float("nan"))
    dm_protein = summary.get("delta_mean_protein", {}).get("spearman_mean", float("nan"))

    verdict = apply_decision_rule(dm_random, dm_protein, per_prot_std)
    summary["verdict"] = verdict
    summary["n_variants"] = len(variants)
    summary["n_proteins"] = len(set(proteins))
    summary["n_clusters"] = n_clusters
    summary["n_seeds"] = N_SEEDS

    print(f"\n{'='*60}")
    print(f"VERDICT: {verdict}")
    print(f"  delta_mean random ρ  : {dm_random:.3f}")
    print(f"  delta_mean protein ρ : {dm_protein:.3f}")
    print(f"  Δ (random − protein) : {dm_random - dm_protein:.3f}")
    print(f"  per-protein ρ std    : {per_prot_std:.3f}")
    print(f"{'='*60}")

    with open(os.path.join(OUT, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nResults written to {OUT}/")


if __name__ == "__main__":
    main()
