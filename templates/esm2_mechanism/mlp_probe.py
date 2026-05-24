"""
MLP probe: gene-split and family-split CV on delta_mean and wt_only.
Answers: does the nonlinear mechanism signal in delta survive family-split?

Usage:
    python mlp_probe.py --run_dir run_0 --model esm2_t33_650M_UR50D
"""
import argparse, json, os
from collections import Counter
import numpy as np
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import roc_auc_score, f1_score
from sklearn.preprocessing import StandardScaler

parser = argparse.ArgumentParser()
parser.add_argument("--run_dir", default="run_0")
parser.add_argument("--model", default="esm2_t33_650M_UR50D")
parser.add_argument("--seed", type=int, default=42)
args = parser.parse_args()

data_dir = os.path.join(args.run_dir, "data")

wt     = np.load(f"{data_dir}/embeddings_wt_{args.model}.npy").astype(np.float64)
mut    = np.load(f"{data_dir}/embeddings_mut_{args.model}.npy").astype(np.float64)
deltas = mut - wt

with open(f"{data_dir}/gerasimavicius_variants.json") as f: variants = json.load(f)
with open(f"{data_dir}/pfam_families.json") as f: pfam_map = json.load(f)
with open(f"{data_dir}/sequences.json") as f: seq_cache = json.load(f)

valid = [v for v in variants if v["uniprot_id"] in seq_cache]
for v in valid: v["label_3class"] = "LOF" if v["mechanism"] in ("HI","AR") else v["mechanism"]
valid = valid[:deltas.shape[0]]
labels = np.array([v["label_3class"] for v in valid])
genes  = np.array([v["gene"] for v in valid])
n = len(labels)
print(f"N={n}, {dict(Counter(labels))}")


def gene_splits(genes, n_folds=5, seed=42):
    u = np.array(sorted(set(genes)))
    np.random.RandomState(seed).shuffle(u)
    return [(np.where(~np.isin(genes, f))[0], np.where(np.isin(genes, f))[0])
            for f in np.array_split(u, n_folds)]


def family_splits(genes, pfam_map, n_folds=5, seed=42):
    g2p = {g: pfam_map[g] for g in np.unique(genes) if pfam_map.get(g)}
    fams = np.array(sorted(set(g2p.values())))
    np.random.RandomState(seed).shuffle(fams)
    splits = []
    for fold_fams in np.array_split(fams, n_folds):
        fs = set(fold_fams)
        te = np.array([genes[i] in g2p and g2p[genes[i]] in fs for i in range(n)])
        tr = np.array([genes[i] in g2p and g2p[genes[i]] not in fs for i in range(n)])
        if tr.sum() >= 10 and te.sum() >= 5:
            splits.append((np.where(tr)[0], np.where(te)[0]))
    print(f"  family splits={len(splits)}, families={len(fams)}")
    return splits


def run_mlp(X, y, splits, tag):
    f1s, gof, dn, lof = [], [], [], []
    for tr, te in splits:
        sc = StandardScaler()
        Xtr = sc.fit_transform(X[tr]); Xte = sc.transform(X[te])
        ytr, yte = y[tr], y[te]
        if len(set(ytr)) < 2: continue
        clf = MLPClassifier(hidden_layer_sizes=(256,), max_iter=300, random_state=args.seed)
        clf.fit(Xtr, ytr)
        pred = clf.predict(Xte); proba = clf.predict_proba(Xte)
        f1s.append(f1_score(yte, pred, average="macro", zero_division=0))
        for i, cls in enumerate(clf.classes_):
            yb = (yte == cls).astype(int)
            if yb.sum() > 0 and (1-yb).sum() > 0:
                a = roc_auc_score(yb, proba[:, i])
                if cls == "GOF": gof.append(a)
                elif cls == "DN": dn.append(a)
                elif cls == "LOF": lof.append(a)
    result = {
        "macro_f1": float(np.mean(f1s)),
        "auroc_GOF": float(np.mean(gof)) if gof else float("nan"),
        "auroc_DN":  float(np.mean(dn))  if dn  else float("nan"),
        "auroc_LOF": float(np.mean(lof)) if lof else float("nan"),
        "n_folds": len(f1s),
    }
    print(f"{tag}: F1={result['macro_f1']:.3f}  GOF={result['auroc_GOF']:.3f}  "
          f"DN={result['auroc_DN']:.3f}  LOF={result['auroc_LOF']:.3f}  (folds={len(f1s)})")
    return result


gs = gene_splits(genes)
fs = family_splits(genes, pfam_map)

results = {}
print()
results["delta_mean_gene_split"]   = run_mlp(deltas, labels, gs, "delta_mean   gene-split  ")
results["delta_mean_family_split"] = run_mlp(deltas, labels, fs, "delta_mean   family-split")
print()
results["wt_only_gene_split"]      = run_mlp(wt, labels, gs,     "wt_only      gene-split  ")
results["wt_only_family_split"]    = run_mlp(wt, labels, fs,     "wt_only      family-split")

print("\n=== HEADLINE ===")
dg = results["delta_mean_gene_split"]["macro_f1"]
df = results["delta_mean_family_split"]["macro_f1"]
print(f"delta_mean:  gene-split F1={dg:.3f} -> family-split F1={df:.3f}  (delta={dg-df:+.3f})")
wg = results["wt_only_gene_split"]["macro_f1"]
wf = results["wt_only_family_split"]["macro_f1"]
print(f"wt_only:     gene-split F1={wg:.3f} -> family-split F1={wf:.3f}  (delta={wg-wf:+.3f})")
print(f"Chance: 0.333")

out_path = os.path.join(args.run_dir, "mlp_probe_results.json")
with open(out_path, "w") as f:
    json.dump(results, f, indent=2)
print(f"\nResults written to {out_path}")
