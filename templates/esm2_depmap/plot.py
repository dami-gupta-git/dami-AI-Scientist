import json
import os
import os.path as osp

import matplotlib.pyplot as plt
import numpy as np

folders = [f for f in os.listdir("./") if f.startswith("run") and osp.isdir(f)]

for folder in folders:
    fp = osp.join(folder, "final_info.json")
    if not osp.exists(fp):
        continue
    with open(fp) as f:
        results = json.load(f)
    means = results.get("means", results)

    # Plot 1: Null distribution vs observed Mantel r
    null_path = osp.join(folder, "null_dist.npy")
    if osp.exists(null_path):
        null = np.load(null_path)
        r_obs = means.get("mantel_r", 0)
        p_val = means.get("mantel_p_value", 1)
        plt.figure(figsize=(8, 5))
        plt.hist(null, bins=50, color="steelblue", alpha=0.7, label="Null distribution")
        plt.axvline(r_obs, color="red", linewidth=2, label=f"Observed r = {r_obs:.4f} (p={p_val:.4f})")
        plt.xlabel("Mantel r (Spearman)")
        plt.ylabel("Count")
        plt.title("Mantel Test: ESM-2 Distance vs Essentiality Distance")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(f"{folder}_mantel_null.png")
        plt.close()

    # Plot 2: Scatter sample of pairwise distances
    esm2_path = osp.join(folder, "esm2_dist.npy")
    ess_path = osp.join(folder, "ess_dist.npy")
    if osp.exists(esm2_path) and osp.exists(ess_path):
        esm2_dist = np.load(esm2_path)
        ess_dist = np.load(ess_path)
        n = esm2_dist.shape[0]
        idx = np.triu_indices(n, k=1)
        esm2_flat = esm2_dist[idx]
        ess_flat = ess_dist[idx]

        # Sample 10k pairs for scatter
        rng = np.random.RandomState(0)
        sample = rng.choice(len(esm2_flat), size=min(10000, len(esm2_flat)), replace=False)

        plt.figure(figsize=(7, 6))
        plt.scatter(esm2_flat[sample], 1 - ess_flat[sample], alpha=0.1, s=1, color="steelblue")
        plt.xlabel("ESM-2 Cosine Distance")
        plt.ylabel("Essentiality Correlation (Pearson)")
        plt.title(f"ESM-2 Distance vs Essentiality Correlation\n(Mantel r={means.get('mantel_r', 0):.4f}, p={means.get('mantel_p_value', 1):.4f})")
        plt.grid(True, alpha=0.2)
        plt.tight_layout()
        plt.savefig(f"{folder}_scatter.png")
        plt.close()

    print(f"{folder}: Mantel r={means.get('mantel_r', 'N/A')}, p={means.get('mantel_p_value', 'N/A')}, n_genes={means.get('n_genes', 'N/A')}")
