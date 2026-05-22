import json
import os
import os.path as osp

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np

folders = [f for f in os.listdir("./") if f.startswith("run") and osp.isdir(f)]
final_results = {}
results_info = {}

for folder in folders:
    fp = osp.join(folder, "final_info.json")
    if not osp.exists(fp):
        continue
    with open(fp) as f:
        final_results[folder] = json.load(f)

    rp = osp.join(folder, "all_results.npy")
    if not osp.exists(rp):
        continue
    results_dict = np.load(rp, allow_pickle=True).item()
    seeds = sorted({k.split("_")[0] for k in results_dict if "train_log" in k})
    train_logs = [results_dict[f"{s}_train_log"] for s in seeds]
    val_logs = [results_dict[f"{s}_val_log"] for s in seeds]
    epochs = [e["epoch"] for e in train_logs[0]]
    results_info[folder] = {
        "epochs": epochs,
        "val_auroc": np.mean([[e["auroc"] for e in vl] for vl in val_logs], axis=0),
        "val_auroc_sterr": np.std([[e["auroc"] for e in vl] for vl in val_logs], axis=0) / max(np.sqrt(len(seeds)), 1),
    }

labels = {"run_0": "ESM-2 650M (protein embeddings)"}


def color_palette(n):
    cmap = plt.get_cmap("tab20")
    return [mcolors.rgb2hex(cmap(i)) for i in np.linspace(0, 1, n)]


runs = list(labels.keys())
colors = color_palette(len(runs))

# AUROC curve with Evo2 reference line
plt.figure(figsize=(8, 5))
for i, run in enumerate(runs):
    if run not in results_info:
        continue
    info = results_info[run]
    plt.plot(info["epochs"], info["val_auroc"], label=labels[run], color=colors[i])
    plt.fill_between(info["epochs"],
                     info["val_auroc"] - info["val_auroc_sterr"],
                     info["val_auroc"] + info["val_auroc_sterr"],
                     color=colors[i], alpha=0.2)
plt.axhline(0.78, linestyle="--", color="orange", label="Evo2-7B DNA baseline (0.78)")
plt.axhline(0.5, linestyle=":", color="black", label="Chance (0.50)")
plt.title("ESM-2 vs Evo2: Validation AUROC on DNA Repair vs TSG")
plt.xlabel("Epoch")
plt.ylabel("AUROC")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("val_auroc_comparison.png")
plt.close()

# Comparison bar chart: ESM-2 vs Evo2 on key metrics
if final_results:
    evo2_reference = {
        "best_val_auroc": 0.78,
        "silhouette_score": 0.044,
        "distance_ratio": 0.963,
        "kmeans_ari": -0.001,
    }
    run_names = [r for r in runs if r in final_results]
    metrics = ["best_val_auroc", "silhouette_score", "distance_ratio"]
    metric_labels = ["Best AUROC", "Silhouette Score", "Distance Ratio (within/between)"]

    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    for ax, metric, mlabel in zip(axes, metrics, metric_labels):
        esm2_vals = [final_results[r]["means"].get(metric, 0) for r in run_names]
        esm2_errs = [final_results[r]["stderrs"].get(metric, 0) for r in run_names]
        evo2_val = evo2_reference[metric]

        x = np.arange(len(run_names) + 1)
        bar_labels = [labels[r] for r in run_names] + ["Evo2-7B DNA (reference)"]
        bar_vals = esm2_vals + [evo2_val]
        bar_errs = esm2_errs + [0]
        bar_colors = colors[:len(run_names)] + ["#ff7f0e"]

        ax.bar(x, bar_vals, yerr=bar_errs, color=bar_colors, capsize=5)
        ax.set_xticks(x)
        ax.set_xticklabels(bar_labels, rotation=15, ha="right", fontsize=8)
        ax.set_title(mlabel)
        ax.set_ylabel(mlabel)
        ax.grid(True, alpha=0.2, axis="y")

    plt.suptitle("ESM-2 Protein vs Evo2 DNA: DNA Repair vs Tumor Suppressor Classification")
    plt.tight_layout()
    plt.savefig("esm2_vs_evo2_comparison.png")
    plt.close()
