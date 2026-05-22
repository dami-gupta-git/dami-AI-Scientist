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
        "val_accuracy": np.mean([[e["accuracy"] for e in vl] for vl in val_logs], axis=0),
    }

labels = {"run_0": "Evo2-7B (layer 28) — Glycolysis vs T-cell"}


def color_palette(n):
    cmap = plt.get_cmap("tab20")
    return [mcolors.rgb2hex(cmap(i)) for i in np.linspace(0, 1, n)]


runs = list(labels.keys())
colors = color_palette(len(runs))

# AUROC curve
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
plt.axhline(0.78, linestyle="--", color="gray", label="DNA repair vs TSG baseline (0.78)")
plt.axhline(0.5, linestyle=":", color="black", label="Chance (0.50)")
plt.title("Validation AUROC: Glycolysis vs T-cell Signaling")
plt.xlabel("Epoch")
plt.ylabel("AUROC")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("val_auroc.png")
plt.close()

# Final bar comparison
if final_results:
    run_names = [r for r in runs if r in final_results]
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for ax, metric, mlabel in zip(axes,
                                   ["final_val_accuracy", "final_val_auroc", "final_val_f1"],
                                   ["Accuracy", "AUROC", "F1"]):
        means = [final_results[r]["means"][metric] for r in run_names]
        errs = [final_results[r]["stderrs"][metric] for r in run_names]
        ax.bar([labels[r] for r in run_names], means, yerr=errs, color=colors[:len(run_names)], capsize=5)
        if metric == "final_val_auroc":
            ax.axhline(0.78, linestyle="--", color="gray", label="Prior baseline")
        ax.set_ylim(0, 1)
        ax.set_title(mlabel)
        ax.set_ylabel(mlabel)
        ax.tick_params(axis="x", rotation=15)
    plt.suptitle("Validation: Glycolysis vs T-cell Signaling")
    plt.tight_layout()
    plt.savefig("final_metrics_bar.png")
    plt.close()
