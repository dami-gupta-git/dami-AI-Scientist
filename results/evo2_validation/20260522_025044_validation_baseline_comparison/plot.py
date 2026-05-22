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

labels = {
    "run_0": "Baseline (layer 28, seed 0)",
    "run_1": "Layer 28 (seed 0) + Centroid Analysis",
    "run_2": "Layer 28 (5 seeds)",
    "run_3": "Layer 14 (seed 0)",
    "run_4": "Layer 40 (seed 0)",
}


def color_palette(n):
    cmap = plt.get_cmap("tab20")
    return [mcolors.rgb2hex(cmap(i)) for i in np.linspace(0, 1, n)]


runs = list(labels.keys())
colors = color_palette(len(runs))

# Extract centroid similarity data if available
centroid_data = {}
for run in runs:
    if run in final_results and "centroid_cosine_similarity" in final_results[run]["means"]:
        centroid_data[run] = {
            "similarity": final_results[run]["means"]["centroid_cosine_similarity"],
            "distance": final_results[run]["means"]["centroid_cosine_distance"],
            "within_glycolysis": final_results[run]["means"]["within_glycolysis_cosine_sim"],
            "within_tcell": final_results[run]["means"]["within_tcell_cosine_sim"],
        }

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

# Final bar comparison - Best AUROC
if final_results:
    run_names = [r for r in runs if r in final_results]
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Best validation AUROC
    ax = axes[0]
    means = [final_results[r]["means"]["best_val_auroc"] for r in run_names]
    errs = [final_results[r]["stderrs"]["best_val_auroc"] for r in run_names]
    bars = ax.bar(range(len(run_names)), means, yerr=errs, color=colors[:len(run_names)], capsize=5)
    ax.axhline(0.78, linestyle="--", color="red", linewidth=2, label="DNA repair vs TSG (0.78)")
    ax.axhline(0.5, linestyle=":", color="black", linewidth=1, label="Chance (0.50)")
    ax.set_ylim(0, 1)
    ax.set_title("Best Validation AUROC", fontsize=12, fontweight='bold')
    ax.set_ylabel("AUROC", fontsize=11)
    ax.set_xticks(range(len(run_names)))
    ax.set_xticklabels([labels[r] for r in run_names], rotation=45, ha='right', fontsize=9)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis='y')
    
    # Final validation AUROC
    ax = axes[1]
    means = [final_results[r]["means"]["final_val_auroc"] for r in run_names]
    errs = [final_results[r]["stderrs"]["final_val_auroc"] for r in run_names]
    bars = ax.bar(range(len(run_names)), means, yerr=errs, color=colors[:len(run_names)], capsize=5)
    ax.axhline(0.5, linestyle=":", color="black", linewidth=1, label="Chance (0.50)")
    ax.set_ylim(0, 1)
    ax.set_title("Final Validation AUROC", fontsize=12, fontweight='bold')
    ax.set_ylabel("AUROC", fontsize=11)
    ax.set_xticks(range(len(run_names)))
    ax.set_xticklabels([labels[r] for r in run_names], rotation=45, ha='right', fontsize=9)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis='y')
    
    # Final validation accuracy
    ax = axes[2]
    means = [final_results[r]["means"]["final_val_accuracy"] for r in run_names]
    errs = [final_results[r]["stderrs"]["final_val_accuracy"] for r in run_names]
    bars = ax.bar(range(len(run_names)), means, yerr=errs, color=colors[:len(run_names)], capsize=5)
    ax.set_ylim(0, 1)
    ax.set_title("Final Validation Accuracy", fontsize=12, fontweight='bold')
    ax.set_ylabel("Accuracy", fontsize=11)
    ax.set_xticks(range(len(run_names)))
    ax.set_xticklabels([labels[r] for r in run_names], rotation=45, ha='right', fontsize=9)
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.suptitle("Glycolysis vs T-cell Signaling: Performance Across Experiments", fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig("final_metrics_bar.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: final_metrics_bar.png")

# Centroid similarity analysis plot
if centroid_data:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Between-class vs within-class similarity
    ax = axes[0]
    run_names = list(centroid_data.keys())
    x = np.arange(len(run_names))
    width = 0.25
    
    between_class = [centroid_data[r]["similarity"] for r in run_names]
    within_glycolysis = [centroid_data[r]["within_glycolysis"] for r in run_names]
    within_tcell = [centroid_data[r]["within_tcell"] for r in run_names]
    
    ax.bar(x - width, between_class, width, label='Between-class (centroid)', color='red', alpha=0.7)
    ax.bar(x, within_glycolysis, width, label='Within glycolysis', color='blue', alpha=0.7)
    ax.bar(x + width, within_tcell, width, label='Within T-cell', color='green', alpha=0.7)
    
    ax.set_ylabel('Cosine Similarity', fontsize=11)
    ax.set_title('Embedding Geometry: Between vs Within-Class Similarity', fontsize=12, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels([labels[r] for r in run_names], rotation=45, ha='right', fontsize=9)
    ax.legend(fontsize=9)
    ax.set_ylim(0.94, 1.0)
    ax.grid(True, alpha=0.3, axis='y')
    ax.axhline(1.0, linestyle='--', color='gray', linewidth=0.5)
    
    # Centroid distance
    ax = axes[1]
    distances = [centroid_data[r]["distance"] for r in run_names]
    bars = ax.bar(range(len(run_names)), distances, color=colors[:len(run_names)])
    ax.set_ylabel('Cosine Distance', fontsize=11)
    ax.set_title('Centroid Cosine Distance (1 - similarity)', fontsize=12, fontweight='bold')
    ax.set_xticks(range(len(run_names)))
    ax.set_xticklabels([labels[r] for r in run_names], rotation=45, ha='right', fontsize=9)
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add text annotations
    for i, (bar, dist) in enumerate(zip(bars, distances)):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.0001,
                f'{dist:.4f}', ha='center', va='bottom', fontsize=8)
    
    plt.suptitle('Unsupervised Embedding Analysis: Lack of Geometric Separation', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig("centroid_analysis.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: centroid_analysis.png")

# Layer comparison plot (runs 1, 3, 4)
layer_runs = {"run_1": ("Layer 28", 28), "run_3": ("Layer 14", 14), "run_4": ("Layer 40", 40)}
layer_runs_available = {k: v for k, v in layer_runs.items() if k in final_results}

if len(layer_runs_available) >= 2:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    run_names = list(layer_runs_available.keys())
    layer_names = [layer_runs_available[r][0] for r in run_names]
    layer_nums = [layer_runs_available[r][1] for r in run_names]
    
    # Best AUROC by layer
    ax = axes[0]
    aurocs = [final_results[r]["means"]["best_val_auroc"] for r in run_names]
    ax.plot(layer_nums, aurocs, 'o-', linewidth=2, markersize=10, color='steelblue')
    ax.axhline(0.78, linestyle="--", color="red", linewidth=2, label="DNA repair vs TSG (0.78)")
    ax.axhline(0.5, linestyle=":", color="black", linewidth=1, label="Chance (0.50)")
    ax.set_xlabel('Evo2 Layer', fontsize=11)
    ax.set_ylabel('Best Validation AUROC', fontsize=11)
    ax.set_title('Performance Across Evo2 Layers', fontsize=12, fontweight='bold')
    ax.set_xticks(layer_nums)
    ax.set_xticklabels(layer_names)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0.4, 0.85)
    
    # Centroid distance by layer
    ax = axes[1]
    if all(r in centroid_data for r in run_names):
        distances = [centroid_data[r]["distance"] for r in run_names]
        ax.plot(layer_nums, distances, 'o-', linewidth=2, markersize=10, color='coral')
        ax.set_xlabel('Evo2 Layer', fontsize=11)
        ax.set_ylabel('Centroid Cosine Distance', fontsize=11)
        ax.set_title('Geometric Separation Across Layers', fontsize=12, fontweight='bold')
        ax.set_xticks(layer_nums)
        ax.set_xticklabels(layer_names)
        ax.grid(True, alpha=0.3)
        
        # Add value annotations
        for layer, dist in zip(layer_nums, distances):
            ax.text(layer, dist + 0.00005, f'{dist:.4f}', ha='center', va='bottom', fontsize=9)
    
    plt.suptitle('Layer Analysis: Consistent Lack of Functional Encoding', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig("layer_comparison.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: layer_comparison.png")

print("\nAll plots generated successfully!")
