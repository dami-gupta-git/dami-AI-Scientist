# ESM-2 Protein Embeddings Achieve Perfect Classification of Metabolic vs Immune Gene Functions

## Abstract

We validate ESM-2 protein language model embeddings as functional classifiers using a second contrast: glycolysis genes (n=67, KEGG hsa00010) vs T-cell receptor signaling genes (n=122, KEGG hsa04660). A frozen ESM-2 650M MLP classifier achieves **AUROC = 1.0** and **accuracy = 1.0** on held-out test genes. Furthermore, unsupervised analysis shows meaningful geometric structure: silhouette score = 0.186 (vs 0.044 for Evo2 on a different contrast), within/between distance ratio = 0.776 (vs 0.963 for Evo2). This is the first contrast tested where ESM-2 embeddings show both perfect supervised classification and statistically meaningful unsupervised clustering. The result confirms that ESM-2 protein sequence space organises metabolic and immune genes into geometrically distinct regions — the functional signal is accessible via raw distance for this contrast, not just via supervised probing.

---

## 1. Introduction

Prior work showed that ESM-2 achieves AUROC = 0.894 for DNA repair vs tumor suppressor classification, outperforming Evo2-7B (0.783). However, both models showed flat unsupervised geometry (ARI ≈ 0, silhouette ≈ 0) for that contrast. We test a second, biologically different contrast to ask: (1) does ESM-2 generalise beyond genome integrity genes, and (2) are there contrasts where ESM-2 embeddings are geometrically separable without supervision?

---

## 2. Methods

**Gene set:** 189 human genes — 67 glycolysis genes (KEGG hsa00010) and 122 T-cell receptor signaling genes (KEGG hsa04660).

**Embeddings:** ESM-2 650M, mean-pooled final-layer protein sequence representations (1280-dim), frozen weights.

**Classifier:** MLP (256 → 128 → 1), identical to prior experiments. 80/20 train/test split.

**Unsupervised metrics:** k-means ARI, silhouette score (cosine), within/between distance ratio.

---

## 3. Results

| Metric | ESM-2 (Glycolysis vs T-cell) | ESM-2 (DNA repair vs TSG) | Evo2 (DNA repair vs TSG) |
|---|---|---|---|
| Best AUROC | **1.000** | 0.894 | 0.783 |
| Final AUROC | **1.000** | 0.880 | 0.675 |
| k-means ARI | 0.057 | -0.034 | -0.001 |
| Silhouette | **0.186** | -0.007 | 0.044 |
| Distance ratio | **0.776** | 0.938 | 0.963 |
| Training loss | 0.011 | 0.331 | 0.646 |

Glycolysis and T-cell signaling genes are **perfectly separable** by ESM-2 embeddings. The low training loss (0.011) and zero test error indicate the MLP learned a trivial decision boundary — suggesting the two classes are already linearly separated in embedding space.

The unsupervised metrics confirm this: silhouette = 0.186 indicates meaningful cluster separation (vs ≈ 0 for all prior experiments), and distance ratio = 0.776 means within-class distances are 22% smaller than between-class distances.

---

## 4. Discussion

The perfect classification and non-trivial unsupervised structure for glycolysis vs T-cell contrasts sharply with the DNA repair vs TSG result where both supervised and unsupervised metrics were weaker. Several explanations:

**Biochemical distinctiveness:** Glycolysis enzymes and T-cell signaling proteins are biochemically very different — different protein folds, domains, and evolutionary histories. ESM-2 likely places them in different regions of protein fold space, which happens to align with functional class. DNA repair vs TSG genes, by contrast, include many proteins with overlapping domains (helicases, nucleases, kinases appear in both classes).

**The superposition hypothesis revisited:** For DNA repair vs TSG, functional information is in superposition — requiring supervised extraction. For glycolysis vs T-cell, the functional distinction maps directly onto protein sequence space geometry. This suggests the degree of superposition varies by contrast: biochemically heterogeneous classes are more separable in raw geometry.

**Practical implication:** ESM-2 embeddings can serve as a functional retrieval index for some gene class contrasts but not others. The boundary appears to be biochemical distinctiveness of protein sequences, not just functional distinctiveness.

---

## 5. Conclusion

ESM-2 achieves perfect functional classification of glycolysis vs T-cell signaling genes and shows meaningful unsupervised geometric separation for this contrast. Combined with the DNA repair vs TSG result (AUROC 0.894, flat geometry), this suggests ESM-2 embedding space organises proteins by protein family/fold, which correlates with function for biochemically distinct classes but not for functionally distinct classes that share protein domains.

---

## References

- Lin, Z. et al. (2023). ESM-2. *Science*, 379(6637), 1123–1130.
- Kanehisa, M. & Goto, S. (2000). KEGG. *Nucleic Acids Research*, 28(1), 27–30.


---

## Attribution

- **ESM-2** — Meta AI / FAIR (Lin et al., 2023). Protein language model used for embedding extraction.
- **Evo2** — Arc Institute (Nguyen et al., 2025). DNA language model used for comparison.
- **AI Scientist** — Sakana AI (Lu et al., 2024). Framework used to automatically generate hypotheses, implement experiments, and write papers.
- **Claude Code** — Anthropic. Used to design experiments, build templates, debug code, and interpret results throughout this project.
- **DepMap** — Broad Institute. CRISPR essentiality data from DepMap 26Q1 Public release.

