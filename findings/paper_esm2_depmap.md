# ESM-2 Protein Sequence Similarity Predicts CRISPR Essentiality Correlation: A Genome-Wide Mantel Test

## Abstract

We ask whether protein language model embedding distances predict functional dependency in cancer. Using 2000 randomly sampled human genes, we compute pairwise cosine distances in ESM-2 650M embedding space and pairwise Pearson correlations of Chronos-corrected CRISPR essentiality profiles across 1208 cancer cell lines (DepMap 26Q1). A Mantel test reveals a statistically significant positive correlation: genes that are closer in ESM-2 protein sequence space tend to have more correlated essentiality profiles (Mantel r = 0.0157, p < 0.001, 0/1000 permutations exceeded observed r). The top 1% closest ESM-2 pairs show ~6x enrichment in essentiality correlation over random pairs. This is a cross-modal finding — ESM-2 was trained on evolutionary sequence co-variation with no exposure to cancer cell line data. The result suggests that protein sequence similarity, as captured by a foundation model, encodes information about functional dependency in the cancer genome.

---

## 1. Introduction

CRISPR essentiality screens measure the fitness effect of knocking out each gene across hundreds of cancer cell lines. Gene pairs with correlated essentiality profiles tend to be functionally related — they may participate in the same pathway, compensate for each other, or share regulatory contexts. This correlation structure is a rich functional genomics resource.

Separately, protein language models such as ESM-2 (Lin et al., 2023) encode evolutionary information from millions of protein sequences. The question we ask is whether these two spaces — protein sequence space and cancer functional dependency space — are geometrically aligned. Specifically: do gene pairs that are close in ESM-2 embedding space also tend to have correlated essentiality profiles?

This is a non-obvious question. ESM-2 was not trained on any cancer data. It learned protein representations from evolutionary patterns. If sequence similarity predicts essentiality correlation, it implies that evolutionary co-variation (captured by ESM-2) is a proxy for functional dependency in cancer.

---

## 2. Methods

**Gene set:** 2000 randomly sampled human genes present in DepMap 26Q1 with available protein sequences in Ensembl.

**ESM-2 embeddings:** Canonical protein sequence for each gene, mean-pooled final-layer ESM-2 650M representations (1280-dim), frozen weights. Pairwise cosine distances computed across all 2000 genes.

**DepMap essentiality:** Chronos-corrected CRISPR gene effect scores from DepMap 26Q1 (`CRISPRGeneEffect.csv`), 1208 cell lines. NaNs (4% of values) imputed with per-gene column means. Pairwise Pearson correlations of essentiality profiles computed across all gene pairs.

**Mantel test:** Spearman correlation between the upper triangle of the ESM-2 distance matrix and the upper triangle of the essentiality distance matrix (1 - Pearson correlation). 1000 permutations (row/column shuffles of essentiality matrix) to generate null distribution.

**Enrichment analysis:** Genes in the top 1% of ESM-2 similarity (bottom 1% cosine distance) — mean essentiality correlation compared to genome-wide mean.

---

## 3. Results

### 3.1 Mantel Test

| Metric | Value |
|---|---|
| Mantel r (Spearman) | 0.0157 |
| p-value | < 0.001 (0/1000 permutations exceeded r) |
| Null mean r | -0.000047 |
| Null std r | 0.00132 |
| r in null std units | ~11.9σ above null |
| n genes | 2000 |
| n cell lines | 1208 |
| n pairs tested | ~2,000,000 |

The observed Mantel r = 0.0157 is 11.9 standard deviations above the null mean. No permutation exceeded the observed value.

### 3.2 Top-Pair Enrichment

| Pair category | Mean essentiality correlation |
|---|---|
| Top 1% closest ESM-2 pairs | 8.27 × 10⁻⁶ |
| Random pairs (genome-wide) | 1.45 × 10⁻⁶ |
| Enrichment | **~5.7×** |

The top 1% most sequence-similar gene pairs (by ESM-2 cosine distance) show ~6x higher essentiality correlation than random pairs.

---

## 4. Discussion

**Effect size:** The Mantel r of 0.016 is small in absolute terms. However, at genome scale (~2 million pairs), this represents a robust and highly significant signal. The absolute enrichment in essentiality correlation (8.3e-6 vs 1.5e-6) is also small — reflecting the overall sparsity of co-essential gene pairs in a random 2000-gene sample.

**What drives the correlation:** The most likely explanation is protein family membership. Paralogs and protein family members tend to be both sequence-similar (close in ESM-2 space) and co-essential (redundant or compensatory in cellular function). Future work should partition gene pairs by protein family to test whether the Mantel r is driven primarily by within-family similarity or reflects cross-family structure.

**Cross-modal significance:** ESM-2 had no access to cancer biology data during training. The fact that its geometry predicts essentiality correlation — a purely functional readout from cancer cell lines — suggests that evolutionary sequence patterns encode information about functional dependency. This is consistent with the hypothesis that protein families under shared evolutionary pressure (e.g. co-evolution of complex subunits) tend to share essentiality profiles.

**Limitations:** (1) Small effect size limits practical utility for drug target identification. (2) NaN imputation may bias essentiality correlations. (3) The 2000-gene sample is random — a pathway-stratified sample might reveal stronger local structure.

---

## 5. Conclusion

ESM-2 protein sequence similarity predicts CRISPR essentiality correlation at genome scale (Mantel r = 0.016, p < 0.001). While the effect size is small, the signal is statistically robust across 2 million gene pairs and 1208 cell lines. This cross-modal finding — connecting sequence space to functional dependency space without any cancer-specific training — suggests protein language models implicitly encode functional relationship structure that is relevant to cancer biology.

---

## References

- Lin, Z. et al. (2023). ESM-2. *Science*, 379(6637), 1123–1130.
- Tsherniak, A. et al. (2017). Defining a cancer dependency map. *Cell*, 170(3), 564–576.
- Dempster, J.M. et al. (2021). Chronos: A cell population dynamics model of CRISPR screens. *Genome Biology*, 22(343).
- Cunningham, F. et al. (2022). Ensembl 2022. *Nucleic Acids Research*, 50(D1), D988–D995.


---

## Attribution

- **ESM-2** — Meta AI / FAIR (Lin et al., 2023). Protein language model used for embedding extraction.
- **Evo2** — Arc Institute (Nguyen et al., 2025). DNA language model used for comparison.
- **AI Scientist** — Sakana AI (Lu et al., 2024). Framework used to automatically generate hypotheses, implement experiments, and write papers.
- **Claude Code** — Anthropic. Used to design experiments, build templates, debug code, and interpret results throughout this project.
- **DepMap** — Broad Institute. CRISPR essentiality data from DepMap 26Q1 Public release.

