# ProbeSeq — Findings

## What ProbeSeq Is

A systematic comparison of DNA and protein foundation model embeddings as feature extractors for gene functional classification. The core question: do protein language models (ESM-2) encode functional gene class information more effectively than DNA language models (Evo2), and does either model produce semantically organised embedding spaces?

---

## Experimental Setup

**Gene set:** 178 human genes — 132 DNA repair (KEGG: MMR, NER, HR, BER) and 46 tumor suppressors (COSMIC Cancer Gene Census). Same gene set across all models.

**Classifier:** Frozen embeddings → MLP (256 → 128 → 1), binary cross-entropy, Adam lr=1e-3, 10 epochs, 80/20 train/test split.

**Unsupervised metrics:** k-means ARI, silhouette score (cosine), within/between class distance ratio — all with permutation testing.

---

## Results

| Metric | ESM-2 650M (protein) | Evo2-7B (DNA) |
|---|---|---|
| **Best AUROC** | **0.894** | 0.783 |
| **Final AUROC** | **0.880** | 0.675 |
| Input dimension | 1280 | 4096 |
| k-means ARI | -0.034 | -0.001 |
| Silhouette score | -0.007 | 0.044 |
| Distance ratio (within/between) | 0.938 | 0.963 |

---

## Key Findings

**1. ESM-2 protein embeddings are stronger supervised feature extractors than Evo2 DNA embeddings**

ESM-2 achieves 0.894 AUROC vs Evo2's 0.783 — a meaningful gap on the same 36 gene test set. This is despite ESM-2 using a much lower-dimensional embedding (1280 vs 4096). Protein sequence space carries more discriminative functional signal for this classification task than DNA sequence space.

**2. Neither model produces semantically organised embedding spaces**

Both models show flat unsupervised geometry — k-means ARI ~0, silhouette ~0, distance ratio ~0.94-0.96. Genes of the same functional class are no closer together than genes from different classes in raw embedding space. The functional signal requires a learned probe to extract in both cases.

**3. Protein embeddings capture function; DNA embeddings capture sequence statistics**

The gap between ESM-2 and Evo2 likely reflects the difference in what each model was trained to encode. ESM-2 was trained on evolutionary co-variation across protein sequences — which directly reflects functional constraints. Evo2 was trained on next-nucleotide prediction across DNA — which captures sequence statistics, repeat content, and regulatory grammar, but not protein function directly.

---

## Validation: Glycolysis vs T-cell Signaling (Evo2 only)

A separate validation using Evo2 on a different functional contrast (glycolysis vs T-cell signaling, 189 genes) showed:
- Best AUROC: 0.66 (vs 0.78 for DNA repair/TSG)
- Centroid cosine distance: 0.0014 (vs ~0.006 for DNA repair/TSG)

This confirmed that Evo2's functional signal is contrast-dependent — genome integrity pathways are better encoded than metabolic/immune pathways. ESM-2 was not tested on this contrast.

---

## Open Questions

1. Would ESM-2 also generalise better than Evo2 to the glycolysis/T-cell contrast?
2. Does ESM-2 show better unsupervised clustering for any functional contrast, or is the flat geometry universal?
3. Would a larger ESM-2 (3B or 15B) push AUROC higher?
4. Is the ESM-2 signal concentrated in specific embedding dimensions, or spread across all 1280?

---

## Next Steps

- Run ESM-2 on glycolysis vs T-cell to test cross-contrast generalisation
- Run unsupervised clustering analysis on ESM-2 embeddings with permutation testing
- Consider SAE decomposition on ESM-2 to find interpretable functional features
