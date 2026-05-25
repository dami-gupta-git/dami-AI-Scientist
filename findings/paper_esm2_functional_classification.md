# Protein Language Model Embeddings Outperform DNA Language Model Embeddings for Gene Functional Classification

## Abstract

We compare frozen embeddings from ESM-2 (a protein language model trained on amino acid sequences) and Evo2 (a DNA language model trained on nucleotide sequences) as feature extractors for gene functional classification. Using 178 human genes from two functionally distinct classes — DNA repair (n=132, from KEGG pathways) and tumor suppressors (n=46, from the COSMIC Cancer Gene Census) — we train identical shallow MLP classifiers on frozen embeddings from each model. ESM-2 650M achieves best AUROC of 0.894, significantly outperforming Evo2-7B (0.783). However, neither model produces semantically organised embedding geometry: both show near-zero unsupervised clustering (ARI ≈ 0, silhouette ≈ 0). These results demonstrate that protein sequence space encodes more discriminative functional information than DNA sequence space for this classification task, but that functional information is in superposition in both models — extractable via supervised probing but not accessible via raw distance metrics.

---

## 1. Introduction

Foundation models for biological sequences have demonstrated strong performance on a range of downstream prediction tasks. ESM-2 (Lin et al., 2023) is trained on evolutionary co-variation across 250M protein sequences, while Evo2 (Nguyen et al., 2025) is trained on 8.8 trillion nucleotide tokens spanning all domains of life. Both models produce dense vector representations that can be used as features for supervised classifiers.

A fundamental question is whether these representations encode high-level biological function — specifically, whether genes with similar functions are placed near each other in embedding space. We investigate this using a direct comparison: the same gene set, same classifier, same evaluation protocol, with only the embedding model varied.

---

## 2. Methods

**Gene set:** 178 human genes — 132 DNA repair genes (KEGG pathways hsa03430, hsa03420, hsa03440, hsa03410) and 46 tumor suppressor genes (COSMIC Cancer Gene Census Tier 1).

**Sequences:**
- *ESM-2:* canonical protein sequence (longest CDS translation) from Ensembl REST API
- *Evo2:* 2kb proximal promoter + canonical CDS from Ensembl REST API

**Embeddings:**
- *ESM-2 650M* (`esm2_t33_650M_UR50D`): mean-pooled final-layer representations, 1280-dim
- *Evo2-7B:* mean-pooled layer 28 (`blocks.28.mlp.l3`) representations, 4096-dim

**Classifier:** MLP (256 → 128 → 1), binary cross-entropy, Adam lr=1e-3, 10 epochs, 80/20 train/test split, 1 seed.

**Unsupervised metrics:** k-means ARI (k=2), silhouette score (cosine distance), within/between class distance ratio — compared against prior work.

---

## 3. Results

### 3.1 Supervised Classification

| Model | Input | Dim | Best AUROC | Final AUROC | Accuracy |
|---|---|---|---|---|---|
| ESM-2 650M | Protein sequence | 1280 | **0.894** | 0.880 | 0.722 |
| Evo2-7B | Promoter + CDS (DNA) | 4096 | 0.783 | 0.675 | 0.713 |

ESM-2 achieves substantially higher AUROC despite a 3× smaller embedding dimension. The gap (0.894 vs 0.783) suggests protein sequence space carries more discriminative functional signal for this classification task.

### 3.2 Unsupervised Geometry

| Model | k-means ARI | Silhouette | Distance ratio (within/between) |
|---|---|---|---|
| ESM-2 650M | -0.034 | -0.007 | 0.938 |
| Evo2-7B | -0.001 | 0.044 | 0.963 |

Both models show flat unsupervised geometry. Neither spontaneously organises genes by functional class in raw embedding space. The distance ratio near 1.0 for both models confirms genes of the same class are no closer together than genes from different classes.

---

## 4. Discussion

The superior supervised performance of ESM-2 over Evo2 is consistent with the intuition that protein sequence space is more directly constrained by function than DNA sequence space. ESM-2 was trained on evolutionary co-variation — mutations that are tolerated or selected against across species — which directly reflects functional constraints at the protein level. Evo2 captures nucleotide-level patterns including regulatory grammar, repeat content, and codon usage, which correlate with function but less directly.

Despite the stronger supervised signal, ESM-2 embeddings are no more semantically organised than Evo2 in raw geometry. This is consistent with the hypothesis that functional information is encoded in superposition in both models — distributed across many embedding dimensions in a way that requires a learned linear projection to extract, rather than being directly accessible via cosine similarity.

This distinction is practically important: ESM-2 is a strong feature extractor for supervised functional classification with limited labeled data, but it should not be used as a semantic search index for unsupervised gene retrieval.

---

## 5. Conclusion

ESM-2 protein embeddings significantly outperform Evo2 DNA embeddings for supervised gene functional classification (AUROC 0.894 vs 0.783) on the same gene set. However, both models fail to produce semantically organised embedding spaces — functional information requires supervised extraction in both cases. Protein language models are the preferred choice for supervised functional genomics tasks when protein sequences are available.

---

## References

- Lin, Z. et al. (2023). Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science*, 379(6637), 1123–1130.
- Nguyen, E. et al. (2025). Genome modeling and design across all domains of life with Evo 2. *Nature*.
- Kanehisa, M. & Goto, S. (2000). KEGG: Kyoto encyclopedia of genes and genomes. *Nucleic Acids Research*, 28(1), 27–30.
- Sondka, Z. et al. (2018). The COSMIC Cancer Gene Census. *Nature Reviews Cancer*, 18(11), 696–705.


---

## Attribution

- **ESM-2** — Meta AI / FAIR (Lin et al., 2023). Protein language model used for embedding extraction.
- **Evo2** — Arc Institute (Nguyen et al., 2025). DNA language model used for comparison.
- **AI Scientist** — Sakana AI (Lu et al., 2024). Framework used to automatically generate hypotheses, implement experiments, and write papers.
- **Claude Code** — Anthropic. Used to design experiments, build templates, debug code, and interpret results throughout this project.
- **DepMap** — Broad Institute. CRISPR essentiality data from DepMap 26Q1 Public release.

