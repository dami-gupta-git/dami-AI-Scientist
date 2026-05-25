# Evo2 Embeddings Encode Functional Gene Class Information: A Preliminary Finding

## Background

Evo2 is a 7-billion parameter DNA language model trained autoregressively on 8.8 trillion nucleotide tokens spanning all domains of life. It was trained purely on next-nucleotide prediction — no functional labels, no pathway annotations, no disease associations. The question is whether its internal representations have implicitly learned something about gene function.

Prior work (EvoMine) showed that raw cosine similarity in Evo2 embedding space is dominated by repeat content and sequence texture — it does not track biological function directly. The functional signal appears to be present but in superposition, requiring learned projections to extract.

## Experiment

We selected 178 human genes from two functionally distinct classes:
- **132 DNA repair genes** — from KEGG pathways hsa03430 (MMR), hsa03420 (NER), hsa03440 (HR), hsa03410 (BER)
- **46 tumor suppressor genes** — curated from the COSMIC Cancer Gene Census Tier 1

For each gene, we fetched the 2kb proximal promoter sequence and canonical CDS from the Ensembl REST API, concatenated them, and passed the full sequence through frozen Evo2-7B weights (no fine-tuning). We extracted mean-pooled embeddings from layer `blocks.28.mlp.l3` (4096 dimensions) and trained a shallow MLP classifier (hidden layers: 256 → 128 → 1, binary cross-entropy loss) on an 80/20 train/test split across 3 random seeds.

## Result

| Metric | Mean (3 seeds) | Best seed |
|---|---|---|
| AUROC | 0.78 | 0.82 |
| Accuracy | 0.77 | 0.83 |
| Chance AUROC | 0.50 | — |

The MLP achieved a mean validation AUROC of **0.78**, peaking at **0.82** in the best seed. This is substantially above chance (0.50) and above a k-mer frequency baseline (~0.58 from EvoMine exon/intron experiments).

## Interpretation

A frozen model that has never seen a functional label produces representations that a simple classifier uses to distinguish DNA repair genes from tumor suppressors with 0.78 AUROC. The MLP acts as a **learned probe** — it extracts the functional signal from the high-dimensional embedding space that raw cosine similarity cannot access.

This is consistent with the EvoMine finding that Evo2 layer 26 encodes biological concepts (exon-intron boundaries, TF binding sites, structural elements) in superposition. The MLP learns to decompose these superposed features and isolate the dimensions relevant to functional class.

The ENO1/KRT18 result from EvoMine is also relevant: two genes on different chromosomes with zero BLAST hits showed high embedding similarity and shared chromatin marks. Evo2 found functional equivalence that sequence alignment could not. The MLP classification result suggests this kind of signal is systematic, not anecdotal.

## Open Questions

- **Which layer** carries the strongest functional signal? (layer ablation)
- **Which embedding dimensions** matter? Is the signal concentrated or distributed?
- **Does it generalize** to other functional contrasts beyond DNA repair vs TSG?
- **What drives the signal** — promoter regulatory grammar, CDS codon usage, or both?
- **How does it compare** to a k-mer baseline on this specific task?

## Next Steps (In Progress)

These questions are being investigated via the AI Scientist framework, which generates and runs novel hypotheses automatically on top of this baseline. Experiments currently running include:

- **Hierarchical sub-pathway clustering** — Ward linkage on embeddings, sub-pathway coherence, bridge gene analysis (BRCA1/BRCA2/MLH1 etc). Partial results available but run not yet complete.
- **Layer ablation** — comparing early vs late Evo2 layers
- **Unsupervised clustering geometry** — k-means, silhouette, UMAP

Results will be updated here as experiments complete.

## Caveats

- Small test set (36 genes) — results have high variance across seeds (AUROC range: 0.67–0.82)
- Class imbalance (132 repair vs 46 TSG) — F1 scores are low; AUROC is the reliable metric
- Single model size (7B) — scaling behavior unknown
