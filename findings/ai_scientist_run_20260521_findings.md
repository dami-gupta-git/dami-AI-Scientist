# AI Scientist Run — Findings (May 21-22, 2026)

## Overview

We ran the AI Scientist framework on a new experiment template (`evo2_function`) designed to test whether frozen Evo2-7B embeddings encode functional gene class information. The pipeline automatically generated hypotheses, implemented them, ran experiments, and attempted to write papers. This documents the raw experimental findings.

**Dataset:** 178 human genes — 132 DNA repair (KEGG: MMR, NER, HR, BER pathways) and 46 tumor suppressors (COSMIC Cancer Gene Census). For each gene: 2kb proximal promoter + canonical CDS fetched from Ensembl, passed through frozen Evo2-7B, mean-pooled embeddings extracted from layer `blocks.28.mlp.l3` (4096 dimensions).

---

## Baseline Result

**Experiment:** Frozen Evo2-7B embeddings → MLP classifier (256 → 128 → 1) trained to distinguish DNA repair vs tumor suppressor genes. 80/20 train/test split, 3 random seeds.

| Metric | Mean (3 seeds) | Best seed |
|---|---|---|
| AUROC | 0.78 | 0.82 |
| Accuracy | 0.77 | 0.83 |
| F1 | 0.21 | 0.42 |

**Interpretation:** A frozen model with no functional supervision produces embeddings that a simple MLP can use to classify gene function with 0.78 AUROC. The low F1 reflects class imbalance (132 vs 46 genes) — AUROC is the reliable metric here.

---

## Experiment 1: Unsupervised Functional Clustering

**Hypothesis:** Do Evo2 embeddings spontaneously cluster by functional class without supervision?

**Method:** K-means (k=2) on embeddings, Adjusted Rand Index (ARI) and Normalized Mutual Information (NMI) against true labels, permutation null (1000 shuffles), within/between class cosine distances, silhouette score.

| Metric | Value | p-value | Interpretation |
|---|---|---|---|
| K-means ARI | ~0.00 | 0.516 | Not significant |
| NMI | ~0.00 | — | Negligible |
| Within-class distance | 0.059 | — | — |
| Between-class distance | 0.061 | — | — |
| Within/between ratio | 0.96 | — | Near 1.0 = no separation |
| Silhouette score | 0.044 | 0.132 | Not significant |

**Conclusion: Raw Evo2 embeddings do NOT cluster by functional class.** The geometry is essentially random with respect to function. This directly replicates the EvoMine finding that raw cosine similarity is not organized by biological function.

**Why this matters:** The MLP achieves 0.78 AUROC on the same embeddings. The functional signal exists but is in superposition — it requires a learned projection to extract, not raw distance.

---

## Experiment 2: Hierarchical Sub-pathway Structure

**Hypothesis:** Do DNA repair sub-pathways (MMR, NER, HR, BER) form coherent sub-clusters within the embedding space? Do "bridge genes" (genes annotated to both DNA repair and tumor suppression) sit geometrically between the two class centroids?

**Method:** Ward linkage hierarchical clustering, cophenetic correlation, per-sub-pathway mean within vs between cosine distances with permutation test (1000 iterations), bridge gene centroid distance analysis.

### Hierarchical Structure
- **Cophenetic correlation: 0.45** — moderate hierarchy exists but not strong (1.0 would be perfect)

### Sub-pathway Coherence

| Sub-pathway | N genes | Within dist | Between dist | p-value | Significant? |
|---|---|---|---|---|---|
| MMR | 23 | 0.064 | 0.062 | 0.625 | No |
| NER | 48 | 0.063 | 0.060 | 0.746 | No |
| HR | 32 | 0.063 | 0.060 | 0.635 | No |
| **BER** | **29** | **0.041** | **0.054** | **0.015** | **Yes ✅** |

**BER (base excision repair) is the only sub-pathway that forms a significantly coherent sub-cluster.** BER genes are notably more similar to each other (within distance 0.041) than to other repair genes (between distance 0.054). This is biologically interesting — BER is mechanistically distinct from the other repair pathways, operating on small chemical lesions rather than helix-distorting damage.

### Bridge Gene Analysis

Bridge genes identified (annotated to both DNA repair and tumor suppression):
**BRCA1, BRCA2, MLH1, MSH2, MSH6, PMS2, PALB2, ATM** (8 genes)

| Comparison | Mean cosine distance |
|---|---|
| Bridge genes → DNA repair centroid | 0.050 |
| Bridge genes → TSG centroid | 0.045 |
| Other TSGs → DNA repair centroid | 0.032 |
| Other TSGs → TSG centroid | 0.027 |

**t-statistic: 1.77, p=0.084** — trending but not significant (n=8 bridge genes is small).

Bridge genes are closer to the TSG centroid than to DNA repair centroid, but notably farther from both centroids than other TSGs are from their own centroid. This suggests bridge genes occupy an intermediate embedding position, consistent with their dual functional annotation — but the result needs more genes to reach significance.

---

## Summary of Key Findings

1. **Supervised signal is real:** MLP probe achieves 0.78 AUROC, confirming Evo2 encodes functional information extractable by learned projection.

2. **Unsupervised signal is absent:** Raw embedding geometry does not cluster by functional class (ARI ≈ 0, p=0.516). Confirms EvoMine finding — raw cosine similarity does not track function. The AI Scientist-generated paper frames this clearly: *"Evo2 embeddings contain features useful for supervised functional prediction, but do not encode functional semantics geometrically."*

3. **BER is anomalous:** Base excision repair genes (n=29) form a significantly coherent sub-cluster (p=0.015) — the only DNA repair sub-pathway that does. Within-class distance 0.041 vs between-class distance 0.054. MMR (p=0.625), NER (p=0.746), HR (p=0.635) show no significant coherence. Mechanism unclear — may reflect BER's mechanistic distinctiveness (small non-helix-distorting lesions vs. bulky adducts/DSBs) mapping to distinctive sequence features.

4. **Bridge genes trend toward intermediate positions:** BRCA1, BRCA2, MLH1, MSH2, MSH6, PMS2, PALB2, ATM — 8 genes annotated to both DNA repair and tumor suppression — sit closer to the TSG centroid (0.045) than to DNA repair centroid (0.050), with p=0.084. Underpowered with n=8 but directionally consistent with dual functional identity.

5. **The paper produced is high quality:** The AI-generated paper correctly frames the supervised vs. unsupervised distinction, uses appropriate statistics, and draws honest conclusions. Notable omission: the BER sub-pathway finding (from the hierarchical experiment) was not included since it ran as a separate experiment.

---

## Validation Experiment: Glycolysis vs T-cell Signaling

**Purpose:** Test whether the 0.78 AUROC generalises to a completely different functional contrast with no overlap with DNA repair or cancer biology.

**Dataset:** 189 human genes — 67 glycolysis (KEGG hsa00010) and 122 T-cell receptor signaling (KEGG hsa04660). Same pipeline: Evo2-7B layer 28, mean-pooled embeddings.

| Metric | Glycolysis vs T-cell | DNA repair vs TSG |
|---|---|---|
| Best AUROC | 0.66 | 0.78 |
| Centroid cosine distance | 0.0014 | ~0.006 |
| Within-class cosine sim | ~0.955 | ~0.941 |

**The supervised signal is weaker and the unsupervised geometry is more compressed.** Glycolysis and T-cell signaling genes are nearly indistinguishable in raw Evo2 embedding space (centroid distance 0.0014 vs 0.006 for DNA repair/TSG).

**Interpretation:** The 0.78 AUROC for DNA repair vs TSG is **not fully general** — Evo2 encodes some functional contrasts better than others. Genome integrity pathways (DNA repair, tumor suppression) appear to be specifically well-represented, likely because mutations in these genes have been under strong purifying selection across mammals, leaving clear evolutionary sequence signatures. Glycolysis genes are more conserved at the protein level but less distinctive at the regulatory/promoter level, making them harder to separate in embedding space.

---

## Overall Story

1. Evo2 embeddings encode functional information — but only extractable with supervision (MLP probe), not raw geometry
2. The signal strength is contrast-dependent: genome integrity (0.78 AUROC) >> metabolic/immune (0.66 AUROC)
3. This suggests Evo2 captures evolutionary constraint signals rather than general functional semantics
4. Raw cosine similarity is useless for functional retrieval in both contrasts — confirms EvoMine

---

## What This Means for EvoMine

The core EvoMine hypothesis — that functional information is in superposition in Evo2 embeddings and requires extraction — is confirmed by both contrasts. The additional finding is that the superposition is uneven: genome integrity genes carry stronger functional signal, possibly because their sequence signatures are more evolutionarily distinctive.

The BER finding is the most unexpected result and worth following up: why does BER cluster when MMR/NER/HR don't? Possible explanations:
- BER genes have distinctive promoter features (e.g., higher CpG density, specific TF binding patterns)
- BER operates on small non-helix-distorting lesions — mechanistic distinctiveness may map to sequence distinctiveness
- BER genes may share codon usage or expression patterns that Evo2 encodes

---

## Caveats

- Small test sets (36 genes for DNA repair/TSG, 38 for glycolysis/T-cell) — high variance
- Single Evo2 layer probed (blocks.28.mlp.l3) — layer ablation not completed
- Paper writeup phase did not complete for most experiments due to pipeline bugs
- The `run_0` results for several experiments contain baseline numbers due to Aider cache reuse — only `run_1`+ results are from the actual modified experiments
