# EvoMine — Experimental Findings

## What EvoMine Is

EvoMine treats Evo2's intermediate-layer embeddings as a **searchable retrieval index of functional biology**. Instead of using embeddings for classification (the standard approach), the hypothesis is that cosine similarity in the latent space reflects functional relationships — a "semantic BLAST." For variants of uncertain significance (VUS), this reframes interpretation: rather than scoring a variant in isolation, retrieve its functional neighbors and use their known properties as evidence.

The foundational question: is the functional information in Evo2's embedding space **geometrically accessible** via simple cosine distance, or is it in superposition and requires decomposition to extract?

---

## What Was Tested

### Setup
- 25 human genes, 5 functional categories
- 512bp sliding windows across 10kb gene regions
- Evo2-7B embeddings extracted at multiple intermediate layers
- Pairwise cosine similarity computed with and without repeat filtering
- Compared against BLAST sequence similarity as a baseline

### Key Finding 1: Raw Cosine Similarity Is Repeat-Dominated

Most strong embedding similarities were driven by common repetitive elements (Alu repeats). After filtering, the majority of high-similarity pairs between different genes showed **zero detectable sequence similarity** (BLAST hits = 0) — but this was not evidence of functional similarity, it was noise from shared repeat content.

**Conclusion: Raw cosine similarity does not reliably track biological function.**

### Key Finding 2: The Signal Exists But Is Rare

One clean, interpretable example emerged after repeat filtering:

- A promoter window in **VIM** (vimentin, chr10) showed cosine similarity of **0.948** with a promoter window in **DES** (desmin, chr2)
- Zero sequence homology (BLAST hits = 0)
- Both are active promoters in muscle and connective tissue cells
- Both share CTCF and POLR2A binding
- Both show high cross-species conservation (phastCons 0.84–0.89) — independently of Evo2
- VIM and DES are well-known co-expressed genes that function together

This shows Evo2 can surface functional equivalence that sequence alignment cannot see. But this signal is rare and only surfaces after heavy filtering.

### Key Finding 3: Learned Probes Outperform Raw Similarity

Evo2 layer 26 embeddings achieved **0.875 AUC** for exon/intron classification vs a **0.583** k-mer + GC baseline. Linear regression consistently outperformed other approaches in downstream tasks.

**Conclusion: The functional signal is present but in superposition. Learned projections (linear probes, MLP classifiers) extract it; raw cosine similarity cannot.**

This is consistent with Goodfire + Arc Institute's finding via SAEs that Evo2 layer 26 encodes biological concepts (exon-intron boundaries, TF binding sites, protein structural elements) in superposition.

---

## What Was Not Tested

- None of the 6 use cases (VUS resolution, cross-gene analogy, non-coding triage, drug repurposing, evolutionary outlier detection, tumor signature fingerprinting) were directly evaluated
- The scalar vs. vector variant effect distinction — existing tools collapse variant impact to a single log-likelihood delta; the vector approach preserving directional disruption information — remains untested
- The BRCA1 saturation mutagenesis retrieval experiment (FAISS index + k-NN pathogenicity prediction) was designed but not run
- The agentic vs. fixed pipeline comparison was designed but not run

---

## Open Questions

1. Does repeat-filtering + learned projection recover reliable functional retrieval, or does the signal remain too noisy?
2. Which Evo2 layer carries the most linearly extractable functional signal?
3. Does the vector (directional) variant effect contain information the scalar log-likelihood misses?
4. Can latent-space k-NN outperform sequence-space k-NN for VUS pathogenicity prediction?
5. Does N-masking repeats preserve the functional signal or destroy it?

---

## Related Work

The AI Scientist experiment (separate project) extends EvoMine's probe finding to gene-level functional classification:
- 178 human genes (132 DNA repair, 46 tumor suppressors)
- Frozen Evo2-7B layer 28 embeddings + MLP classifier
- Baseline AUROC: **0.78** (confirming learned probes extract functional signal)
- Further experiments in progress: layer ablation, sub-pathway clustering, bridge gene analysis

See [evo2_function_classification_finding.md](evo2_function_classification_finding.md) for details.
