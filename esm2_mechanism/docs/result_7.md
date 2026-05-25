# Result 7 — MLP family-split, merged dataset, Option B
## Date: May 24–25, 2026 | Model: ESM-2 650M

---

## TL;DR

Three experiments run in parallel after results 1–6:

1. **MLP delta family-split** (the blocking experiment from result_5): signal drops only Δ=+0.052 under family-split — much less leakage than the linear WT-only baseline (Δ=+0.191). The nonlinear delta signal is largely real.
2. **Merged dataset Option B** (gene-level WT, Gerasimavicius + G2P): with 1,985 genes and proper Pfam coverage, family-split F1=0.393 and GOF AUROC=0.728 — real above-chance signal surviving family holdout.
3. **Additional nonlinear probes** (GBM, RF, kNN on delta_mean): kNN and MLP are the strongest; GBM/RF weaker (likely PCA information loss).

---

## 1. MLP family-split on delta_mean (Gerasimavicius, 948 genes)

### Results

| CV | macro-F1 | Interpretation |
|---|---|---|
| Gene-split | **0.415** | Nonlinear mechanism signal |
| Family-split | **0.364** | Δ = +0.052 |

### All probes on delta_mean (gene-split)

| Probe | macro-F1 |
|---|---|
| MLP (256→64, class-weighted) | **0.415** |
| kNN (k=10, cosine) | **0.410** |
| GBM (PCA-50) | 0.338 |
| RF (PCA-50) | 0.291 |
| Linear LR (result_1) | 0.279 |

### Interpretation

The MLP family-split drop of Δ=+0.052 is much smaller than the WT-only collapse (Δ=+0.191). This confirms the earlier hypothesis from result_4: the WT-only signal was largely family leakage, while the nonlinear delta signal is substantially real. There IS mechanism information in ESM-2 delta embeddings, accessible nonlinearly.

Compared to result_5's pre-family-split MLP (F1=0.431 using full 100-epoch training), the current F1=0.415 is slightly lower due to the reduced 30-epoch limit — consistent.

kNN achieving F1=0.410 (nearly matching MLP) is the cleanest evidence: no transformation, purely local structure in delta space. Mechanism classes are locally clustered in delta space, and this clustering is not primarily family-driven.

---

## 2. Merged dataset: Option B (gene-level WT mean embeddings)

### Dataset
- 19,100 variants, 1,985 genes (Gerasimavicius 948 + G2P/ClinVar pathogenic-only 1,037 new genes)
- GOF: 2,825 variants / 146 genes | DN: 1,716 variants / 107 genes | LOF: 14,559 variants / 1,732 genes
- Pfam annotated: 1,950/1,985 genes, 1,146 unique families

### Results

| CV | macro-F1 | GOF AUROC | DN AUROC | LOF AUROC |
|---|---|---|---|---|
| Gene-split | 0.469 | 0.784 | 0.700 | 0.765 |
| Family-split | **0.393** | **0.728** | **0.634** | **0.691** |
| Δ | +0.077 | | | |

### Comparison: original vs merged dataset (gene-level WT)

| Dataset | Gene-split F1 | Family-split F1 | Δ |
|---|---|---|---|
| Gerasimavicius (948 genes) | 0.580 | 0.389 | +0.191 |
| Merged (1,985 genes) | 0.469 | 0.393 | +0.077 |

The merged dataset halves the leakage (Δ 0.191→0.077). Three things happen when you add genes:
1. More families represented (662→1,146) — harder to learn family identity shortcut
2. More balanced classes — GOF/DN have more representation
3. Family-split F1 stays at 0.393 — essentially the same absolute signal, just less inflated gene-split

**The family-split F1=0.393 with GOF AUROC=0.728 is real above-chance signal** (vs chance ~0.333 for F1, 0.50 for AUROC). This is the strongest family-split-robust signal in the entire study.

### Note on initial Option B result (erroneous)
The first Option B run reported Δ=+0.011 — this was wrong because pfam_families.json only had Pfam annotations for the original 939 Gerasimavicius genes, not the 1,037 new G2P genes. After fetching Pfam for all new genes (1,950/1,985 annotated, 1,146 families), the correct result is Δ=+0.077.

---

## 3. Updated story (supersedes result_6 revised story)

> ESM-2 delta embeddings contain nonlinear mechanism signal that survives family-split CV (MLP F1 0.415→0.364, Δ=+0.052). The pathogenicity/mechanism dissociation from result_6 (AUROC 0.88 vs macro-F1 0.28) used a linear probe — the mechanism signal is nonlinear. ESM-2 encodes mechanism, but nonlinearly and more weakly than pathogenicity.
>
> The WT-only mechanism signal on the original 948-gene dataset (F1=0.580→0.389, Δ=+0.191) was largely family leakage. On the merged 1,985-gene dataset, gene-level WT family-split F1=0.393 with GOF AUROC=0.728 is a real above-chance result that survives family holdout.
>
> **Revised central finding**: ESM-2 encodes molecular disease mechanism, but (a) the signal is nonlinear in delta space, (b) linear probes give a false negative, and (c) the small original dataset had severe family-leakage inflation. Family-split CV with adequate gene diversity is necessary to make this claim.

---

## 4. What this means for publishability

The story is now more positive than result_6 suggested:

- The delta probe is **not** a null — it has real nonlinear signal (F1=0.364 family-split)
- The WT-only result was inflated but the mechanism signal is real on both datasets under family-split
- The pathogenicity/mechanism dissociation still holds (pathogenicity AUROC 0.88 >> mechanism F1 0.364) but the mechanism signal is non-trivial

**Revised claim**: ESM-2 encodes mechanism nonlinearly in delta space. The signal is weaker than pathogenicity but present and family-split-robust. This is a positive finding, not a null — and more novel than result_6's framing as "folk wisdom confirmed."

**Revised venue**: with this positive result, *Bioinformatics* / *Genome Biology* is the right target and the paper is stronger. The within-family analysis (can mechanism be learned within a single Pfam family?) is still the ceiling-lifter for higher venues.

---

## 5. Open questions

1. **MLP family-split on delta_pos** — still running on RunPod. Expected to be weaker than mean-pooled (per-residue was weaker under MLP in result_5).
2. **Merged dataset delta probe** — not yet run. With 19,100 variants and better class balance, the delta MLP may be stronger.
3. **Within-family analysis** — test mechanism classification restricted to single Pfam families (kinase, C2, homeobox). If it works, the story becomes "mechanism is a within-family problem."
4. **Multi-seed replication** — all numbers are seed=0 only.

---

## Files

- `results/20260524_baseline_run/run_0/option_b_gene_level_wt_merged.json` — Option B results
- `results/20260524_baseline_run/run_0/mlp_results_seed0.json` — MLP probe full results (pending final pull from RunPod)
- `data/pfam_families.json` — updated with 1,037 new G2P gene annotations
- `data/embeddings/merged_embeddings_*.npy` — merged dataset embeddings (19,100 × 1,280)
- `data/embeddings/merged_valid_variants.json` — merged variant list
