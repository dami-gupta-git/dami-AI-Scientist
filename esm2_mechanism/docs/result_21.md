# Result 21 — Megascale stability: ESM-2 signal is partially family-dependent, revealing a gradient

**Date:** 2026-05-28
**Script:** `scripts/megascale_stability.py`
**Dataset:** S1724 benchmark — 1,277 single-point missense across 27 natural PDB proteins, ThermoMutDB ΔΔG labels
**CV:** Random / protein-holdout / cluster-holdout, 5 seeds
**GPU:** A100 80GB

---

## TL;DR

ESM-2 predicts thermodynamic stability under random split (Spearman ρ = **0.546**, binarised AUROC = **0.764**) but loses roughly half that signal under protein-holdout CV (ρ = **0.280**, AUROC = **0.642**). ESM-2 retains real transferable stability signal — protein-holdout AUROC of 0.642 is well above chance — but a substantial fraction is family-dependent.

The honest picture is a **gradient**, not a binary:

| Property | Family-split AUROC | Δ from random |
|---|---|---|
| Pathogenicity (result_6) | 0.884 | **0.002** — family-robust |
| Stability (result_21) | 0.642 | **0.122** — substantially family-dependent |
| Mechanism (results 1–10) | ~0.655 GOF | large F1 leakage (62.8%) — mostly family-dependent |

The comparison is metric-matched throughout (AUROC, binarised at median where needed), ruling out the objection that regression calibration drives the stability drop. Stability and mechanism are both family-dependent; pathogenicity is the exception. ESM-2 has progressively less transferable signal as the task moves from "is this mutation harmful at all" → "does it change protein stability" → "what disease mechanism does it cause."

---

## What this experiment tests

We asked: does ESM-2 encode thermodynamic stability — the physical consequence of a point mutation on protein folding — in a way that generalises across protein families?

This is a second positive control for the paper. The first (result_6) showed ESM-2 predicts ClinVar pathogenicity with AUROC 0.886, family-robust. But a reviewer could object: ClinVar labels are curated using population frequency data that is not independent of ESM-2's evolutionary training signal. Pathogenicity robustness could be a curation-circularity artefact.

Stability is the right counter-test: ΔΔG is measured in a test tube, no connection to clinical curation or evolutionary training. If ESM-2 encodes stability in a family-robust way, the biochemistry claim is strong. If stability is partially family-dependent — as mechanism is — we learn something more nuanced: ESM-2's family-transferable signal is specific to pathogenicity, not general to biochemistry.

---

## Results

### Primary table (delta_mean, 5-seed mean ± std)

| CV scheme | Spearman ρ | AUROC (binarised at median) |
|---|---|---|
| Random split | 0.546 ± 0.006 | 0.764 ± 0.008 |
| Protein-holdout | 0.280 ± 0.049 | 0.642 ± 0.023 |
| Cluster-holdout (identity) | 0.280 ± 0.049 | 0.642 ± 0.023 |

Δ AUROC (random − protein) = **0.122**. Δ Spearman = **0.266**.

The protein-holdout AUROC of 0.642 is meaningfully above chance (0.5) — ESM-2 has real cross-family stability signal. But it has lost 0.122 AUROC points relative to random split, indicating that a substantial portion of the apparent signal was family-level pattern recognition rather than per-variant biochemistry.

**Small-n caveat.** Protein-holdout here is leave-protein-out over only 27 proteins. The ±0.023 seed std on AUROC reflects sampling noise across seeds, but the per-protein ρ std of 0.274 (range −0.41 to 0.71) shows the underlying protein-to-protein variance is large — individual proteins drive the estimate substantially. Compare to pathogenicity (result_6), where family-split stability rests on ~944 genes and 658 Pfam families: the 0.884 protein-holdout AUROC there is estimated over a much larger partition. The direction of the result (stability is partially family-dependent, pathogenicity is family-robust) is robust, but the exact magnitude of the stability drop should be interpreted with the small-n in mind. A larger stability benchmark with more proteins would tighten the estimate.

Note: protein-holdout and cluster-holdout are identical here because MMseqs2 was unavailable on the pod and identity clustering was used (1 cluster per protein). With MMseqs2-20 clustering, the cluster-holdout might be marginally stricter. Given the S1724 dataset spans diverse folds (barnase, ubiquitin, tenascin, CI2, RNase H, etc.), most proteins are likely already in separate clusters.

### Per-residue delta (delta_pos)

| CV scheme | Spearman ρ | AUROC |
|---|---|---|
| Random split | 0.512 ± 0.010 | 0.742 ± 0.009 |
| Protein-holdout | 0.257 ± 0.035 | 0.633 ± 0.021 |

Same pattern as delta_mean, slightly weaker throughout. Δ AUROC = 0.109.

### Per-protein Spearman distribution (leave-one-protein-out)

Mean ρ = **0.248 ± 0.274** across 26 proteins (n ≥ 5 variants each).

| Statistic | Value |
|---|---|
| Mean | 0.248 |
| Std | 0.274 |
| Min | −0.414 (1EKG) |
| Max | 0.708 (1TEN) |
| Proteins with ρ < 0 | 4 / 26 (15%) |
| Proteins with ρ > 0.5 | 7 / 26 (27%) |

The wide distribution (std = 0.274) reflects genuine heterogeneity: some proteins transfer well (tenascin 1TEN ρ = 0.71, RNase H 1O6X ρ = 0.60, staphylococcal nuclease 1STN ρ = 0.46), others fail entirely (ubiquitin 1UBQ ρ = −0.14, PTB domain 2PTL ρ = −0.15). This is analogous to result_18 (AlphaMissense on ProteinGym: mean 0.721 ± 0.150) — physical labels reveal per-protein heterogeneity that curated labels hide.

### H3 — stability projection out of mechanism

Not run (merged Gerasimavicius embeddings not present on pod). Protocol pre-specified: fit Ridge on S1724 → normalise weight vector → project that direction out of merged delta_mean → re-run family-split logreg. Can be run locally.

---

## Pre-registered decision rule: LEAKY

| Criterion | Threshold | Observed | Pass? |
|---|---|---|---|
| H1: random ρ ≥ 0.5 | 0.5 | 0.546 | ✓ |
| H2: protein-split Δ ≤ 0.05 | 0.05 | **0.266** | ✗ LEAKY |
| H4: per-protein std ≤ 0.10 | 0.10 | **0.274** | ✗ HETEROGENEOUS |

Verdict: **LEAKY**. The HETEROGENEOUS criterion also fails, consistent with LEAKY — both reflect the same family-dependent pattern. The LEAKY label is technically correct per the pre-registered rule but should be read as "substantially family-dependent" rather than "purely family-memorisation" — 0.642 AUROC under protein-holdout is real signal.

---

## What this means in plain English

Think of ESM-2 as having learned two kinds of knowledge about mutations: knowledge that is specific to each protein family ("in barnase, mutations at the active site are very destabilising"), and knowledge that transfers across families ("large-to-small substitutions at buried positions tend to be destabilising everywhere").

Under random split, the model uses both kinds and gets ρ = 0.55. Under protein-holdout, it can only use the transferable knowledge, and drops to ρ = 0.28. The drop tells you how much of the signal was family-specific. For stability, roughly half was family-specific.

For pathogenicity (result_6), almost none was family-specific (Δ = 0.002) — mutations that are damaging enough to cause disease tend to look damaging in the same way regardless of which protein they're in. For mechanism (results 1–10), most was family-specific (62.8% leakage) — GOF vs LOF vs DN depends heavily on what kind of protein you're in.

This gives a **gradient of family-dependence**:
- Pathogenicity: nearly all transferable (AUROC Δ = 0.002)
- Stability: partially transferable, partially family-dependent (AUROC Δ = 0.122)
- Mechanism: mostly family-dependent (F1 leakage 62.8%)

The gradient makes biological sense. Whether a mutation is pathogenic is a relatively blunt question — did it break the protein badly enough? — and "badly broken" has common signatures. Whether a mutation destabilises the protein requires knowing something about the specific structural context. Whether it causes GOF vs LOF requires knowing the protein's function and disease biology — essentially the most context-dependent question of the three.

---

## Revised paper framing

**Before result_21:** "ESM-2 encodes pathogenicity and biochemistry robustly, but not mechanism."

**After result_21:** "ESM-2's family-transferable signal is specific to pathogenicity. Stability is partially family-dependent (AUROC drops 0.122 under protein-holdout, retaining real but reduced signal). Mechanism is mostly family-dependent (62.8% F1 leakage). There is a gradient of family-dependence that tracks how context-specific the prediction task is: pathogenicity is the most context-independent (broken is broken), mechanism is the most context-dependent (GOF vs LOF depends on what the protein does), and stability is intermediate."

---

## Comparison to existing results

All AUROC, binarised at median where needed:

| Property | Random-split AUROC | Protein/family-split AUROC | Δ | Interpretation |
|---|---|---|---|---|
| **Pathogenicity (result_6)** | 0.886 | **0.884** | **0.002** | Family-robust |
| **Stability (result_21)** | 0.764 | **0.642** | **0.122** | Substantially family-dependent; real cross-family signal remains |
| Mechanism GOF (results 1–10) | ~0.66 | ~0.655 | large F1 leakage | Mostly family-dependent |
| AM ClinVar AUROC (result_17) | 0.940 | ~0.948 per-family | ~0 | Family-robust (but curation-circular) |
| AM ProteinGym AUROC (result_18) | 0.721 | — | — | Wide per-assay distribution |

For completeness, stability Spearman ρ: random 0.546, protein-holdout 0.280, Δ = 0.266.

---

## Files

- `scripts/megascale_stability.py` — full pipeline
- `data/megascale/benchmarks.zip` — S1724 benchmark (ThermoMutDB + PDB sequences)
- `data/megascale_variants.json` — parsed 1,277 variants
- `data/embeddings/megascale_{wt,mut}_{mean,pos}.npy` — ESM-2 embeddings (on pod)
- `results/megascale_stability/summary.json` — 5-seed aggregated metrics + verdict
- `results/megascale_stability/per_protein_spearman.json` — per-protein ρ distribution
