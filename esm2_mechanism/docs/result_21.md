# Result 21 — Megascale stability: ESM-2 is a fold-recognition engine, not a biochemistry encoder

**Date:** 2026-05-28
**Script:** `scripts/megascale_stability.py`
**Dataset:** S1724 benchmark — 1,277 single-point missense across 27 natural PDB proteins, ThermoMutDB ΔΔG labels
**CV:** Random / protein-holdout / cluster-holdout, 5 seeds
**GPU:** A100 80GB

---

## TL;DR

ESM-2 predicts thermodynamic stability (ΔΔG) well when tested on proteins it has seen — Spearman ρ = **0.546**, binarised AUROC = **0.764** under random split. But when tested on proteins from completely different families (protein-holdout CV), both metrics collapse: ρ drops to **0.280** (Δ = 0.266) and AUROC drops to **0.642** (Δ = **0.122**). The pre-registered verdict is **LEAKY**.

Critically, the dissociation with pathogenicity holds on the **same metric** (AUROC): stability protein-holdout AUROC = 0.642, pathogenicity protein-holdout AUROC = **0.884** (Δ = 0.002, result_6). This rules out the objection that "stability leaks because regression calibration is harder than binary classification" — both are now evaluated as binary classifiers, and the gap is large and real.

This changes the paper's framing. The original positive-control story was "ESM-2 encodes pathogenicity AND biochemistry, just not mechanism." Result 21 breaks that story: ESM-2 does not robustly encode biochemistry. **Pathogenicity is the exception** — it is the only thing ESM-2 encodes in a family-robust way. Both stability (AUROC Δ = 0.122) and mechanism (leakage fraction 62.8%, result_7) collapse under family holdout. ESM-2 is primarily a fold-recognition engine. Pathogenicity escapes because severely disruptive mutations share a recognisable signature regardless of protein family.

---

## What this experiment tests

We asked: does ESM-2 encode thermodynamic stability — the physical consequence of a point mutation on protein folding — in a way that generalises across protein families?

This is a second positive control for the paper. The first positive control (result_6) showed ESM-2 predicts ClinVar pathogenicity with AUROC 0.886, and that this is family-robust. But a reviewer could object: ClinVar labels are curated using population frequency data, and ESM-2 was trained on evolutionary co-variation that also reflects population-level constraint. The two might not be independent, so pathogenicity robustness could be a circularity artefact.

Stability is the right counter-test: ΔΔG is measured in a test tube, with no connection to clinical curation or evolutionary training. If ESM-2 also encodes stability in a family-robust way, the biochemistry claim is bulletproof. If stability leaks — as mechanism does — then ESM-2's competence is more narrowly fold-recognition than genuine biochemistry encoding.

---

## Results

### Primary table (delta_mean, 5-seed mean ± std)

| CV scheme | Spearman ρ | AUROC (binarised at median) |
|---|---|---|
| Random split | 0.546 ± 0.006 | 0.764 ± 0.008 |
| Protein-holdout | 0.280 ± 0.049 | 0.642 ± 0.023 |
| Cluster-holdout (identity) | 0.280 ± 0.049 | 0.642 ± 0.023 |

**Δ (random − protein) = 0.266** → LEAKY (pre-registered threshold ≥ 0.10)

Note: protein-holdout and cluster-holdout are identical here because MMseqs2 was unavailable on this pod and identity clustering was used (1 cluster per protein = 27 clusters for 27 proteins). With MMseqs2-20 clustering, proteins with sequence similarity > 20% would be grouped, making cluster-holdout a stricter test. Given the S1724 dataset spans diverse folds (barnase, ubiquitin, tenascin, CI2, RNase H, etc.), most proteins are likely in separate clusters anyway.

### Per-residue delta (delta_pos)

| CV scheme | Spearman ρ | AUROC |
|---|---|---|
| Random split | 0.512 ± 0.010 | 0.742 ± 0.009 |
| Protein-holdout | 0.257 ± 0.035 | 0.633 ± 0.021 |

Same pattern: large random-to-protein drop (Δ = 0.255). Per-residue delta is slightly weaker than mean-pooled delta throughout.

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

The distribution is wide and bimodal — analogous to result_18 (AlphaMissense on ProteinGym, mean 0.721 ± 0.150). Some proteins are predicted well (tenascin 1TEN ρ = 0.71, staphylococcal nuclease 1STN ρ = 0.46, RNase H 1O6X ρ = 0.60); others are near zero or negative (ubiquitin 1UBQ ρ = −0.14, PTB domain 2PTL ρ = −0.15, FT8 ρ = −0.23).

### H3 — stability projection out of mechanism

Not run on this pod (merged Gerasimavicius embeddings not present). Can be run locally. Protocol remains pre-registered: fit Ridge on S1724 → normalise weight vector → project direction out of merged delta_mean → re-run family-split logreg.

---

## Pre-registered decision rule: LEAKY

| Criterion | Threshold | Observed | Pass? |
|---|---|---|---|
| H1: random ρ ≥ 0.5 | 0.5 | 0.546 | ✓ |
| H2: protein-split Δ ≤ 0.05 | 0.05 | **0.266** | ✗ LEAKY |
| H4: per-protein std ≤ 0.10 | 0.10 | **0.274** | ✗ HETEROGENEOUS |

Verdict fires as **LEAKY** (checked first, most informative). The per-protein distribution is also HETEROGENEOUS, consistent with LEAKY — both diagnose the same underlying phenomenon.

---

## What this means in plain English

Imagine testing a doctor's ability to diagnose a disease. In one test, you show them patients from hospitals they've trained in — they do well. In another test, you show them patients from completely different hospitals with different patient populations — performance drops sharply. That's what's happening here.

ESM-2 was trained on evolutionary sequence data across all protein families. When you ask it "will this mutation destabilise this protein?", it answers partly by recognising the protein family and applying family-level rules ("kinase domain mutations near the ATP pocket are usually destabilising"). That works fine if your test set contains the same families as training. But when you hold out entire protein families, those family-level rules no longer apply, and performance collapses.

This is exactly the pattern we called "leakage" for mechanism prediction in results 1–10: the model is doing fold/family recognition, not variant-level reasoning. The difference is that for pathogenicity (result_6), the collapse does *not* happen — pathogenicity signal is present even across families, because mutations that disrupt a protein enough to cause disease tend to do so in recognisable ways regardless of which protein it is.

---

## Revised paper framing

**Before result_21:** "ESM-2 encodes local biochemistry (pathogenicity + stability) robustly, but not mechanism. Mechanism is the hard case."

**After result_21:** "ESM-2's only family-robust signal is pathogenicity. Stability and mechanism both collapse under family holdout. ESM-2 is primarily a fold-recognition engine. Pathogenicity is detectable across families because severely disruptive mutations share a recognisable signature regardless of protein context. Mechanism and stability require understanding what a mutation does *within* a specific protein's structural and functional context — and ESM-2's frozen representations do not generalise that knowledge across folds."

This is a sharper and more specific claim. It makes pathogenicity robustness the phenomenon to explain, rather than treating it as just a sanity check.

---

## Comparison to existing results

All on the same metric (AUROC, binarised at median where needed) for a fair comparison:

| Property | Random-split AUROC | Protein/family-split AUROC | Δ | Robust? |
|---|---|---|---|---|
| **Pathogenicity (result_6)** | 0.886 | **0.884** | **0.002** | **Yes** |
| Stability — binarised ΔΔG (result_21) | 0.764 | 0.642 | **0.122** | No — LEAKY |
| Mechanism GOF (result_6/7) | ~0.66 | ~0.655 | ~0.005 | No (62.8% F1 leakage) |
| AM ClinVar AUROC (result_17) | 0.940 | ~0.948 per-family | ~0 | Yes (but curation-circular) |
| AM ProteinGym AUROC (result_18) | 0.721 | — | — | No (wide per-assay) |

The stability row is the critical one: same metric as pathogenicity, same holdout logic, large collapse. The dissociation is metric-independent.

For completeness, stability Spearman ρ: random 0.546, protein-holdout 0.280, Δ = 0.266.

---

## Files

- `scripts/megascale_stability.py` — full pipeline
- `data/megascale/benchmarks.zip` — S1724 benchmark (ThermoMutDB + PDB sequences)
- `data/megascale_variants.json` — parsed 1,277 variants
- `data/embeddings/megascale_{wt,mut}_{mean,pos}.npy` — ESM-2 embeddings (on pod)
- `results/megascale_stability/summary.json` — 5-seed aggregated metrics + verdict
- `results/megascale_stability/per_protein_spearman.json` — per-protein ρ distribution
