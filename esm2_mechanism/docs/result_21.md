# Result 21 — Megascale stability as a second ESM-2 positive control

**Date:** 2026-05-27 (pre-registration)
**Status:** Not yet run — waiting on result_20 (perturbation scan) before launching
**Script:** `scripts/megascale_stability.py`
**Dataset:** S1724 benchmark — 1,277 single-point missense across 27 natural PDB proteins, ThermoMutDB ΔΔG labels
**CV:** Random / protein-holdout / MMseqs2-cluster-holdout, 5 seeds
**Outputs:** `results/megascale_stability/{summary,per_protein_spearman,h3_stability_projection}.json`

---

## Motivation

Result_6 gives one positive control: pathogenicity (AUROC 0.886, family-split-stable). But ClinVar curation and population-frequency labels are not fully independent (result_17/18 caveat). A second positive control with physical-ground-truth labels and no curation circularity is needed to rule out "the pipeline works only on labels that leak from training data."

S1724 is the right dataset: 27 real PDB proteins, ThermoMutDB curated single-point ΔΔG from the original biophysics literature, no overlap with ESM-2 training labels.

**Sequencing note:** This experiment is held pending result_20 (perturbation scan, 567k probes). If result_20 passes G1 (F1 > 0.368), the paper acquires a positive mechanism story and result_21 becomes defensive reviewer-satisfaction. If result_20 fails G1, the mechanism null hardens and result_21 becomes load-bearing for the "pipeline is sound" argument. Decision to launch will be made after result_20 reads out.

---

## Pre-registered hypotheses

**H1 — Stability is encoded.**
ESM-2 delta_mean predicts S1724 ΔΔG at Spearman ρ ≥ 0.5 under random split.
Prior: very likely to pass. ThermoMPNN (supervised) gets ~0.7; ESM-2 zero-shot log-likelihood correlates with DMS fitness; the random-split number is essentially a sanity check.

**H2 — Stability is family-robust.**
ρ drops by ≤ 0.05 (absolute) under protein-holdout CV.
Prior: uncertain. This is the informative test. If ESM-2's stability signal is partly family-memorisation (same as mechanism leakage), the drop will be ≥ 0.10 (LEAKY verdict) — a finding that would substantially reshape the project's central claim.

**H3 — Stability direction does not rescue mechanism prediction.**
Protocol (fully pre-specified):
1. Fit Ridge on S1724 delta_mean → ΔΔG.
2. Extract stability projection vector v = normalised Ridge weight vector.
3. Compute residuals: delta_mean_residual = delta_mean − (delta_mean · v)v  (remove one direction).
4. Re-run family-split logistic regression on residuals, merged dataset, 5 seeds.
5. H3 passes if projected F1 ≤ baseline F1 + 0.01 — the stability direction carries no mechanism signal beyond the raw delta.

**H4 — Per-protein ρ distribution is tight (std ≤ 0.10).**
If std ≥ 0.15, verdict is HETEROGENEOUS — analogous to result_18 (AM on ProteinGym). This would mean ESM-2's stability signal is protein-specific, not uniformly competent, and the curation-vs-physical-label distinction matters for per-stratum variance.

---

## Decision table (ordered by informativeness)

| Verdict | Condition | Interpretation | Paper impact |
|---|---|---|---|
| **LEAKY** | ρ ≥ 0.5, protein-split Δ ≥ 0.10 | Stability signal partly family-memorisation — same shortcut as mechanism | Major: central claim must distinguish "ESM-2 learns fold identity" from "ESM-2 learns biochemistry" |
| **HETEROGENEOUS** | ρ ≥ 0.5, Δ ≤ 0.05, per-protein std ≥ 0.15 | Works on average, fails on some proteins — physical labels reveal per-protein variance that curation hides | Moderate: sharpens result_17/18 comparison; per-stratum evaluation is the right two-axis view |
| **ROBUST** | ρ ≥ 0.5, Δ ≤ 0.05, per-protein std ≤ 0.10 | Expected; second positive control holds | Strengthens "pipeline is sound": ESM-2 encodes biochemistry (pathogenicity AND stability), not mechanism |
| **WEAK** | ρ 0.3–0.5 | Partial stability signal | Moderately informative; report and stop |
| **NULL** | ρ < 0.3 | Very unexpected | Would undermine framing: ESM-2 competence may be curation-specific, not biochemical |

---

## Dataset note: S1724 vs main Megascale CSV

The plan originally referenced the full Megascale CSV (~2.3M rows, synthetic mini-proteins from Rocklin/Baker de novo design). S1724 is the curated benchmark from `benchmarks.zip` — 27 **natural PDB proteins** with literature ΔΔG values. This matters for the Pfam-coverage concern: all 27 S1724 proteins are real, annotated domains; there is no mini-protein / no-Pfam subset to worry about. The protein-holdout split is clean throughout.

If the result warrants it, a supplementary run on the full Megascale CSV (restricted to the dmsv7 Cho 2026 recalibrated subset, natural proteins only, MMseqs2-cluster-split) can be added as an appendix.

---

## Files (to be created on run)

- `results/megascale_stability/summary.json` — 5-seed aggregated metrics + verdict
- `results/megascale_stability/per_protein_spearman.json` — leave-one-protein-out ρ per protein
- `results/megascale_stability/h3_stability_projection.json` — H3 baseline vs projected F1
- `data/embeddings/megascale_{wt,mut}_{mean,pos}.npy` — cached ESM-2 embeddings
