# Publication plan — bioRxiv with versioned releases

The plan is to post a minimal v1 quickly (priority date + early feedback), then add controls and generalisation in v2 and v3 as experiments complete. Each version is a self-contained scientific claim; later versions strengthen rather than replace earlier ones.

---

## v1 — minimal frozen-probe GOF finding (target: ~1 week)

### Working title
*"Frozen ESM-2 embeddings retain a cross-family gain-of-function signal without fine-tuning"*

### Abstract (draft)
ESM-Effect / ESMGain (Tang et al. 2025) showed that fine-tuning ESM-2 enables prediction of gain-of-function effects in deep mutational scanning data. We show that this signal is already present in **frozen ESM-2 embeddings**, recoverable by a simple probe over mean-pooled per-gene representations at the clinical disease level. On 1,985 human disease genes drawn from Gerasimavicius et al. 2022 and Gene2Phenotype, a logistic regression probe achieves GOF one-vs-rest AUROC 0.73–0.80 under family-disjoint cross-validation, while dominant-negative and loss-of-function classes do not exceed AUROC 0.55 and 0.69 respectively. The selectivity to GOF survives cross-validation in which entire Pfam families are held out, suggesting ESM-2 has implicitly learned sequence properties characteristic of gain-of-function-prone proteins that generalise beyond protein family identity — without requiring task-specific fine-tuning.

### What's in v1

| Section | Content |
|---|---|
| Methods | ESM-2 650M, mean-pooled per-variant or per-gene embeddings, logistic regression + MLP probes, 5-fold gene-split AND family-split CV |
| Dataset | Gerasimavicius (948 genes) + merged with G2P/ClinVar pathogenic (1,985 genes total) |
| Headline table | Per-class AUROC (GOF / DN / LOF) under gene-split and family-split CV |
| Brief family-split justification | One paragraph explaining why family-disjoint CV matters (proteins cluster by family; family correlates with mechanism). No deep clustering analysis. |
| Contrast with ESMGain | One paragraph: ESMGain fine-tunes; we don't. Different and complementary claims. |
| Honest scope statement | "On Gerasimavicius + G2P, ESM-2 650M, frozen mean-pooled, gene-disjoint and family-disjoint CV" |

### What's deliberately NOT in v1

- Pathogenicity positive control
- Family-clustering quantification (k-NN purity, family probe, 74.8% within-family agreement)
- Pfam-coverage methodological cautionary tale
- Pathogenicity-mechanism dissociation framing
- Multi-seed replication
- Second model (SaProt / ESM-3)
- DDG2P replication
- Within-family analysis
- Stability-subspace projection
- Per-residue delta exploration

### Why v1 is safe to post

- Headline claim is **specific and bounded**: "frozen-probe GOF AUROC 0.73–0.80 family-split on this dataset with this probe"
- ESMGain is cited as the closest prior work and the contribution is framed as complementary, not competing
- Numbers reported with std across folds; nothing oversold
- No claim that ESM-2 "encodes mechanism" (broad); only that frozen embeddings retain cross-family GOF signal (narrow)
- Reproducibility: link to `esm2_mechanism/` scripts and JSON outputs

### Page target
**4–6 pages** (workshop-paper size). Abstract + intro + methods + one results table + one figure + brief discussion.

### The one figure
Bar chart: per-class AUROC (GOF / DN / LOF) × CV scheme (gene-split / family-split) × dataset (Gerasimavicius / merged). Six bars per panel, two panels by probe type (logreg / MLP). Makes the GOF survival vs DN+LOF collapse visually obvious.

### What v1 must NOT claim
- "ESM-2 encodes mechanism" (too broad)
- "We provide a novel methodology" (without the diagnostic sections, this is weak)
- "Our results contradict X" (we contradict no one explicitly at v1)
- "Family-split CV is the necessary diagnostic" (saved for v2)
- Anything about pathogenicity (saved for v2)

### Prior-work positioning (related work table for v1 intro)

The 2025–2026 GOF/LOF prediction literature is dense. The contribution sits in a specific gap none of these works occupy: **frozen mean-pooled ESM-2 embeddings + simple linear probe + family-disjoint CV on a large human disease gene set, with no fine-tuning, no extra modalities, and no architectural complexity.**

| Paper | Year | Key features | How it differs from this work |
|---|---|---|---|
| **ESM-Effect / ESMGain** (Glaser et al.) | 2025 | Fine-tunes ESM-2 (35M); emphasizes GOF with rBME metric | Strongly advocates fine-tuning over frozen embeddings; shows static/frozen embeddings underperform on their setup. Focus on DMS datasets, not large-scale human disease genes with family-disjoint CV. **Our finding (frozen embeddings retain GOF signal) is a direct complement: the signal is already there before fine-tuning, at least on disease-level mechanism labels.** |
| **Prediction of GOF/LOF/Neutral in Missense Variants** (Oliveira et al.) | 2025 | ESM-2 + classical ML (RF, XGBoost, LR) for GOF/LOF/Neutral | Uses embeddings but not frozen mean-pooled + LR only; reports lower GOF performance; does not run family-disjoint CV on a large disease gene set. |
| **ClearVariant** | 2025 | Attention-based model on ESM-2 for GOF/LOF | Requires a deeper learned architecture; our claim is that the signal is recoverable by a simple logistic regression on the frozen mean-pool — no architecture engineering required. |
| **PreMode** (Zhong et al.) | 2025 | Multimodal mode-of-action: ESM-2 + structure + MSA | Adds structure and MSA modalities. Our setup is intentionally minimal — frozen ESM-2 alone — which strengthens the interpretive claim about what ESM-2 representations themselves contain. |
| **ESMRank** | 2026 | ESM-2 embeddings + multimodal features for pathogenicity | Looks at GOF genes but does not evaluate variant-level GOF prediction with family-disjoint CV. |

**The specific frozen-mean-pool-plus-linear-probe-with-family-split setup is, to our knowledge, novel.** All four 2025 papers either fine-tune (ESMGain, ClearVariant), use heavier classical ML pipelines (Oliveira), or combine ESM-2 with additional modalities (PreMode, ESMRank). None of them isolate the question *"what mechanism information is already in frozen ESM-2 representations, without fine-tuning or auxiliary features, after blocking the family-recognition shortcut?"*

### What this means for v1's framing

v1 is positioned as **an interpretability/representation-analysis result**, not a competitive predictor. We are not claiming our setup beats ESMGain or ClearVariant on prediction accuracy; we are characterising what ESM-2's pretrained representation already contains. That framing is robust because:

- ESMGain shows fine-tuning helps → consistent with our finding that frozen signal exists but is modest
- Oliveira shows simpler ML pipelines underperform → consistent with the idea that current setups don't fully exploit what's in the embeddings
- PreMode shows adding modalities helps → consistent with frozen-ESM-2-only being a floor, not a ceiling
- Our contribution: characterise that floor precisely, with the leakage diagnostic that none of the above runs

The honest one-line positioning: *"While recent work focuses on building better mechanism predictors (via fine-tuning, additional modalities, or deeper architectures), we ask what mechanism information is already encoded in frozen ESM-2 representations, evaluated with family-disjoint cross-validation."*

---

## v2 — add controls and methodological framing (target: ~3 weeks after v1)

### What's added in v2

1. **Pathogenicity positive control** — same pipeline, ClinVar pathogenic vs benign, 17,236 variants, AUROC 0.88 family-split-stable. Establishes the pipeline can extract a signal when one is present, anchoring the modest mechanism numbers as a real ceiling rather than a pipeline failure.

2. **Pathogenicity–mechanism dissociation** as a second framing. The same embeddings encode whether a mutation is damaging strongly and linearly (AUROC 0.88) but how it acts only weakly and nonlinearly (macro-F1 floor ~0.39 family-split). This is a controlled side-by-side comparison that prior work (PreMode, AlphaMissense paper, LoGoFunc, Badonyi & Marsh) states qualitatively but does not demonstrate.

3. **Family-clustering quantification** — k=5 family purity 26× chance, 50-way family probe 27× majority baseline, 74.8% within-family mechanism agreement. Gives the family-split CV justification quantitative teeth.

4. **Pfam-coverage methodological note** — the worked example showing how silent CV failure (Δ=+0.011 inflated to +0.077 when annotations were extended) makes leakage diagnostics non-trivial to apply correctly. Useful cautionary tale.

5. **Multi-seed replication** — five seeds on all headline numbers, tighten confidence intervals.

### What v2's title becomes
*"Frozen ESM-2 embeddings encode pathogenicity broadly but mechanism narrowly: a controlled dissociation and family-split CV diagnostic"*

### Page target
**8–12 pages.** Adds one section on the positive control, one on the family clustering, one on the dissociation framing, and a discussion of methodology.

### v2's stronger claims (now defensible because of controls)
- "The mechanism null is a real absence of signal, not a pipeline failure (positive control AUROC 0.88)"
- "Family-split CV is necessary to detect family-recognition shortcuts; we provide quantified worked examples"
- "The pathogenicity-mechanism dissociation in PLM embeddings is real, controlled, and quantified"

---

## v3 — generalisation across models and datasets (target: ~2–3 months after v1)

### What's added in v3

1. **Second mechanism dataset (DDG2P)** — replicate the GOF-survives-family-split finding on an independent ~2,000-gene mechanism dataset curated by EBI G2P. If the pattern replicates, generalisation is established.

2. **Second model: SaProt first, ESM-3 second** — structure-aware PLMs as steelmen. Run in this order for principled reasons:
   - **SaProt** adds *structure tokens only* (foldseek 3Di). If SaProt recovers DN/LOF mechanism where ESM-2 fails, structure is the missing ingredient. If SaProt also fails, the negative result becomes much stronger.
   - **ESM-3** adds *structure + function tokens* (GO terms etc.). If ESM-3 recovers mechanism but SaProt didn't, the cause is function-token supervision, not structure. If both fail, the negative claim generalises across pretraining objectives.
   - **Why SaProt before ESM-3:** SaProt is fully open-weight (ESM-3 is partially gated through EvolutionaryScale), cheaper per embedding, and cleanly isolates *structure* as the variable. ESM-3 alone would conflate two pretraining differences (structure + function) and not tell us which mattered.
   - **Why both, not just one:** SaProt + ESM-3 together let us answer "which pretraining ingredient (if any) recovers mechanism information," which is more decisive than either alone.

3. **Within-family mechanism analysis** — test whether mechanism is learnable inside a single Pfam family. The potential positive flip side: if mechanism IS recoverable within a homologous family, the field has been measuring the wrong problem. Reframes the paper from "PLMs don't do mechanism cross-proteome" to "mechanism prediction is a within-family problem."

4. **Evo2 comparison** — if relevant. Tests whether genomic-context features (paralogs, dosage) recover signal that pure-protein features miss.

5. **Per-class PR-AUC and calibration** — for clinical relevance.

### What v3's title becomes
Depends on outcomes:
- If SaProt + ESM-3 + DDG2P all confirm the GOF-selective pattern: *"Cross-family disease mechanism is selectively encoded for gain-of-function in protein language models: a multi-model, multi-dataset characterisation"*
- If within-family analysis is positive: *"Disease mechanism prediction from protein language models is a within-family problem"*
- If SaProt recovers mechanism but ESM-2 doesn't: *"Structure-aware protein language models recover disease mechanism that sequence-only models miss"*
- If ESM-3 recovers mechanism but SaProt doesn't: *"Function-aware pretraining is necessary for mechanism encoding in protein language models"*

### Page target
**15–25 pages.** Full peer-reviewed paper. Target *Bioinformatics* or *Genome Biology* methodological note; possibly *Nat Methods* if within-family or SaProt results are strong.

---

## Versioning strategy notes

### Why v1 must be narrower than what you eventually want to claim

bioRxiv versions are a public scientific record. If v1 says X and v3 walks back to Y, readers cite v1's X. Worst case: v1 gets cited as a retraction. Best practice: v1 scope is a strict subset of v3 scope, so each version *strengthens* rather than *replaces* the prior.

### What stays the same across all versions

- The GOF AUROC 0.73–0.80 family-split number
- The DN/LOF chance-level family-split numbers
- The contrast with ESMGain (frozen vs fine-tuned)
- The methods (mean-pooled ESM-2 650M, logreg + MLP, gene-split + family-split CV)

These are the load-bearing claims. They're already supported by the data in `results/20260524_baseline_run/run_0/`.

### What might change across versions

- Multi-seed replication might shift point estimates ±0.02–0.04 — flag this in v1 by reporting fold std
- Re-running MLP delta on merged dataset with corrected Pfam might shift the MLP numbers (currently pending) — v1 should report Gerasimavicius MLP results only
- DDG2P / SaProt could falsify the GOF claim — v3 will report honestly either way

### What gets cut entirely if results don't cooperate

- If SaProt + DDG2P don't replicate the GOF survival, the paper ends at v2 (with a discussion noting the boundary)
- If within-family analysis is positive but cross-family isn't, v3 reframes around within-family
- If a 2026 paper preempts the frozen-probe GOF finding before v1 posts, the paper ends at v1 as a "we also observed this" preliminary note

---

## Practical timeline

| Stage | Time | Blocker |
|---|---|---|
| v1 writing | 3–5 days | None — all data in hand |
| v1 figure | 1 day | None |
| v1 polish + post | 1 day | None |
| **v1 live** | **~1 week from start** | |
| v2 — pathogenicity & clustering writeup | 5–7 days | Multi-seed reruns (cheap on cached embeddings) |
| v2 figure additions | 1 day | None |
| **v2 live** | **~3 weeks from start** | |
| v3 — DDG2P embedding extraction | ~3 days GPU | RunPod availability |
| v3 — SaProt embedding extraction | ~3 days GPU | RunPod availability (run before ESM-3 — cleaner control) |
| v3 — ESM-3 embedding extraction | ~5–7 days GPU | RunPod availability + EvolutionaryScale access for larger ESM-3 variants |
| v3 — within-family analysis | 1 day | After embeddings |
| v3 writing & figures | 2–3 weeks | Above experiments must complete |
| **v3 live** | **~2–3 months from start** | Tighter timeline if ESM-3 deferred to v4 |

---

## Files needed before v1 can post

| File | Status |
|---|---|
| Headline numbers JSON | ✓ `results/20260524_baseline_run/run_0/option_b_gene_level_wt_merged.json` |
| MLP probe results | ✓ `results/20260524_baseline_run/run_0/mlp_results_seed0.json` |
| Per-class AUROCs all conditions | ✓ in JSONs above |
| Multi-seed for v1 numbers | ✗ Need to run 5 seeds — cheap |
| Figure | ✗ Need to write `plot_publication_v1.py` |
| LaTeX template | ✗ Pick one (`bioRxiv-style.cls` works) |
| ESMGain citation context | ✗ Read ESMGain paper carefully and confirm the framing contrast |

---

## Single-line summary

**v1 = frozen-probe GOF cross-family signal (narrow, fast, safe). v2 = + positive control + leakage diagnostic. v3 = + generalisation across models and datasets. Each version is a publishable scientific record; later versions strengthen, never replace.**
