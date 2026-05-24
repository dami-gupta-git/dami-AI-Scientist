# esm2_mechanism — results index

Six `result_*.md` files written across May 23–24, 2026. Read in the order below for the coherent narrative arc. Result 3 is superseded by results 4–6 and should be read with that caveat.

---

## Reading order

### 1. `result_1.md` — Baseline: linear probe sets up the puzzle
**Script:** `experiment.py` · **Run:** May 23–24, gene-split CV on 10,231 variants (948 genes)
**Headline numbers:** Linear delta_mean macro-F1 = **0.279** (chance). WT-only mysteriously = **0.580**. Per-residue delta (0.373) beats mean-pooled.
**What it concludes:** Linear delta probe finds no mechanism signal; WT-only's 0.58 is the suspicious number that needs explaining.
**Open question after reading:** *Is WT-only's 0.58 a real gene-level mechanism signal, or homology leakage from gene-split CV?*

### 2. `result_2.md` — Gene-split vs family-split baselines
**Script:** `family_split_baselines.py` · Same embeddings as result 1, 8 feature sets × 2 CV schemes
**Headline numbers:** WT-only macro-F1 **collapses** 0.580 → 0.389 under family-split (Δ = +0.191). Delta probe stays flat (no signal to lose). Surprising survivor: GOF AUROC = **0.801** under family-split WT-only.
**What it concludes:** Most of the WT-only "mechanism signal" was paralog leakage. AlphaMissense carries zero mechanism information. The one thread worth pulling: cross-family GOF signal in WT.
**Open question after reading:** *Is there any nonlinear mechanism signal in the delta the linear probe couldn't see?*

### 3. `result_3.md` — MLP probe on delta (FIRST pass) ⚠ SUPERSEDED
**Script:** `experiment_mlp.py` · MLP probe (256→64, dropout 0.3) on delta_mean
**Headline numbers:** MLP delta_mean macro-F1 = **0.414** (up from linear 0.279). GOF AUROC 0.729, LOF AUROC 0.727 — both cross pre-registered 0.72 "meaningful" threshold.
**What it concludes (at the time):** Mechanism is nonlinearly encoded in delta; revises the original claim accordingly.
**⚠ Why this is superseded:** Result 5 reframes the same experiment as having two competing explanations (real signal vs residual family leakage) and lands on family leakage as more parsimonious. Result 4 quantifies the family-leakage mechanism. Result 6 confirms via positive control. The "mechanism is nonlinearly encoded" conclusion here is no longer the favoured interpretation.

### 4. `result_4.md` — Family clustering: the causal explanation
**Script:** `family_clustering.py` · Pfam clustering analysis on WT, mut, delta embeddings
**Headline numbers:** k=5 family purity = **26× chance** (z = +78). 50-way family probe accuracy = **0.587** (27× majority baseline). **74.8%** of disease genes share their family's majority mechanism. Delta embeddings strip most family signal but retain a residual (z = +18).
**What it concludes:** The 0.58 WT-only signal is fully explained by family recognition × family-mechanism correlation. No mechanism learning required. Includes a novelty assessment (2/5 — folk wisdom).

### 5. `result_5.md` — Nonlinear probes (MLP/kNN/GBM/RF), CAUTIOUS revisit
**Script:** `experiment_mlp.py` extended · 4 nonlinear probes on delta_mean and delta_pos
**Headline numbers:** MLP delta_mean macro-F1 = 0.431, kNN = 0.410, GBM = 0.336, RF = 0.292. delta_pos shows no nonlinear lift. DN stays at chance across all probes (AUROC ~0.52–0.57).
**What it concludes:** MLP lift is real but its cause is unresolved — could be nonlinear mechanism signal (A) OR nonlinear recovery of residual family clustering (B). Lands on B as more parsimonious and identifies the resolving experiment: **MLP under family-split CV** (still pending).
**Relationship to result 3:** Same experiment, more probes, more honest framing. Read this instead of result 3 if you only have time for one.

### 6. `result_6.md` — Pathogenicity positive control: the final answer
**Script:** `pathogenicity_control.py` · 17,236 ClinVar pathogenic/benign variants across 944 genes / 658 Pfam families
**Headline numbers:** delta_mean MLP AUROC = **0.878** for pathogenicity. Gene-split → family-split drop = **0.002** (essentially zero). WT-only barely above chance for pathogenicity (AUROC ~0.54–0.60).
**What it concludes:** Pipeline is sound (positive control passes). The asymmetry IS the central finding: **ESM-2 encodes whether a mutation matters, not how it acts.** Includes calibrated novelty assessment showing this is folk wisdom that has been formally motivated (PreMode, AlphaMissense paper, LoGoFunc, Badonyi & Marsh) but never demonstrated as a controlled side-by-side comparison.

---

## The coherent story across all 6

1. **(1)** Linear probe says delta has no mechanism signal; WT-only mysteriously has some.
2. **(2)** Family-split CV says WT-only's signal mostly evaporates — family leakage suspected.
3. **(3 / 5)** Nonlinear probes lift delta from 0.28 → 0.42 — but this is more likely residual family signal than mechanism learning.
4. **(4)** Family clustering quantified — 26× chance, 75% within-family mechanism agreement — provides the causal explanation.
5. **(6)** Positive control: same pipeline gets AUROC 0.88 on pathogenicity, family-split-stable. **Asymmetry is the publishable finding**: pathogenicity yes, mechanism no.

---

## Supporting docs (not part of the result series)

- `EXPERIMENT.md` — Pre-registration document (hypothesis, predictions, pre-registered thresholds)
- `explain.txt` — Plain-English explanation of the experiment design
- `progress_notes.md` — Running log of decisions and observations across the project
- `../scripts/README.md` — What each script does and when to use it

---

## Companion data

Each result references metric files under `../results/20260524_baseline_run/run_0/`:

| Result | JSON file |
|---|---|
| 1 | `final_info_seed0.json`, `detailed_results_seed0.json` |
| 2 | `family_split_baselines.json` |
| 3, 5 | `mlp_results_seed0.json`, `mlp_probe_results.json` |
| 4 | `family_clustering.json` |
| 6 | `pathogenicity_control.json` (also at `../data/pathogenicity_control.json`) |

Embeddings under `../data/embeddings/`:
- `embeddings_{wt,mut}{,_pos}_esm2_t33_650M_UR50D.npy` — Gerasimavicius (used by results 1–5)
- `emb_{wt,mut}_mean_pathogenicity_esm2_t33_650M_UR50D_n17259.npy` — ClinVar pathogenicity set (used by result 6)

---

## Publishability summary (from result 6)

- **bioRxiv preprint:** publishable as is, framed as methodological consolidation (controlled demonstration + leakage diagnostic + MissION reconciliation), not as a discovery paper.
- ***Bioinformatics* / *Genome Biology* methodological note:** realistic with current evidence.
- ***Nat Methods* / *Nat Commun*:** not realistic from this evidence alone — PreMode and AlphaMissense are too directly adjacent. Would require adding the within-family mechanism analysis (positive flip side) and structure-aware model replication (SaProt / ESM-3).

---

## Highest-priority next experiments (still open)

1. **MLP under family-split CV** on WT-only AND delta — directly resolves the explanation A vs B question in result 5
2. **DDG2P replication** — second dataset (~2,000 genes)
3. **SaProt or ESM-3 replication** — the structure-aware steelman
4. **Within-family mechanism analysis** — test whether mechanism is learnable inside a single Pfam family (potential positive flip side)
