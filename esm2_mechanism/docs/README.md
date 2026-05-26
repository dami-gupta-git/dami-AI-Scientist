# esm2_mechanism — results index

**This is a standalone research project, not an AI Scientist run.** All experiments were designed and executed manually. The code lives in `esm2_mechanism/scripts/` and was run directly on RunPod (A100 80GB) or local CPU. The project will likely move to its own repository.

Eleven `result_*.md` files written across May 23–25, 2026. Read in the order below for the coherent narrative arc. Results 3 and 5 are superseded by result 7.

---

## Current state (as of result 11 pilot)

Two threads are live:

1. **Results 1–10 (frozen ESM-2 + supervised probes):** mechanism floor under family-split CV is ~0.36–0.40 macro-F1 with the contrastive probe (result 9), and clan-holdout shows ~half the family-split signal is fold memorisation (result 10). Pathogenicity positive control AUROC 0.88 (result 6) confirms the pipeline is sound. DN remains the hardest class throughout.
2. **Experiment 11 (gene-level proteome features, no ESM-2):** Stage 0 pilot (result 11) on 4 gene-level features (gnomAD pLI, LOEUF, mis_z + Ensembl paralog count) achieves macro-F1 = **0.417 ± 0.009** family-split over 5 seeds, all 5/5 returning STRONG_SIGNAL by the pre-registered rule. Pilot is a sanity check only; Phase 1 with ~30 proteome features is the actual experiment (sibling agent currently handling).

Comparisons between the two threads — and to prior gene-level mechanism predictors (Badonyi & Marsh 2024) — are deliberately deferred to Phase 1, when feature sets are comparable. See `plan_experiment.md` for the staged plan.

---

## Reading order

### 1. `result_1.md` — Baseline: linear probe sets up the puzzle
**Script:** `experiment.py` · **Run:** May 23–24, gene-split CV on 10,231 variants (948 genes)
**Headline numbers:** Linear delta_mean macro-F1 = **0.279** (chance). WT-only mysteriously = **0.580**. Per-residue delta (0.373) beats mean-pooled.
**What it concludes:** Linear delta probe finds no mechanism signal; WT-only's 0.58 is the suspicious number that needs explaining.
**Open question after reading:** *Is WT-only's 0.58 a real gene-level mechanism signal, or homology leakage from gene-split CV?*

### 2. `result_2.md` — Gene-split vs family-split baselines
**Script:** `family_split_baselines.py` · Same embeddings as result 1, 8 feature sets × 2 CV schemes
**Headline numbers:** WT-only macro-F1 **collapses** 0.580 → 0.389 under family-split (Δ = +0.191). Delta probe stays flat. GOF AUROC = **0.801** under family-split WT-only.
**What it concludes:** Most of WT-only signal is paralog leakage. AlphaMissense carries zero mechanism information.
**Open question after reading:** *Is there any nonlinear mechanism signal in the delta the linear probe couldn't see?*

### 3. `result_3.md` — MLP probe on delta (FIRST pass) ⚠ SUPERSEDED BY RESULT 7
**Script:** `experiment_mlp.py` · MLP probe (256→64, dropout 0.3) on delta_mean, gene-split only
**Headline numbers:** MLP delta_mean macro-F1 = **0.414**.
**⚠ Why this is superseded:** No family-split CV — the 0.414 is gene-split only and includes family leakage. Result 7 provides the correct family-split number (0.364) and calibration against chance.

### 4. `result_4.md` — Family clustering: the causal explanation
**Script:** `family_clustering.py` · Pfam clustering analysis on WT, mut, delta embeddings
**Headline numbers:** k=5 family purity = **26× chance** (z = +78). 50-way family probe = **27× majority baseline**. **74.8%** of genes share their family's majority mechanism.
**What it concludes:** WT-only signal explained by family recognition × family-mechanism correlation. Includes novelty assessment (2/5 — folk wisdom, but not yet demonstrated as a controlled comparison).

### 5. `result_5.md` — Nonlinear probes (MLP/kNN/GBM/RF) ⚠ PARTIALLY SUPERSEDED BY RESULT 7
**Script:** `experiment_mlp.py` extended · 4 probes on delta_mean and delta_pos, gene-split only
**Headline numbers:** MLP = 0.431, kNN = 0.410, GBM = 0.336, RF = 0.292.
**⚠ Limitation:** Gene-split only. Result 7 provides the family-split calibration showing 62% of the gene-split lift is leakage.

### 6. `result_6.md` — Pathogenicity positive control
**Script:** `pathogenicity_control.py` · 17,236 ClinVar pathogenic/benign variants, 944 genes
**Headline numbers:** Pathogenicity MLP AUROC = **0.878**, family-split Δ = **0.002**.
**What it concludes:** Pipeline is sound. Pathogenicity AUROC 0.88 (family-split-stable) vs mechanism floor ~0.39 (family-split) — **the dissociation is sharper than originally framed** (see result 7 for correction).

### 7. `result_7.md` — Full calibration: all numbers, honest framing
**Scripts:** `experiment_mlp.py` with family-split, `build_merged_dataset.py`, Option B gene-level WT
**Headline numbers:**
- MLP delta gene-split **0.415**, family-split **0.364** (+0.031 above chance; 62% of lift is leakage)
- Gene-level WT merged dataset family-split **0.393**, GOF AUROC **0.728**
- Always-predict-LOF baseline: **0.279** (Gerasimavicius), **0.311** (gene-level merged)
- Family-split floor **~0.39** consistent across 3 different setups (per-variant/gene-level, 2 datasets, linear/MLP)
**What it concludes:** The ~0.39 floor is real but small. The pathogenicity-mechanism dissociation is sharper than result 6 suggested. The GOF AUROC (0.73–0.80 family-split) is the strongest individual signal. See PUBLISH.md for v1 paper plan built around this.

### 8. `result_8.md` — Within-family mechanism (first pass)
**Script:** ad-hoc analysis on cached Gerasimavicius embeddings · Local CPU, seed=42
**Headline numbers:** Within-family gene-split CV on the 5 largest Pfam families. **PF00520 (ion channel) delta F1=0.407, AUROC=0.659 (2-class GOF/DN)** — the most interpretable result. Other families largely at chance due to tiny sample sizes (6–12 genes).
**What it concludes:** Directional signal that mechanism is partially learnable within a homologous subfamily; consistent with MissION-style findings. Not publishable at single seed + small N; merits replication on merged dataset with MLP.

### 9. `result_9.md` — Contrastive metric learning recovers cross-family signal
**Script:** `contrastive_mechanism.py` · A100 80GB, seed=0
**Headline numbers:** Supervised contrastive projection head (1280→256→64, TripletMarginLoss, positives = same mechanism / different family) pushes family-split macro-F1 from MLP's 0.364 to **0.397** on Gerasimavicius (+0.033) and to **0.387** on merged (+0.035 above MLP floor). Lift is equal under gene-split (+0.060) and family-split (+0.059) — the critical diagnostic that the recovered signal is *not* leakage. LOF benefits most; **DN stays flat (+0.012 Geras, −0.025 merged)**.
**What it concludes:** Frozen ESM-2 delta does encode small cross-family mechanism signal not accessible to a standard MLP — but only for LOF. DN remains essentially absent.

### 10. `result_10.md` — Clan-level holdout: partial generalisation, not pure memorisation
**Script:** `clan_holdout.py` · Local CPU, seed=0
**Headline numbers:** Leave-one-Pfam-clan-out evaluation across 21 qualifying clans gives MLP macro-F1 = **0.299 ± 0.076** — below family-split floor (0.352) but above majority (0.254). Per-class AUROCs (GOF 0.597 / DN 0.575 / LOF 0.636) confirm real cross-fold signal. **Approximately half the family-split mechanism signal is clan-level memorisation; the remainder is genuine cross-fold generalisation.** Heterogeneous across clans (Cupin F1=0.536; Ion_channel F1=0.190).
**What it concludes:** The ~0.36 family-split floor is roughly half real, half fold-memorisation. Mechanism is more readable from sequence in architectures with stereotyped structural mechanisms (cupins, death domains, GPCRs) and unreadable in plastic repeat proteins (ankyrins, EF-hands).

### 11. `result_11.md` — Stage 0 pilot: 4 gene-level features predict mechanism under family-split CV
**Script:** `proteome_pilot.py` · Local CPU, seeds 0–4 (5-seed replication)
**Headline numbers:** Logistic regression on 4 public gene-level features (pLI, LOEUF, mis_z, paralog_count) under family-split CV on 1,234 genes / 725 families achieves macro-F1 = **0.417 ± 0.009** (+0.122 above majority 0.295). Per-class AUROCs (mean ± std): GOF **0.686 ± 0.011**, **DN 0.687 ± 0.009**, LOF **0.735 ± 0.001** — balanced and tight across seeds.
**What it concludes:** Stage 0 sanity check passes. Public gene-level features carry meaningful mechanism signal under the project's family-split CV, robust to seed choice. 5/5 seeds returned STRONG_SIGNAL by the pre-registered rule. Proceed to Phase 1 (full ~30-feature pull). Comparisons to prior gene-level classifiers and to ESM-2 results 1–10 are deferred to Phase 1 when the feature set is comparable.

---

## The coherent story across all 11

1. **(1)** Linear probe at chance on delta; WT-only at 0.58 — suspicious.
2. **(2)** Family-split: WT-only collapses (0.58→0.39). Delta stays flat. GOF AUROC 0.80 survives.
3. **(4)** Family clustering explains the collapse: ESM-2 encodes family identity, family correlates with mechanism.
4. **(6)** Positive control: pathogenicity AUROC 0.88, family-split-stable. Pipeline works when signal is there.
5. **(7)** Full calibration: ~0.39 family-split floor across all setups. 62% of gene-split signal is leakage. Dissociation with pathogenicity is **sharper** than originally thought.
6. **(8)** Within ion-channel family, mechanism is learnable (PF00520 GOF/DN AUROC 0.659) — directional, small N.
7. **(9)** Contrastive metric learning lifts cross-family floor from 0.364 → 0.397 with equal gene-split / family-split deltas. LOF benefits; DN does not.
8. **(10)** Clan holdout: ~half the family-split signal is fold memorisation, ~half is genuine. Heterogeneous across architectures.
9. **(11)** Stage 0 pilot for the proteome-features thread. Four gene-level features (no ESM-2) reach macro-F1 0.417 ± 0.009 and DN AUROC 0.687 ± 0.009 over 5 seeds. Sanity check only; the real comparison lives in Phase 1 of Experiment 11.

**Publication plan:** `PUBLISH.md` — v1 plan still based on results 1–10 (frozen-pLM mechanism characterisation). Reframing decisions wait on Phase 1 of Experiment 11.

---

## Supporting docs

- `EXPERIMENT.md` — Pre-registration document (original hypothesis and thresholds)
- `PUBLISH.md` — Publication plan: v1/v2/v3 versioned bioRxiv strategy
- `plan_experiment.md` — **Experiment 11 plan: per-variant ESM-2 + gene-level proteome features** (staged execution: pilot → V2 → V3 → V4, with pre-registered decision rules and family-mean-centering baked in)
- `explain.txt` — Plain-English explanation of the experiment design
- `progress_notes.md` — Running log of decisions, bugs fixed, observations
- `../scripts/README.md` — What each script does and when to use it

---

## Companion data

| Result | JSON file |
|---|---|
| 1 | `results/20260524_baseline_run/run_0/final_info_seed0.json` |
| 2 | `results/20260524_baseline_run/run_0/family_split_baselines.json` |
| 3, 5 | `results/20260524_baseline_run/run_0/mlp_results_seed0.json` |
| 4 | `results/20260524_baseline_run/run_0/family_clustering.json` |
| 6 | `results/20260524_baseline_run/run_0/pathogenicity_control.json` |
| 7 | `results/20260524_baseline_run/run_0/option_b_gene_level_wt_merged.json` + merged dataset MLP |
| 8 | `results/20260524_baseline_run/run_0/within_family_analysis.json` |
| 9 | `results/20260524_baseline_run/run_0/contrastive_results_{geras,merged}_seed0.json` |
| 10 | `results/20260524_baseline_run/run_0/clan_holdout_results_seed0.json` |
| 11 | `results/proteome_pilot/pilot_results.json` |

Embeddings under `data/embeddings/`:
- `embeddings_{wt,mut}{,_pos}_esm2_t33_650M_UR50D.npy` — Gerasimavicius (results 1–5, 7)
- `merged_embeddings_{wt,mut}_{mean,pos}.npy` — merged 1,985-gene dataset (results 7, 9)
- `emb_{wt,mut}_mean_pathogenicity_*.npy` — ClinVar pathogenicity set (result 6)

Proteome features (result 11 pilot):
- `data/gene_features_pilot.tsv` — per-gene feature table (4 features)
- `data/cache/proteome_pilot/gnomad_v4.1_constraint.tsv` — cached gnomAD bulk constraint
- `data/cache/proteome_pilot/paralogs/` — per-gene Ensembl REST cache

---

## Highest-priority next experiments

1. **Experiment 11 Phase 1** — full proteome feature pull (HPA tissue specificity, PaxDb abundance, BioPlex interactome, Mathieson half-life, PhosphoSitePlus PTM density, ClinGen dosage sensitivity). See `plan_experiment.md`. Pilot has already returned STRONG_SIGNAL with just 4 features; the question is how much further V2 pushes, especially DN.
2. **Multi-seed replication** — pilot is seed=0 only; all earlier results are seed=0 or seed=42 only. Five seeds blocking before any posting.
3. **The figure** — bar chart of per-class AUROC × CV scheme × probe/feature class, now updated to include result_11's gene-level result side-by-side with ESM-2 results 7–10.
4. **Within-family MLP on merged dataset** — completes result 8's directional finding with the larger gene set and MLP probe. ~$5, an afternoon.
5. **End-to-end LoRA contrastive fine-tune of ESM-2** — only worth doing if Phase 1 of Experiment 11 shows that ESM-2 + proteome (V3) genuinely beats proteome alone (V2). Otherwise sequence is dispensable.
