# Scripts

Listed by project arc. See `docs/README.md` for the results these scripts produced.

## Arc 1 — Frozen ESM-2 characterisation (results 1–10)

**`experiment.py`**
Full baseline pipeline. Downloads Gerasimavicius dataset, fetches UniProt sequences and Pfam families, extracts ESM-2 650M embeddings, fits stability subspace, runs linear probe with gene-split and family-split CV, baselines, and probe direction orthogonality. Writes `run_0/final_info.json`. Requires GPU for embedding extraction. *Results 1–2, 4, 7.*

**`experiment_mlp.py`**
Nonlinear probes on cached delta embeddings. PyTorch MLP (1280→256→64→3) with dropout, class weighting, and gene-split early stopping. Also runs GBM, RF (PCA-50), and kNN probes. Use `--family_split` for family-split CV. CPU-only (reads cached .npy files). Produces `mlp_results_seed{N}.json`. *Results 3, 5, 7.*

**`family_clustering.py`**
Diagnostic: quantifies how strongly ESM-2 embeddings cluster by Pfam family. Computes silhouette, k-NN family purity (k=5,10) vs shuffled null, within/between cosine distance ratio, and a linear family-probe accuracy on gene-level embeddings. Produces `family_clustering.json`. *Result 4.*

**`family_split_baselines.py`**
All baselines (WT-only, mutant-only, WT+mutant concat, delta mean, delta per-residue, one-hot AA, FoldX ΔΔG, AlphaMissense) under both gene-split and family-split CV. Prints Δ(gene−family) per feature — positive values indicate homology leakage. Produces `family_split_baselines.json`. *Result 2.*

**`pathogenicity_control.py`**
Positive control: predicts ClinVar pathogenic vs benign on the same gene set using the same pipeline. Three phases: fetch ClinVar variants → extract ESM-2 embeddings (GPU) → run probes under gene-split and family-split. Validates pipeline soundness. Produces `pathogenicity_control.json`. *Result 6.*

**`contrastive_mechanism.py`**
Supervised contrastive projection head (1280→256→64) with TripletMarginLoss and cross-family-only positive sampling. Evaluates via k-NN (k=10, cosine) in projected space. Requires GPU for sensible runtime. *Result 9.*

**`clan_holdout.py`**
Leave-one-Pfam-clan-out CV on cached embeddings. Uses Pfam-A.clans.tsv.gz to map families → clans. CPU-only. *Result 10.*

**`mut_only_mlp.py`**
MLP on mut-only embeddings (no WT subtraction). Demonstrates that mut-only ≈ WT-only — confirms WT-only signal is gene identity, not mutation-specific. *Supporting analysis for result_2's WT-only interpretation.*

## Arc 2 — Gene-level proteome features (results 11–14)

**`proteome_pilot.py`**
Stage 0 pilot: downloads gnomAD v4.1 constraint TSV, fetches Ensembl paralog counts (resume-safe, rate-limited), trains LogReg + tiny MLP on 4 features (pLI, LOEUF, mis_z, paralog_count) under family-split CV. 5-seed pre-registered decision rule. Produces `pilot_results_summary_5seed.json`. *Result 11.*

**`build_proteome_features.py`**
Builds the 37-column feature matrix from 6 sources: gnomAD constraint, Ensembl paralogs, HPA tissue specificity (categorical mapped to 5-value score), PaxDb abundance, BioPlex PPI degree, ClinGen HI/TS. Computes family-mean-centred residuals and binary missingness indicators. Outputs `data/gene_proteome_features.tsv`, `data/proteome_features_aligned.npy`, `data/proteome_feature_columns.json`. *Result 12.*

**`proteome_mechanism.py`**
V1–V4 modelling: ESM-2 delta MLP (V1), proteome LogReg+MLP (V2), concat MLP (V3), contrastive head on V3 inputs (V4). 5-fold family-split CV, 5 seeds. Produces per-seed JSONs + summary in `results/proteome_mechanism/`. *Result 13 main matrix.*

**`per_gene_ablation.py`**
Two analyses: (a) T2 per-gene scoring by aggregating per-variant predictions to per-gene; (b) T4 feature-class ablation on V2 (drop one class at a time, report ΔF1 and Δ per-class AUROC). Reuses cached V1/V2/V3 models from `proteome_mechanism.py`. *Result 13 addenda.*

**`clinical_utility.py`**
Within ClinGen HI=3 evaluation: LogReg + MLP trained under family-split CV across 5 seeds, two feature sets (FULL 37-dim vs NO-MISS 18-dim, dropping missingness indicators). Reports GOF-vs-LOF and DN-vs-LOF AUROC within HI=3 subset, plus single-feature baselines (pLI, LOEUF, mis_z, paralog_count, PPI_degree) for context. Bootstrap CIs and calibration on seed 0. *Result 14.*

## Arc 3 — Badonyi structural priors + within-family (results 15–16)

**`build_badonyi_features.py`**
Loads Badonyi 2024 PLOS One S3 Table (pDN/pGOF/pLOF for 20,365 human proteins), aligns to merged_gene_list, computes residuals and missingness indicators. Outputs `data/badonyi_features.tsv`, `data/badonyi_features_aligned.npy`. *Result 15 data prep.*

**`badonyi_mechanism.py`**
Six-variant evaluation under family-split: V1 (ESM-2), V2 (proteome), V_bad (Badonyi pDN/pGOF/pLOF only), V2+bad, V1+bad, V_all. LogReg or MLP per variant. 5 seeds. Produces `results/badonyi_mechanism/badonyi_mechanism_summary.json`. *Result 15.*

**`badonyi_leakage_analysis.py`**
Leakage triage for V_bad: stratifies family-split CV by Badonyi-training-set membership (IN vs OUT vs ALL). Compares V2/V_bad/V2+bad performance across regimes. *Result 15 Appendix A.*

**`fetch_uniprot_sequences.py`**
Fetches missing UniProt sequences (1,035 entries not in `data/sequences.json`) via UniProt REST in batches of 100. Resume-safe. *Result 15 Appendix B prerequisite.*

**`mmseqs_cluster_holdout.py`**
MMseqs2 cluster-split CV: holds out entire sequence-similarity clusters at 20% identity, 20% coverage (matching Saadat & Fellay 2025). Re-runs V1/V2/V_bad/V2+bad/V_all. Requires MMseqs2 installed (`brew install mmseqs2`) and the cluster mapping at `data/mmseqs_clusters.json`. *Result 15 Appendix B.*

**`within_family_mechanism.py`**
Leave-one-gene-out CV within each Pfam family that qualifies (≥6 genes, ≥2 mechanism classes). Compares 5 feature sets: raw proteome, residual proteome (family-mean-centred), raw Badonyi, residual Badonyi, combined residuals. 24 families, 238 genes. *Result 16.*

**`badonyi_holdout_survival.py`**
Tests whether Badonyi's *raw published* SVM predictions (pDN/pGOF/pLOF) survive family-split and MMseqs2-20 cluster-split. No re-training — uses the published predictions directly. Pre-registered decision rules from `docs/plan_badonyi.md`. Also stratifies IN vs OUT of Badonyi's training set. *Result 16 addendum.*

## Data utilities

**`fetch_clinvar_variants.py`**
ClinVar pathogenic/likely-pathogenic missense variants for a gene list. Resume-safe, rate-limited to ≤3 NCBI req/s. Outputs `clinvar_variants.tsv`. *Used to extend the merged dataset beyond Gerasimavicius.*

**`build_merged_dataset.py`**
Builds the merged 1,985-gene / 19,100-variant dataset from Gerasimavicius + G2P + ClinVar. *Result 7.*

**`extract_merged_embeddings.py`**
Extracts ESM-2 650M WT and mutant embeddings for the merged variant list. Requires GPU. Outputs `data/embeddings/merged_embeddings_*.npy`.

## Visualisation

**`plot.py`**
Matplotlib plots from `final_info.json`: AUROC bar chart per class, probe direction cosine matrix heatmap, variance-explained bar chart. Run as `python plot.py run_0`. *Result 1 era; needs updating to cover results 11–16 if used for the v2 figure.*
