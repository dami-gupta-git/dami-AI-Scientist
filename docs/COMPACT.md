# Session Compact — May 21-22, 2026

## Who You Are Working With
Dami Gupta — computational biology researcher exploring whether genomic/protein foundation model embeddings encode functional information. Working directory: `/Users/dgupta/code/downloads/AIScientist/AI-Scientist`. Prior project: EvoMine (`/Users/dgupta/code/portfolio/EvoMine/`). Prefers concise responses, no emojis, no trailing summaries. Does not want layer ablation experiments. Wants to use AI Scientist framework to generate and run hypotheses automatically.

---

## What This Project Is

We are using the **AI Scientist** framework (github.com/SakanaAI/AI-Scientist, forked to github.com/dami-gupta-git/dami-AI-Scientist) to automatically generate research hypotheses, implement them, run experiments, and write papers. The framework uses:
- **Aider** to modify experiment code based on LLM-generated ideas
- **Claude** (claude-sonnet-4-5 — the only model that works with Dami's API key) to generate ideas and write papers
- **RunPod** GPU (L40S 46GB) to run experiments

The core research question: **do protein/DNA foundation model embeddings encode functional gene class information, and if so, how?**

---

## Repo State

- **Fork:** `https://github.com/dami-gupta-git/dami-AI-Scientist.git`
- **Local working directory:** `/Users/dgupta/code/downloads/AIScientist/AI-Scientist`
- **Current local branch:** `esm2-depmap`
- **RunPod directory:** `/evo2/dami-AI-Scientist`
- **RunPod branch:** `esm2-depmap`
- `results/` is gitignored — results pulled manually
- Model must be `claude-sonnet-4-5` (older claude-3-5-sonnet models 404 on Dami's API key)

### Branches
| Branch | Purpose |
|---|---|
| `evo2-template` | Evo2 DNA embeddings — DNA repair vs TSG classification |
| `evo2-validation` | Evo2 DNA embeddings — glycolysis vs T-cell validation |
| `esm2-function` | ESM-2 protein embeddings — DNA repair vs TSG (same gene set) |
| `esm2-depmap` | **Currently running** — ESM-2 vs DepMap Mantel test |

---

## RunPod SSH Access

```bash
# Write key to temp file first:
cat > /tmp/runpod_key << 'EOF'
-----BEGIN OPENSSH PRIVATE KEY-----
b3BlbnNzaC1rZXktdjEAAAAABG5vbmUAAAAEbm9uZQAAAAAAAAABAAAAMwAAAAtzc2gtZW
QyNTUxOQAAACCoiNs7CclxszcKTUjWu6ZUowMMUyhCBQJIWS7gR4/dowAAAKj2aahM9mmo
TAAAAAtzc2gtZWQyNTUxOQAAACCoiNs7CclxszcKTUjWu6ZUowMMUyhCBQJIWS7gR4/dow
AAAECuj0NFcVrvrIYO+oGUdIR9AXX3X8VkRIf2CrtmkM6SfKiI2zsJyXGzNwpNSNa7plSj
AwxTKEIFAkhZLuBHj92jAAAAHmRndXB0YUBNYWMuaHNkMS5tYS5jb21jYXN0Lm5ldAECAw
QFBgc=
-----END OPENSSH PRIVATE KEY-----
EOF
chmod 600 /tmp/runpod_key

# Connect:
ssh -o StrictHostKeyChecking=no -i /tmp/runpod_key root@103.196.86.40 -p 38349

# Run commands remotely:
sshpass -p "1234567" ssh -o StrictHostKeyChecking=no -i /tmp/runpod_key root@103.196.86.40 -p 38349 "your command"
```

Note: sshpass password auth doesn't work (RunPod blocks it). Use the key only. sshpass is installed locally via brew.

---

## Currently Running on RunPod

**tmux session:** `esm2depmap`  
**Script:** `/evo2/dami-AI-Scientist/templates/esm2_depmap/experiment.py --out_dir run_0`  
**Log:** `/tmp/esm2_depmap_run.log`  

**What it does:**
1. ✅ Downloaded DepMap 26Q1 CRISPRGeneEffect.csv (421MB) to `run_0/data/`
2. ⏳ Fetching protein sequences for 2000 random human genes from Ensembl (~60 min total)
3. ⬜ Extract ESM-2 650M embeddings (~30 min on L40S)
4. ⬜ Compute pairwise cosine distances (ESM-2) and Pearson correlations (essentiality)
5. ⬜ Mantel test: Spearman correlation between the two distance matrices (1000 permutations)
6. ⬜ Write `run_0/final_info.json` with Mantel r and p-value

**Check progress:**
```bash
sshpass -p "1234567" ssh -o StrictHostKeyChecking=no -i /tmp/runpod_key root@103.196.86.40 -p 38349 "ls -lh /evo2/dami-AI-Scientist/templates/esm2_depmap/run_0/data/"
```
When `dataset_2000.json` appears, Ensembl fetch is done. When `embeddings_esm2_t33_650M_UR50D_*.npy` appears, embeddings are done. When `final_info.json` appears, it's complete.

**The key result to look for:**
```bash
cat /evo2/dami-AI-Scientist/templates/esm2_depmap/run_0/final_info.json
```
`mantel_r` > 0 and `mantel_p_value` < 0.05 = ESM-2 sequence space correlates with functional dependency space.

---

## Experimental Results So Far

### evo2_function (DNA repair vs TSG, 178 genes, Evo2-7B layer 28, DNA sequences)
- **Best AUROC: 0.78** | Silhouette: 0.044 (p=0.132, NS) | ARI: ~0 (p=0.516, NS)
- Supervised signal real; unsupervised geometry flat
- BER sub-pathway coherence p=0.015 (only significant sub-pathway of 4)
- Bridge genes (BRCA1/2, MLH1, MSH2, MSH6, PMS2, PALB2, ATM): closer to TSG centroid, p=0.084

### evo2_validation (glycolysis vs T-cell, 189 genes, Evo2-7B layer 28)
- **Best AUROC: 0.66** | Centroid distance: 0.0014
- Weaker than DNA repair/TSG — genome integrity better encoded than metabolic/immune

### esm2_function (DNA repair vs TSG, 178 genes, ESM-2 650M, protein sequences)
- **Best AUROC: 0.894** | Silhouette: -0.007 | ARI: -0.034 | Distance ratio: 0.938
- Significantly stronger supervised signal than Evo2 (0.89 vs 0.78)
- Unsupervised geometry still flat

### AI Scientist generated paper (complete PDF)
- `results/evo2_function/20260522_001652_unsupervised_functional_clustering/unsupervised_functional_clustering.pdf`
- Title: "Evo2 Embeddings Require Supervision to Extract Functional Information"
- High quality — correct framing, real stats, honest negative result

---

## Key Findings / Story

1. **Both Evo2 and ESM-2 encode functional information — but only extractable via supervised probing, not raw geometry.** ARI~0, silhouette~0 in all experiments.
2. **ESM-2 protein embeddings >> Evo2 DNA embeddings** for supervised classification (0.89 vs 0.78 AUROC on same gene set). Protein sequence space more functionally organised than DNA sequence space.
3. **Evo2's functional signal is contrast-dependent** — genome integrity (0.78) >> metabolic/immune (0.66).
4. **BER is anomalous** — base excision repair genes uniquely form a coherent sub-cluster (p=0.015). Other repair sub-pathways don't.
5. **Core claim (EvoMine framing):** functional information is in superposition in both models — needs learned extraction, not raw cosine similarity.

---

## Bug Fixes Already Applied to Repo

- `launch_scientist.py`: `novel_ideas` KeyError → `.get("novel", True)`
- `launch_scientist.py`: `baseline_results["means"]` KeyError → handles flat structure
- `launch_scientist.py`: `--no-writeup` flag added (writes markdown findings instead of LaTeX+review)
- `ai_scientist/perform_experiments.py`: same `means` KeyError fix
- `ai_scientist/llm.py`: added `claude-sonnet-4-5`, `claude-haiku-4-5`
- RunPod `vortex/model/utils.py`: `weights_only=False` for Evo2 checkpoint loading
- `templates/evo2_function/experiment.py`: cache fallback (copies baseline embeddings to avoid re-fetching)

## Common RunPod Issues
- `--skip-novelty-check` required (no S2 API key or rate limits)
- `--no-writeup` recommended (LaTeX is slow and expensive)
- Must run `python3 experiment.py --out_dir run_0` first before `launch_scientist.py`
- Aider tends to add new gene classes — `prompt.json` has explicit constraints to prevent this
- TexLive install needed for LaTeX: `apt-get install -y texlive-full`

---

## Project Files of Interest

- `/Users/dgupta/code/portfolio/EvoMine/evomine_findings.md` — EvoMine prior results
- `/Users/dgupta/code/portfolio/EvoMine/ai_scientist_run_20260521_findings.md` — today's Evo2 results
- `/Users/dgupta/code/portfolio/EvoMine/probeseq_findings.md` — ESM-2 vs Evo2 comparison
