# Session Compact 2 — May 23, 2026

## Who You Are Working With
Dami Gupta — computational biology researcher. Prefers concise responses, no emojis, no trailing summaries. Does not want layer ablation experiments. Global rules: never fabricate scientific values, no fallback defaults, surface disagreements before acting.

---

## Project State

Using AI Scientist framework to auto-generate and run research hypotheses on Evo2 DNA embeddings. Core question: **do Evo2 embeddings encode functional gene class information, and what type of classifier best extracts it?**

### Repo
- **Local:** `/Users/dgupta/code/downloads/AIScientist/AI-Scientist`
- **Remote:** `https://github.com/dami-gupta-git/dami-AI-Scientist.git`
- **Active branch:** `evo2-xgboost` (also pushed: `evo2-supervised`)
- `results/` is gitignored — pull manually via scp

### RunPod Connections
| Pod | Purpose | SSH |
|---|---|---|
| Main (Evo2) | Experiments + AI Scientist | `ssh -i ~/.ssh/id_runpod root@216.81.245.125 -p 13793` |
| LaTeX pod | Paper writeup only | `ssh -i ~/.ssh/id_runpod root@194.68.245.123 -p 22025` |

Both use `~/.ssh/id_runpod` key (added to `authorized_keys` via web terminal).

---

## What's Running Right Now

1. **LaTeX writeup** (LaTeX pod) — `supervised_probe_comparison` paper, compiling PDF. Check: `tail -f /tmp/writeup3.log`
2. **XGBoost AI Scientist run** — DONE. 3 ideas completed: `xgboost_hyperparameter_tuning`, `model_complementarity_ensemble`, `embedding_geometry_analysis`

---

## Experimental Results This Session

### Experiment: `evo2-supervised` — supervised probe comparison (DNA repair vs TSG, 178 genes, Evo2-7B layer 28)

**Multi-seed (5 seeds) — reliable:**
| Probe | Test AUROC | CV AUROC |
|---|---|---|
| LR | 0.661 | 0.696 ± 0.091 |
| SVM (RBF) | 0.722 | 0.742 ± 0.113 |
| MLP | 0.740 | — |
| Nonlinearity gap (MLP-LR) | +0.078 | — |

**Single seed — noisy, for reference only:**
| Probe | Test AUROC | CV AUROC |
|---|---|---|
| LR | 0.642 | 0.689 ± 0.116 |
| SVM | 0.675 | 0.717 ± 0.146 |
| XGBoost (default) | 0.609 | 0.722 ± 0.100 |
| XGBoost (tuned grid search) | 0.626 | **0.746 ± 0.127** |
| MLP | 0.621 | — |

**Key conclusion:** XGBoost has best CV AUROC (0.746) on small dataset. MLP needs more data. Test set results (n=36) are too noisy for single-seed conclusions — always use multi-seed.

### AI Scientist ideas run this session (evo2-xgboost branch)
- `xgboost_hyperparameter_tuning` — grid search over 81 configs, best CV AUROC 0.746
- `model_complementarity_ensemble` — ensemble of all 4 models
- `embedding_geometry_analysis` — geometry of embedding space

---

## Key Open Issue: Repeat Contamination

From `papers/evomine_findings.md`: raw Evo2 embeddings are dominated by repetitive elements (Alu/LINE). DNA repair genes may be systematically more repeat-rich than TSGs, meaning classifiers could be learning repeat content rather than function.

**Fix:** N-mask repeats before Evo2 embedding extraction, re-run probes, compare AUROC. If AUROC drops → repeats were the signal. If stable → signal is real.

**Status:** Not yet implemented. Dami wants to finish current runs first, then decide.

---

## Branch Structure
| Branch | Experiment | Status |
|---|---|---|
| `evo2-xgboost` | XGBoost + LR + SVM + MLP probe comparison | Done; results pulled locally |
| `evo2-supervised` | LR + SVM + MLP probe comparison | Done; paper being written |
| `evo2-template` | Original Evo2 MLP baseline | Old |
| `esm2-depmap` | ESM-2 vs DepMap Mantel test | Previously run |

---

## Setup Notes (for new pods)

```bash
# Install deps
apt-get update -y && apt-get install -y tmux
pip3 install anthropic openai aider-chat backoff scikit-learn scipy matplotlib pymupdf4llm pypdf google-generativeai xgboost
pip3 install evo2
pip3 install flash-attn --no-build-isolation  # ~15 min compile
# Patch Evo2 checkpoint loader
sed -i 's/weights_only=True/weights_only=False/' \
  /usr/local/lib/python3.11/dist-packages/vortex/model/utils.py

# LaTeX pod only
apt-get install -y texlive texlive-latex-extra texlive-science chktex
dpkg --configure -a  # if install interrupted
```

## Run AI Scientist
```bash
ANTHROPIC_API_KEY=<key> OPENAI_API_KEY=<key> python3 launch_scientist.py \
  --experiment evo2_function \
  --model claude-sonnet-4-5 \
  --skip-novelty-check \
  --no-writeup \
  --num-ideas 2
```

## Run writeup (from repo root)
```bash
python3 -m ai_scientist.perform_writeup \
  --folder results/evo2_function/<folder_name> \
  --model claude-sonnet-4-5
```

---

## Results Location
- **Local:** `results/evo2_function/<timestamp>_<idea_name>/`
- **Pod:** `/workspace/dami-AI-Scientist/results/evo2_function/`
- Pull: `scp -r -i ~/.ssh/id_runpod -P 13793 root@216.81.245.125:/workspace/dami-AI-Scientist/results/evo2_function/<folder> results/evo2_function/`

## API Keys
- `ANTHROPIC_API_KEY` — claude-sonnet-4-5 (only model that works with Dami's key)
- `OPENAI_API_KEY` — for review step
