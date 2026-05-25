# Running AI Scientist on RunPod

## Overview

The workflow is:
1. Generate hypotheses locally (`--ideas-only`)
2. SSH into RunPod, pull the branch, run the baseline experiment to build `run_0/`
3. Run `launch_scientist.py` — it uses Aider + Claude to implement and run each idea
4. Pull results back locally

---

## 1. SSH Key Setup (one-time per pod)

The working key is `~/.ssh/id_runpod_2` (registered as `runpod_2` in RunPod account settings — injected automatically on pod start).

If the key isn't injected on a new pod, run this in the **pod web terminal**:
```bash
echo "ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIMuXdSvZt602kwA42h3d78lyGwRgK35z4TA26mG8qlhw runpod_2" > ~/.ssh/authorized_keys
chmod 600 ~/.ssh/authorized_keys
```

Then connect:
```bash
ssh -o StrictHostKeyChecking=no -i ~/.ssh/id_runpod_2 root@<POD_IP> -p <PORT>
```

---

## 2. One-Time Pod Setup

Run once on a fresh pod:

```bash
# Install system deps
apt-get update -y && apt-get install -y tmux

# Install Python deps
pip3 install -r requirements.txt

# Install flash-attn (required for Evo2, compiles from source — takes ~15 min)
pip3 install flash-attn --no-build-isolation

# Patch Evo2 checkpoint loader (required for PyTorch >= 2.4)
sed -i 's/weights_only=True/weights_only=False/' \
  /usr/local/lib/python3.11/dist-packages/vortex/model/utils.py
```

---

## 3. Clone and Set Up the Repo

```bash
cd /workspace
git clone https://github.com/dami-gupta-git/dami-AI-Scientist.git
cd dami-AI-Scientist
git checkout <branch>   # e.g. evo2-supervised
```

---

## 4. Generate Hypotheses (locally, on your Mac)

```bash
ANTHROPIC_API_KEY=<key> /opt/homebrew/bin/python3.11 launch_scientist.py \
  --experiment evo2_function \
  --ideas-only \
  --skip-novelty-check \
  --num-ideas 3
```

Ideas are saved to `templates/<experiment>/ideas.json`. Review and edit them before pushing.

Push to remote:
```bash
git add templates/<experiment>/ideas.json
git commit -m "Update ideas"
git push
```

---

## 5. Run the Baseline Experiment on RunPod

The baseline must be run first — it fetches gene sequences, extracts embeddings, and writes `run_0/final_info.json` which AI Scientist uses as the starting point.

```bash
# On RunPod, in a tmux session so it survives disconnection
tmux new-session -d -s baseline \
  'cd /workspace/dami-AI-Scientist && \
   git pull && \
   python3 templates/evo2_function/experiment.py --out_dir templates/evo2_function/run_0 \
   2>&1 | tee /tmp/baseline.log; echo DONE >> /tmp/baseline.log'

# Monitor
tail -f /tmp/baseline.log
```

Expected duration: ~1-2 hours (gene fetch ~45 min + Evo2 model download ~30 min + embedding extraction ~30 min).

When complete, `run_0/` will contain:
- `data/dataset.json` — gene sequences and labels (178 genes)
- `data/embeddings_evo2_7b.npy` — Evo2 embeddings (178 × 4096)
- `final_info.json` — baseline metrics (LR/SVM/MLP AUROCs)

Pull the cached data locally to avoid re-fetching next time:
```bash
# On your Mac
scp -i ~/.ssh/id_runpod_2 -P <PORT> \
  root@<POD_IP>:/workspace/dami-AI-Scientist/templates/evo2_function/run_0/data/dataset.json \
  templates/evo2_function/run_0/data/dataset.json
```

---

## 6. Run AI Scientist

```bash
# On RunPod
tmux new-session -d -s aiscientist \
  'cd /workspace/dami-AI-Scientist && \
   ANTHROPIC_API_KEY=<key> OPENAI_API_KEY=<key> \
   python3 launch_scientist.py \
     --experiment evo2_function \
     --model claude-sonnet-4-5 \
     --skip-idea-generation \
     --skip-novelty-check \
     --no-writeup \
   2>&1 | tee /tmp/aiscientist.log; echo DONE >> /tmp/aiscientist.log'

# Monitor
tail -f /tmp/aiscientist.log
```

Key flags:
| Flag | Purpose |
|---|---|
| `--skip-idea-generation` | Use existing `ideas.json` instead of generating new ones |
| `--skip-novelty-check` | Skip Semantic Scholar API (no key needed) |
| `--no-writeup` | Skip LaTeX paper generation (faster) |
| `--num-ideas N` | Only run the first N novel ideas |
| `--model claude-sonnet-4-5` | Required — older claude-3-5-sonnet models return 404 |

To generate papers, omit `--no-writeup` and install LaTeX first:
```bash
apt-get install -y texlive-full   # ~2 GB, takes ~10 min
```

---

## 7. Pull Results Locally

Results are written to `results/<experiment>/<timestamp>_<idea_name>/`.

```bash
# On your Mac — pull all results for an experiment
scp -r -i ~/.ssh/id_runpod_2 -P <PORT> \
  root@<POD_IP>:/workspace/dami-AI-Scientist/results/evo2_function/ \
  results/evo2_function/
```

Each result folder contains:
- `notes.txt` — idea description + per-run results
- `log.txt` — full Aider + experiment log
- `final_info.json` — metrics
- `*.pdf` — paper (if `--no-writeup` was not set)
- `run_*/` — per-seed experiment outputs and plots

---

## 8. Current Experiments

| Branch | Experiment | Status |
|---|---|---|
| `evo2-xgboost` | Evo2 XGBoost probe comparison | In progress |
| `esm2-mechanism` | ESM-2 delta-embeddings — GOF/DN/LOF mechanism geometry | Ready to run |
| `evo2-supervised` | Evo2 supervised probe comparison (LR vs SVM vs MLP) | Previously run |
| `esm2-depmap` | ESM-2 vs DepMap Mantel test | Previously run |

## Running esm2_mechanism

**Before running, set up ideas:**
```bash
# Option A: generate ideas locally first
ANTHROPIC_API_KEY=<key> python3 launch_scientist.py \
  --experiment esm2_mechanism \
  --ideas-only \
  --skip-novelty-check \
  --num-ideas 3

# Option B: use seed_ideas.json directly
cp templates/esm2_mechanism/seed_ideas.json templates/esm2_mechanism/ideas.json
```

**Baseline run on RunPod (~2-2.5 hours on A100):**
```bash
tmux new-session -d -s baseline \
  'cd /workspace/dami-AI-Scientist && \
   git pull && \
   python3 templates/esm2_mechanism/experiment.py --out_dir templates/esm2_mechanism/run_0 \
   2>&1 | tee /tmp/baseline.log; echo DONE >> /tmp/baseline.log'

tail -f /tmp/baseline.log
```

Time breakdown:
- OSF dataset download + parse: ~5 min
- UniProt sequence fetch (~1200 genes): ~20 min
- Pfam family fetch (~1200 genes): ~20 min
- AlphaMissense scores (~8000 variants, rate-limited): ~30-45 min
- ESM-2 650M embedding extraction (~8000 variant pairs): ~30-45 min
- Probes + baselines + orthogonality: ~15 min

**AI Scientist run:**
```bash
tmux new-session -d -s aiscientist \
  'cd /workspace/dami-AI-Scientist && \
   ANTHROPIC_API_KEY=<key> OPENAI_API_KEY=<key> \
   python3 launch_scientist.py \
     --experiment esm2_mechanism \
     --model claude-sonnet-4-5 \
     --skip-idea-generation \
     --skip-novelty-check \
     --no-writeup \
   2>&1 | tee /tmp/aiscientist.log; echo DONE >> /tmp/aiscientist.log'
```

## Current RunPod Connection (May 25 2026)

```bash
ssh -o StrictHostKeyChecking=no -i ~/.ssh/id_runpod_2 root@216.81.245.143 -p 11019
```

A100 SXM4 80GB, 128 cores.
