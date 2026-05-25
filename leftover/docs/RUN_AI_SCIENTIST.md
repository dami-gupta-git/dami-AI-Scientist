# Running AI Scientist on RunPod

## Overview

The workflow is:
1. Generate hypotheses locally (`--ideas-only`)
2. SSH into RunPod, pull the branch, run the baseline experiment to build `run_0/`
3. Run `launch_scientist.py` — it uses Aider + Claude to implement and run each idea
4. Pull results back locally

---

## 1. SSH Key

The working key is `~/.ssh/id_runpod_2` (registered as `runpod_2` in RunPod account settings).

```
ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIMuXdSvZt602kwA42h3d78lyGwRgK35z4TA26mG8qlhw runpod_2
```

RunPod automatically injects this key into `authorized_keys` on pod start — no manual setup needed for new pods as long as `runpod_2` is registered in account settings (Settings → SSH Public Keys).

Connect (substitute IP and port from the RunPod console):
```bash
ssh -o StrictHostKeyChecking=no -i ~/.ssh/id_runpod_2 root@<POD_IP> -p <PORT>
```

SCP (for pulling results):
```bash
scp -i ~/.ssh/id_runpod_2 -P <PORT> root@<POD_IP>:<remote_path> <local_path>
```

---

## 2. One-Time Pod Setup

Run once on a fresh pod:

```bash
# Install system deps
apt-get update -y && apt-get install -y tmux git

# Install Python deps
pip3 install fair-esm scikit-learn scipy

# For Evo2 experiments only:
pip3 install flash-attn --no-build-isolation
sed -i 's/weights_only=True/weights_only=False/' \
  /usr/local/lib/python3.11/dist-packages/vortex/model/utils.py
```

---

## 3. Clone and Set Up the Repo

```bash
git clone https://github.com/dami-gupta-git/dami-AI-Scientist.git /workspace/dami-AI-Scientist
cd /workspace/dami-AI-Scientist
git checkout <branch>
```

---

## 4. Generate Hypotheses (locally, on your Mac)

```bash
ANTHROPIC_API_KEY=<key> python3 launch_scientist.py \
  --experiment <experiment> \
  --ideas-only \
  --skip-novelty-check \
  --num-ideas 3
```

Ideas are saved to `templates/<experiment>/ideas.json`. Review and edit before pushing:
```bash
git add templates/<experiment>/ideas.json && git commit -m "Update ideas" && git push
```

---

## 5. Run the Baseline Experiment on RunPod

The baseline must run first — it fetches data, extracts embeddings, and writes `run_0/final_info.json` which AI Scientist uses as the starting point.

```bash
tmux new-session -d -s baseline \
  'cd /workspace/dami-AI-Scientist && \
   git pull && \
   python3 templates/<experiment>/experiment.py --out_dir templates/<experiment>/run_0 \
   2>&1 | tee /tmp/baseline.log; echo DONE >> /tmp/baseline.log'

# Monitor
ssh -o StrictHostKeyChecking=no -i ~/.ssh/id_runpod_2 root@<POD_IP> -p <PORT> "tail -f /tmp/baseline.log"
```

---

## 6. Run AI Scientist

```bash
tmux new-session -d -s aiscientist \
  'cd /workspace/dami-AI-Scientist && \
   ANTHROPIC_API_KEY=<key> OPENAI_API_KEY=<key> \
   python3 launch_scientist.py \
     --experiment <experiment> \
     --model claude-sonnet-4-5 \
     --skip-idea-generation \
     --skip-novelty-check \
     --no-writeup \
   2>&1 | tee /tmp/aiscientist.log; echo DONE >> /tmp/aiscientist.log'

# Monitor
ssh -o StrictHostKeyChecking=no -i ~/.ssh/id_runpod_2 root@<POD_IP> -p <PORT> "tail -f /tmp/aiscientist.log"
```

Key flags:
| Flag | Purpose |
|---|---|
| `--skip-idea-generation` | Use existing `ideas.json` instead of generating new ones |
| `--skip-novelty-check` | Skip Semantic Scholar API (no key needed) |
| `--no-writeup` | Skip LaTeX paper generation (faster) |
| `--num-ideas N` | Only run the first N ideas |
| `--model claude-sonnet-4-5` | Required — older claude-3-5-sonnet models return 404 |

To generate papers, omit `--no-writeup` and install LaTeX first:
```bash
apt-get install -y texlive-full   # ~2 GB, takes ~10 min
```

---

## 7. Pull Results Locally

```bash
scp -r -i ~/.ssh/id_runpod_2 -P <PORT> \
  root@<POD_IP>:/workspace/dami-AI-Scientist/results/<experiment>/ \
  results/<experiment>/
```

Each result folder contains:
- `notes.txt` — idea description + per-run results
- `log.txt` — full Aider + experiment log
- `final_info.json` — metrics
- `*.pdf` — paper (if `--no-writeup` was not set)
- `run_*/` — per-seed experiment outputs and plots

---

## 8. Experiments

| Branch | Experiment | Status |
|---|---|---|
| `esm2-mechanism` | ESM-2 delta-embeddings — GOF/DN/LOF mechanism geometry | Baseline running |
| `evo2-xgboost` | Evo2 XGBoost probe comparison | In progress |
| `evo2-supervised` | Evo2 supervised probe comparison (LR vs SVM vs MLP) | Previously run |
| `esm2-depmap` | ESM-2 vs DepMap Mantel test | Previously run |

---

## esm2_mechanism specifics

**Ideas setup:**
```bash
# Option A: generate ideas locally
ANTHROPIC_API_KEY=<key> python3 launch_scientist.py \
  --experiment esm2_mechanism --ideas-only --skip-novelty-check --num-ideas 3

# Option B: use seed ideas directly
cp templates/esm2_mechanism/seed_ideas.json templates/esm2_mechanism/ideas.json
```

**Expected baseline duration on A100 (~2-2.5 hours):**
- OSF dataset download + parse: ~5 min
- UniProt sequence fetch (~1200 genes): ~20 min
- Pfam family fetch (~1200 genes): ~20 min
- AlphaMissense scores (~8000 variants, rate-limited at 0.2s/variant): ~30-45 min
- ESM-2 650M embedding extraction (~8000 variant pairs): ~30-45 min
- Probes + baselines + orthogonality: ~15 min

---

## Current Pod (May 24 2026)

```bash
ssh -o StrictHostKeyChecking=no -i ~/.ssh/id_runpod_2 root@154.54.102.28 -p 13732
```

Repo at `/workspace/dami-AI-Scientist`, branch `esm2-mechanism`. MLP probe running in tmux session `mlp`.
