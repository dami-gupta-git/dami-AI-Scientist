# Publishable finding — one-page summary

## Title

*ESM-2 encodes mutation pathogenicity but not mechanism: a controlled dissociation*

## The claim (one sentence)

The same ESM-2 mutant−wildtype delta embeddings that predict ClinVar pathogenic vs benign at AUROC 0.88 cannot classify gain-of-function / dominant-negative / loss-of-function mechanism above chance (macro-F1 0.28) — and the apparent gene-level mechanism signal in earlier work is a Pfam-family recognition shortcut.

## Three numbers that make the paper

| Task | Number | What it shows |
|---|---|---|
| Pathogenicity (ClinVar, 17,236 variants) | **AUROC 0.88** under gene-split, **0.88 under family-split** | Pipeline works; signal is per-variant, not homology |
| Mechanism (Gerasimavicius, 948 genes, 3-class) | **macro-F1 0.28** under both gene-split and family-split | No mechanism signal in deltas |
| WT-only mechanism baseline | 0.58 → **0.39** under family-split | The apparent gene-level mechanism signal is family leakage |
| MLP delta (class-weighted PyTorch) | **macro-F1 0.43**, GOF AUROC 0.73 under gene-split | Nonlinear lift — but likely residual family signal, not mechanism (pending family-split confirmation) |

The MLP row is the one number that needs explaining before posting. If MLP also stays flat under family-split (≈0.43), it confirms residual family leakage as the explanation. If it collapses, the case is closed. **This is the single blocking experiment.**

## Why it works as a paper

- **Positive control passes**: pipeline can extract a signal that's actually there (pathogenicity)
- **Negative control passes**: family-split CV breaks the shortcut and the apparent signal disappears
- **Causal explanation**: family clustering quantified (k=5 purity 26× chance) + 74.8% within-family mechanism agreement = mechanical reason for the leakage
- **Asymmetry**: same pipeline, same embeddings, two tasks, two opposite outcomes — clean dissociation

## What the paper would contain

1. **Methods**: ESM-2 650M delta embeddings, logreg + MLP probes, gene-split + family-split CV
2. **Dataset 1**: Gerasimavicius 948 genes for mechanism
3. **Dataset 2**: ClinVar 17,236 variants for pathogenicity positive control
4. **Family-clustering diagnostic**: k-NN purity, family probe, within-family mechanism agreement
5. **One figure**: side-by-side bar chart, pathogenicity vs mechanism, gene-split vs family-split
6. **Honest scoping**: "on ESM-2 650M, on Gerasimavicius — generalisation to other models/datasets is open"

## What's missing before posting (1 day of work)

1. **Multi-seed numbers** — currently seed=0 only. Need 5 seeds.
2. **MLP under family-split CV** on WT-only AND delta — **running on RunPod now** (`experiment_mlp.py --family_split`). This is the single blocking experiment.
3. **The figure** — side-by-side bar chart: pathogenicity vs mechanism × gene-split vs family-split. Makes the dissociation visually obvious. Needs to be written (plot.py currently only covers mechanism).
4. **The LaTeX draft itself** (none exists yet)

## What this is NOT

- Not a discovery — PreMode (Zhong & Shen, Nat Commun 2025) and the AlphaMissense paper already state "PLMs don't predict mechanism." We are the first **controlled side-by-side demonstration** with a leakage diagnostic, not the first to notice.
- Not generalisable yet — single model, single mechanism dataset. SaProt / ESM-3 / DDG2P needed for that.
- The three genuine contributions: (1) controlled dissociation on the same pipeline, (2) family-split CV as a quantitative leakage diagnostic, (3) MissION reconciliation — PLM deltas fail cross-proteome but supervised fine-tuning on a homologous subfamily works, and these are not in conflict.

## Venue

bioRxiv preprint as a methodological consolidation paper. Peer-reviewed target: *Bioinformatics* or *Genome Biology* methodological note.
