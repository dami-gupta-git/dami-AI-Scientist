# ESM Fine-Tune: LoRA Contrastive Fine-Tuning for Disease Mechanism

## Scientific question

Can LoRA fine-tuning of ESM-2 — with a cross-family supervised contrastive loss — break through the mechanism classification ceiling established by frozen-embedding probing?

**Background:** The esm2_mechanism project established that frozen ESM-2 encodes pathogenicity strongly (AUROC 0.88) but disease mechanism weakly (family-split macro-F1 floor ~0.36–0.40). Result_9 showed that a contrastive projection head on frozen embeddings reaches F1=0.397 under family-split — the current best. That projection head could only rearrange existing directions in the frozen 1280-d space.

**This experiment:** LoRA fine-tuning backprops through ESM-2's attention layers directly, allowing the model to reorganise representations around mechanism rather than just reading off existing directions. If mechanism signal is entangled with family identity in ESM-2's attention layers (result_4: family k-purity z=+18), gradient flow via a family-invariant contrastive loss may disentangle it in a way a frozen projection head cannot.

## Method

### Loss
Supervised contrastive with cross-family triplets (TripletMarginLoss, margin=1.0):
- **Positives:** same mechanism, different Pfam family
- **Negatives:** different mechanism
- **Within-family pairs excluded from positives** — forces family-invariant mechanism structure

### Architecture
- ESM-2 backbone (650M or 3B) with LoRA adapters on q_proj and v_proj
- LoRA rank=8, alpha=16, dropout=0.1
- Delta embedding = mean_pool(mut) − mean_pool(wt) from fine-tuned model
- Evaluation: k-NN (k=10, cosine) + linear probe on delta embeddings

### Cross-validation
- 5-fold gene-split AND 5-fold family-split (same scheme as all prior results)
- LoRA adapters re-trained fresh for each fold (no leakage)

### Dataset
- Merged dataset: 19,100 variants, 1,985 genes, 1,146 Pfam families
- 3-class labels: GOF / HI / LOF (AR variants mapped to LOF/HI per Gerasimavicius labels)

## Runs

### Pod 1 — ESM-2 650M (A100 SXM 80GB)
- SSH: `root@195.26.233.76 -p 21033 -i ~/.ssh/id_runpod_2`
- Command:
  ```bash
  python3 /workspace/lora_contrastive_finetune.py \
    --data_dir /workspace/data \
    --out_dir /workspace/results \
    --model_name esm2_t33_650M_UR50D \
    --lora_rank 8 --lora_alpha 16 \
    --lr 2e-4 --max_epochs 20 --patience 5 \
    --batch_size 16 --grad_accum 4 \
    --embed_batch 16 --max_len 1022 \
    --merged --seed 0
  ```
- Output: `/workspace/results/lora_contrastive_t33_650M_merged_seed0.json`

### Pod 2 — ESM-2 3B (H100 80GB)
- SSH: `root@31.24.80.36 -p 15325 -i ~/.ssh/id_runpod_2`
- Command:
  ```bash
  python3 /workspace/lora_contrastive_finetune.py \
    --data_dir /workspace/data \
    --out_dir /workspace/results \
    --model_name esm2_t36_3B_UR50D \
    --lora_rank 8 --lora_alpha 16 \
    --lr 2e-4 --max_epochs 20 --patience 5 \
    --batch_size 8 --grad_accum 8 \
    --embed_batch 8 --max_len 1022 \
    --merged --seed 0
  ```
- Output: `/workspace/results/lora_contrastive_t36_3B_merged_seed0.json`

## Baselines to beat (from esm2_mechanism)

| Method | Family-split macro-F1 |
|---|---|
| MLP on frozen delta (result_7) | 0.364 |
| Contrastive projection head, frozen (result_9) | 0.397 ← current best |
| **LoRA 650M (this experiment)** | TBD |
| **LoRA 3B (this experiment)** | TBD |

**Threshold:** family-split F1 > 0.427 (result_9 + 0.03) would be a clear improvement.

## Key diagnostic

The critical check is whether the gene-split → family-split drop is smaller for LoRA than for the frozen MLP. If LoRA reduces the leakage fraction, it confirms the fine-tuning is finding family-invariant mechanism signal rather than exploiting family shortcuts via a different route.

| Method | Gene-split F1 | Family-split F1 | Δ (leakage proxy) |
|---|---|---|---|
| MLP frozen (result_7) | 0.415 | 0.364 | 0.051 |
| Contrastive proj frozen (result_9) | 0.470 | 0.397 | 0.073 |
| LoRA 650M | TBD | TBD | TBD |
| LoRA 3B | TBD | TBD | TBD |

## Files

- `scripts/lora_contrastive_finetune.py` — training script
- `docs/EXPERIMENT.md` — this file
- `results/` — JSON outputs pulled from pods after runs complete
