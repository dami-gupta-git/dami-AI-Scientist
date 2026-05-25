"""
LoRA contrastive fine-tuning of ESM-2 for disease mechanism classification.

Unlike result_9 (projection head on FROZEN embeddings), this script backprops
through ESM-2 itself via LoRA adapters, allowing the model to reorganise its
representations around mechanism rather than just reading off existing directions.

Loss: supervised contrastive (cross-family triplets)
  - Positives: same mechanism, DIFFERENT Pfam family
  - Negatives: different mechanism
  - Within-family pairs excluded from positives

Evaluation: linear probe + k-NN on the fine-tuned [CLS]/mean-pool representations,
under both gene-split and family-split 5-fold CV.

Usage:
    python lora_contrastive_finetune.py \
        --data_dir ../data \
        --out_dir ../results/20260524_baseline_run/run_0 \
        --model_name esm2_t33_650M_UR50D \
        --lora_rank 8 \
        --seed 0

    # For ESM-2 3B:
    python lora_contrastive_finetune.py \
        --model_name esm2_t36_3B_UR50D \
        --lora_rank 8 \
        --batch_size 16 \
        --grad_accum 8

Outputs:
    lora_contrastive_results_{model_tag}_seed{seed}.json
"""

import argparse
import json
import os
import warnings
from collections import Counter

import numpy as np
from sklearn.metrics import roc_auc_score, f1_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import LabelEncoder

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Data loading — variants + sequences (no cached embeddings; we re-encode)
# ---------------------------------------------------------------------------

def load_data(data_dir, merged=False):
    mv_path = os.path.join(data_dir, "embeddings", "merged_valid_variants.json")
    if not os.path.exists(mv_path):
        mv_path = os.path.join(data_dir, "merged_valid_variants.json")

    with open(mv_path) as f:
        mv = json.load(f)

    if merged:
        variants = mv
        print(f"Loaded {len(variants)} merged variants")
    else:
        variants = [v for v in mv if v.get("source") == "gerasimavicius"]
        print(f"Loaded {len(variants)} Gerasimavicius variants")

    # Load sequences
    seq_path = os.path.join(data_dir, "sequences.json")
    with open(seq_path) as f:
        sequences = json.load(f)

    labels = np.array([v["label_3class"] for v in variants])
    genes  = np.array([v["gene"] for v in variants])

    print(f"Class distribution: {dict(Counter(labels))}")
    print(f"Unique genes: {len(set(genes))}")
    return variants, labels, genes, sequences


def load_pfam(data_dir, genes):
    pfam_path = os.path.join(data_dir, "pfam_families.json")
    with open(pfam_path) as f:
        pfam_map = json.load(f)
    gene_pfam = np.array([pfam_map.get(g) for g in genes])
    print(f"Pfam coverage: {sum(p is not None for p in gene_pfam)}/{len(genes)}")
    return gene_pfam, pfam_map

# ---------------------------------------------------------------------------
# CV splits (same logic as contrastive_mechanism.py)
# ---------------------------------------------------------------------------

def family_split_cv(genes, gene_pfam, pfam_map, n_folds=5, seed=42):
    gene_to_pfam = {g: pfam_map.get(g) for g in np.unique(genes)}
    annotated_mask = np.array([gene_to_pfam.get(g) is not None for g in genes])
    unique_fams = sorted(set(v for v in gene_to_pfam.values() if v is not None))
    rng = np.random.RandomState(seed)
    fam_arr = np.array(unique_fams)
    rng.shuffle(fam_arr)
    splits = []
    for fold_fams in np.array_split(fam_arr, n_folds):
        fold_set = set(fold_fams)
        te = annotated_mask & np.array([gene_to_pfam.get(g) in fold_set for g in genes])
        tr = annotated_mask & ~te
        if tr.sum() >= 20 and te.sum() >= 10:
            splits.append((np.where(tr)[0], np.where(te)[0]))
    print(f"Family-split: {len(splits)} folds, {len(unique_fams)} families")
    return splits


def gene_split_cv(genes, n_folds=5, seed=42):
    unique_genes = np.array(sorted(set(genes)))
    rng = np.random.RandomState(seed)
    rng.shuffle(unique_genes)
    splits = []
    for fold_genes in np.array_split(unique_genes, n_folds):
        fold_set = set(fold_genes)
        te = np.array([g in fold_set for g in genes])
        tr = ~te
        if tr.sum() >= 20 and te.sum() >= 10:
            splits.append((np.where(tr)[0], np.where(te)[0]))
    return splits

# ---------------------------------------------------------------------------
# Cross-family triplet construction (same logic as result_9)
# ---------------------------------------------------------------------------

def build_cross_family_triplets(labels, gene_pfam, le, max_pairs_per_anchor=8, seed=42):
    rng = np.random.RandomState(seed)
    y = le.transform(labels)
    n = len(labels)
    n_classes = len(le.classes_)

    unique_fams = list({f for f in gene_pfam if f is not None})
    fam_to_int = {f: i for i, f in enumerate(unique_fams)}
    fam_int = np.array([fam_to_int.get(f, -1) for f in gene_pfam], dtype=np.int32)

    by_mech = {c: np.where(y == c)[0] for c in range(n_classes)}
    by_class_arr = {c: np.array(by_mech[c]) for c in range(n_classes)}
    neg_by_class = {c: np.concatenate([by_mech[o] for o in range(n_classes) if o != c])
                    for c in range(n_classes)}

    unique_combos = set((int(y[i]), int(fam_int[i])) for i in range(n))
    combo_pos_pool = {}
    for (c, fam) in unique_combos:
        class_idxs = by_class_arr[c]
        class_fams = fam_int[class_idxs]
        cross_fam_mask = (class_fams != fam) & (class_fams != -1)
        combo_pos_pool[(c, fam)] = class_idxs[cross_fam_mask]

    anchor_list, pos_list, neg_list = [], [], []
    for anchor_i in range(n):
        c = int(y[anchor_i])
        fam = int(fam_int[anchor_i])
        pos_pool = combo_pos_pool.get((c, fam), np.array([], dtype=np.int64))
        neg_pool = neg_by_class[c]
        if len(pos_pool) == 0 or len(neg_pool) == 0:
            continue
        n_pairs = min(max_pairs_per_anchor, len(pos_pool), len(neg_pool))
        anchor_list.append(np.full(n_pairs, anchor_i, dtype=np.int64))
        pos_list.append(pos_pool[rng.randint(0, len(pos_pool), n_pairs)])
        neg_list.append(neg_pool[rng.randint(0, len(neg_pool), n_pairs)])

    anchors   = np.concatenate(anchor_list)
    positives = np.concatenate(pos_list)
    negatives = np.concatenate(neg_list)
    print(f"  Built {len(anchors)} cross-family triplets "
          f"({len(set(anchors.tolist()))} unique anchors)")
    return anchors, positives, negatives

# ---------------------------------------------------------------------------
# ESM-2 tokenisation + embedding extraction
# ---------------------------------------------------------------------------

def load_esm2(model_name, device):
    import esm
    model_map = {
        "esm2_t33_650M_UR50D": esm.pretrained.esm2_t33_650M_UR50D,
        "esm2_t36_3B_UR50D":   esm.pretrained.esm2_t36_3B_UR50D,
        "esm2_t30_150M_UR50D": esm.pretrained.esm2_t30_150M_UR50D,
    }
    if model_name not in model_map:
        raise ValueError(f"Unknown model: {model_name}. Choose from {list(model_map)}")
    print(f"Loading {model_name}...")
    model, alphabet = model_map[model_name]()
    model = model.to(device)
    return model, alphabet


def embed_sequences(model, alphabet, sequences_list, device, batch_size=8, max_len=1022):
    """
    sequences_list: list of (label, seq_str) tuples (ESM batch format)
    Returns mean-pool representations, shape (N, embed_dim).
    """
    import torch
    batch_converter = alphabet.get_batch_converter()
    model.eval()

    all_reps = []
    for i in range(0, len(sequences_list), batch_size):
        batch = sequences_list[i : i + batch_size]
        # Truncate to max_len
        batch = [(lbl, seq[:max_len]) for lbl, seq in batch]
        _, _, tokens = batch_converter(batch)
        tokens = tokens.to(device)

        with torch.no_grad():
            out = model(tokens, repr_layers=[model.num_layers], return_contacts=False)
        reps = out["representations"][model.num_layers]  # (B, L, D)

        # Mean pool over sequence positions (exclude BOS/EOS)
        for j, (_, seq) in enumerate(batch):
            seq_len = min(len(seq), max_len)
            rep = reps[j, 1:seq_len+1, :].mean(0).cpu().float()
            all_reps.append(rep)

    return torch.stack(all_reps)  # (N, D)


def get_variant_sequences(variants, sequences):
    """Return (wt_seq, mut_seq) for each variant."""
    wt_seqs, mut_seqs = [], []
    missing = 0
    for v in variants:
        # sequences.json is keyed by UniProt ID
        seq  = sequences.get(v.get("uniprot_id")) or sequences.get(v.get("gene"))
        if seq is None:
            wt_seqs.append(None)
            mut_seqs.append(None)
            missing += 1
            continue
        pos  = v["aa_pos"] - 1  # 0-indexed
        wt   = v["aa_wt"]
        mut  = v["aa_mut"]
        if pos < 0 or pos >= len(seq):
            wt_seqs.append(None)
            mut_seqs.append(None)
            missing += 1
            continue
        mut_seq = seq[:pos] + mut + seq[pos+1:]
        wt_seqs.append(seq)
        mut_seqs.append(mut_seq)

    if missing:
        print(f"  WARNING: {missing}/{len(variants)} variants missing sequence")
    return wt_seqs, mut_seqs

# ---------------------------------------------------------------------------
# Pure-PyTorch LoRA (no peft) — compatible with fair-esm native interface
# ---------------------------------------------------------------------------

import torch
import torch.nn as nn


class LoRALinear(nn.Module):
    """Wraps an existing nn.Linear with a low-rank LoRA adapter."""
    def __init__(self, linear: nn.Linear, rank: int, alpha: float, dropout: float = 0.0):
        super().__init__()
        self.linear   = linear
        self.rank      = rank
        self.scale     = alpha / rank
        in_f, out_f    = linear.in_features, linear.out_features
        self.lora_A    = nn.Parameter(torch.zeros(rank, in_f))
        self.lora_B    = nn.Parameter(torch.zeros(out_f, rank))
        self.dropout   = nn.Dropout(dropout)
        nn.init.kaiming_uniform_(self.lora_A, a=np.sqrt(5))

    def forward(self, x):
        base = self.linear(x)
        lora = self.dropout(x) @ self.lora_A.T @ self.lora_B.T
        return base + lora * self.scale


def inject_lora(model, rank, alpha, dropout):
    """Replace q_proj and v_proj in every ESM-2 attention layer with LoRALinear."""
    n_injected = 0
    for layer in model.layers:
        attn = layer.self_attn
        for proj_name in ("q_proj", "v_proj"):
            orig = getattr(attn, proj_name)
            lora_layer = LoRALinear(orig, rank, alpha, dropout)
            setattr(attn, proj_name, lora_layer)
            n_injected += 1
    # Freeze all params, then unfreeze only LoRA
    for p in model.parameters():
        p.requires_grad_(False)
    for layer in model.layers:
        attn = layer.self_attn
        for proj_name in ("q_proj", "v_proj"):
            lora_layer = getattr(attn, proj_name)
            lora_layer.lora_A.requires_grad_(True)
            lora_layer.lora_B.requires_grad_(True)
    trainable = sum(p.numel() for p in model.parameters() if p.requires_grad)
    total     = sum(p.numel() for p in model.parameters())
    print(f"  LoRA injected into {n_injected} projections. "
          f"Trainable: {trainable:,} / {total:,} ({100*trainable/total:.2f}%)")
    return model


def train_lora_fold(model, alphabet, wt_seqs, mut_seqs, labels, gene_pfam,
                    train_idx, le, device,
                    lora_rank=8, lora_alpha=16, lora_dropout=0.1,
                    lr=2e-4, max_epochs=20, patience=5,
                    batch_size=16, grad_accum=4, margin=1.0,
                    max_pairs_per_anchor=8, embed_batch=8, max_len=1022,
                    seed=42):
    """
    Fine-tune ESM-2 with LoRA on training fold using cross-family contrastive loss.
    Returns the fine-tuned model.
    """
    import copy
    from torch.utils.data import DataLoader, TensorDataset

    num_layers = model.num_layers

    fold_model = inject_lora(copy.deepcopy(model), lora_rank, lora_alpha, lora_dropout)
    fold_model = fold_model.to(device)

    optimizer = torch.optim.AdamW(
        [p for p in fold_model.parameters() if p.requires_grad],
        lr=lr, weight_decay=1e-2
    )
    triplet_loss_fn = nn.TripletMarginLoss(margin=margin, p=2)

    labels_tr   = labels[train_idx]
    gene_pfam_tr = gene_pfam[train_idx]

    # Build triplet indices (into train_idx)
    anchors, positives, negatives = build_cross_family_triplets(
        labels_tr, gene_pfam_tr, le,
        max_pairs_per_anchor=max_pairs_per_anchor, seed=seed
    )

    # Val split: 15% of triplets
    rng = np.random.RandomState(seed + 99)
    n_val = max(50, len(anchors) // 7)
    val_mask = np.zeros(len(anchors), dtype=bool)
    val_mask[rng.choice(len(anchors), n_val, replace=False)] = True
    tr_mask = ~val_mask

    anc_tr = anchors[tr_mask]; pos_tr = positives[tr_mask]; neg_tr = negatives[tr_mask]
    anc_val = anchors[val_mask]; pos_val = positives[val_mask]; neg_val = negatives[val_mask]

    ds = TensorDataset(torch.tensor(anc_tr), torch.tensor(pos_tr), torch.tensor(neg_tr))
    loader = DataLoader(ds, batch_size=batch_size * grad_accum, shuffle=True)

    # Precompute wt/mut sequences for training variants
    tr_wt  = [wt_seqs[i]  for i in train_idx]
    tr_mut = [mut_seqs[i] for i in train_idx]

    best_val_loss = float("inf")
    patience_count = 0
    best_state = None

    print(f"    LoRA training: {len(anc_tr)} train triplets, {len(anc_val)} val triplets")
    print(f"    Trainable params: {sum(p.numel() for p in fold_model.parameters() if p.requires_grad):,}")

    batch_converter = alphabet.get_batch_converter()

    def encode_batch(seqs_list):
        """
        Encode a list of sequences through fold_model, return mean-pool reps.
        seqs_list: list of seq strings (truncated to max_len).
        Returns tensor of shape (len(seqs_list), embed_dim).
        """
        labeled = [(str(i), s[:max_len]) for i, s in enumerate(seqs_list)]
        _, _, toks = batch_converter(labeled)
        out = fold_model(toks.to(device), repr_layers=[num_layers], return_contacts=False)
        reps = out["representations"][num_layers]
        pooled = []
        for k, (_, s) in enumerate(labeled):
            slen = min(len(s), max_len)
            pooled.append(reps[k, 1:slen+1].mean(0))
        return torch.stack(pooled)  # (B, D)

    def encode_triplet_batch(anc_idx, pos_idx, neg_idx):
        """
        Encode anchor, positive, negative in one batched forward pass.
        Returns (d_anc, d_pos, d_neg) each of shape (B, D), or None on failure.
        """
        B = len(anc_idx)
        # Build concatenated batch: [wt_anc, mut_anc, wt_pos, mut_pos, wt_neg, mut_neg]
        seqs = []
        valid = []
        for j in anc_idx:
            seqs.append(tr_wt[j] if tr_wt[j] else "A")
            seqs.append(tr_mut[j] if tr_mut[j] else "A")
            valid.append(tr_wt[j] is not None and tr_mut[j] is not None)
        for j in pos_idx:
            seqs.append(tr_wt[j] if tr_wt[j] else "A")
            seqs.append(tr_mut[j] if tr_mut[j] else "A")
        for j in neg_idx:
            seqs.append(tr_wt[j] if tr_wt[j] else "A")
            seqs.append(tr_mut[j] if tr_mut[j] else "A")

        labeled = [(str(i), s[:max_len]) for i, s in enumerate(seqs)]
        _, _, toks = batch_converter(labeled)
        out = fold_model(toks.to(device), repr_layers=[num_layers], return_contacts=False)
        reps = out["representations"][num_layers]

        pooled = []
        for k, (_, s) in enumerate(labeled):
            slen = min(len(s), max_len)
            pooled.append(reps[k, 1:slen+1].mean(0))

        # Split: each variant is 2 consecutive entries (wt, mut)
        def delta_block(start):
            ds = []
            for i in range(B):
                wi = start + i * 2
                mi = wi + 1
                ds.append(pooled[mi] - pooled[wi])
            return torch.stack(ds)

        d_anc = delta_block(0)
        d_pos = delta_block(B * 2)
        d_neg = delta_block(B * 4)
        return d_anc, d_pos, d_neg

    # Enable gradient checkpointing to save VRAM
    if hasattr(fold_model, "gradient_checkpointing_enable"):
        fold_model.gradient_checkpointing_enable()

    for epoch in range(max_epochs):
        fold_model.train()
        epoch_loss = 0.0
        n_steps = 0
        optimizer.zero_grad()

        for step, (anc_b, pos_b, neg_b) in enumerate(loader):
            anc_b = anc_b.numpy(); pos_b = pos_b.numpy(); neg_b = neg_b.numpy()

            try:
                d_a, d_p, d_n = encode_triplet_batch(anc_b, pos_b, neg_b)
            except torch.cuda.OutOfMemoryError:
                torch.cuda.empty_cache()
                continue

            loss = triplet_loss_fn(d_a, d_p, d_n) / grad_accum
            loss.backward()
            epoch_loss += loss.item() * grad_accum
            n_steps += 1

            if (step + 1) % grad_accum == 0:
                torch.nn.utils.clip_grad_norm_(fold_model.parameters(), 1.0)
                optimizer.step()
                optimizer.zero_grad()
                torch.cuda.empty_cache()

        if n_steps % grad_accum != 0:
            torch.nn.utils.clip_grad_norm_(fold_model.parameters(), 1.0)
            optimizer.step()
            optimizer.zero_grad()

        # Val loss on capped subset
        fold_model.eval()
        val_losses = []
        val_cap = min(50, len(anc_val))
        with torch.no_grad():
            for i in range(0, val_cap, batch_size):
                ai = anc_val[i:i+batch_size]; pi = pos_val[i:i+batch_size]; ni = neg_val[i:i+batch_size]
                try:
                    d_a, d_p, d_n = encode_triplet_batch(ai, pi, ni)
                    val_losses.append(triplet_loss_fn(d_a, d_p, d_n).item())
                except torch.cuda.OutOfMemoryError:
                    torch.cuda.empty_cache()
        val_loss = float(np.mean(val_losses)) if val_losses else float("inf")

        avg_train = epoch_loss / max(n_steps, 1)
        print(f"    Epoch {epoch+1:3d}  train_loss={avg_train:.4f}  val_loss={val_loss:.4f}")

        if val_loss < best_val_loss - 1e-4:
            best_val_loss = val_loss
            patience_count = 0
            best_state = {k: v.clone() for k, v in fold_model.state_dict().items()}
        else:
            patience_count += 1
            if patience_count >= patience:
                print(f"    Early stop at epoch {epoch+1}")
                break

    if best_state is not None:
        fold_model.load_state_dict(best_state)

    return fold_model

# ---------------------------------------------------------------------------
# Evaluation: extract embeddings from fine-tuned model, then linear + k-NN
# ---------------------------------------------------------------------------

def extract_delta_embeddings(fold_model, alphabet, wt_seqs, mut_seqs, indices,
                              device, embed_batch=8, max_len=1022):
    """Extract delta embeddings for given indices using the fine-tuned model."""
    import torch
    batch_converter = alphabet.get_batch_converter()
    fold_model.eval()

    num_layers = fold_model.num_layers

    all_deltas = []
    valid_flags = []

    for i in range(0, len(indices), embed_batch):
        batch_idx = indices[i:i+embed_batch]
        wt_batch  = []
        mut_batch = []
        batch_valid = []
        for j in batch_idx:
            if wt_seqs[j] is not None and mut_seqs[j] is not None:
                wt_batch.append((str(j), wt_seqs[j][:max_len]))
                mut_batch.append((str(j), mut_seqs[j][:max_len]))
                batch_valid.append(True)
            else:
                batch_valid.append(False)

        if not wt_batch:
            valid_flags.extend(batch_valid)
            continue

        # Concatenate wt and mut in one forward pass
        combined = wt_batch + mut_batch
        _, _, toks = batch_converter(combined)

        with torch.no_grad():
            out = fold_model(toks.to(device), repr_layers=[num_layers], return_contacts=False)
            reps = out["representations"][num_layers]

            n_valid = len(wt_batch)
            k_valid = 0
            for k, is_valid in enumerate(batch_valid):
                if is_valid:
                    wt_seq = wt_batch[k_valid][1]
                    mut_seq = mut_batch[k_valid][1]
                    wt_slen = min(len(wt_seq), max_len)
                    mut_slen = min(len(mut_seq), max_len)
                    wt_rep  = reps[k_valid, 1:wt_slen+1].mean(0)
                    mut_rep = reps[n_valid + k_valid, 1:mut_slen+1].mean(0)
                    all_deltas.append((mut_rep - wt_rep).cpu().float())
                    k_valid += 1
                valid_flags.append(is_valid)

    if not all_deltas:
        return np.zeros((0, 0), dtype=np.float32), np.array(valid_flags)

    delta_matrix = torch.stack(all_deltas).numpy()
    return delta_matrix, np.array(valid_flags)


def evaluate_embeddings(Z_train, Z_test, y_train, y_test, le):
    results = {}

    # k-NN
    k = min(10, len(Z_train) - 1)
    knn = KNeighborsClassifier(n_neighbors=k, metric="cosine")
    knn.fit(Z_train, y_train)
    pred_knn = knn.predict(Z_test)
    results["knn_macro_f1"] = float(f1_score(y_test, pred_knn, average="macro", zero_division=0))

    raw_proba = knn.predict_proba(Z_test)
    train_cls = list(knn.classes_)
    all_cls   = list(le.classes_)
    proba = np.zeros((len(Z_test), len(all_cls)), dtype=np.float32)
    for ti, ci in enumerate(train_cls):
        ai = all_cls.index(le.classes_[ci])
        proba[:, ai] = raw_proba[:, ti]
    for ai, cls_str in enumerate(all_cls):
        ci = le.transform([cls_str])[0]
        y_bin = (y_test == ci).astype(int)
        if y_bin.sum() > 0 and (1 - y_bin).sum() > 0 and proba[:, ai].std() > 0:
            results[f"knn_auroc_{cls_str}"] = float(roc_auc_score(y_bin, proba[:, ai]))

    # Linear probe
    mu = Z_train.mean(0); std = Z_train.std(0) + 1e-8
    Z_tr_n = (Z_train - mu) / std
    Z_te_n = (Z_test  - mu) / std
    lr_clf = LogisticRegression(max_iter=1000, C=0.1, multi_class="multinomial", solver="lbfgs")
    try:
        lr_clf.fit(Z_tr_n, y_train)
        pred_lr = lr_clf.predict(Z_te_n)
        results["linear_macro_f1"] = float(f1_score(y_test, pred_lr, average="macro", zero_division=0))
        lr_proba = lr_clf.predict_proba(Z_te_n)
        for ai, cls_str in enumerate(all_cls):
            if ai < lr_proba.shape[1]:
                ci = le.transform([cls_str])[0]
                y_bin = (y_test == ci).astype(int)
                if y_bin.sum() > 0 and (1 - y_bin).sum() > 0:
                    results[f"linear_auroc_{cls_str}"] = float(
                        roc_auc_score(y_bin, lr_proba[:, ai]))
    except Exception as e:
        print(f"    Linear probe failed: {e}")

    return results

# ---------------------------------------------------------------------------
# Main CV loop
# ---------------------------------------------------------------------------

def run_cv_lora(variants, labels, genes, gene_pfam, pfam_map, sequences,
                le, splits, split_name, base_model, alphabet, device,
                lora_rank=8, lora_alpha=16, lora_dropout=0.1,
                lr=2e-4, max_epochs=20, patience=5,
                batch_size=16, grad_accum=4, margin=1.0,
                embed_batch=8, max_len=1022, seed=42):

    wt_seqs, mut_seqs = get_variant_sequences(variants, sequences)
    y = le.transform(labels)

    fold_results = []

    for fold_i, (train_idx, test_idx) in enumerate(splits):
        labels_tr   = labels[train_idx]
        gene_pfam_tr = gene_pfam[train_idx]
        y_tr = y[train_idx]
        y_te = y[test_idx]

        if len(set(y_tr)) < 2 or len(set(y_te)) < 2:
            print(f"  Fold {fold_i+1}: skipped (missing class in train or test)")
            continue

        print(f"\n  Fold {fold_i+1}/{len(splits)} [{split_name}]  "
              f"train={len(train_idx)} test={len(test_idx)}")
        print(f"    train classes: {dict(Counter(labels_tr))}")

        # Fine-tune LoRA on training fold
        fold_model = train_lora_fold(
            base_model, alphabet, wt_seqs, mut_seqs,
            labels, gene_pfam, train_idx, le, device,
            lora_rank=lora_rank, lora_alpha=lora_alpha, lora_dropout=lora_dropout,
            lr=lr, max_epochs=max_epochs, patience=patience,
            batch_size=batch_size, grad_accum=grad_accum, margin=margin,
            max_pairs_per_anchor=8, embed_batch=embed_batch, max_len=max_len,
            seed=seed + fold_i
        )

        # Extract delta embeddings from fine-tuned model
        print("    Extracting train embeddings...")
        Z_tr, tr_valid = extract_delta_embeddings(
            fold_model, alphabet, wt_seqs, mut_seqs,
            train_idx, device, embed_batch, max_len
        )
        print("    Extracting test embeddings...")
        Z_te, te_valid = extract_delta_embeddings(
            fold_model, alphabet, wt_seqs, mut_seqs,
            test_idx, device, embed_batch, max_len
        )

        # Filter to valid only
        y_tr_valid = y_tr[tr_valid]
        y_te_valid = y_te[te_valid]

        if len(Z_tr) < 10 or len(Z_te) < 5:
            print(f"    Fold {fold_i+1}: too few valid variants after embedding, skipping")
            continue

        fm = evaluate_embeddings(Z_tr, Z_te, y_tr_valid, y_te_valid, le)
        fold_results.append(fm)

        print(f"    k-NN  macro_f1={fm.get('knn_macro_f1', float('nan')):.3f}  "
              f"GOF={fm.get('knn_auroc_GOF', float('nan')):.3f}  "
              f"DN={fm.get('knn_auroc_DN', float('nan')):.3f}  "
              f"LOF={fm.get('knn_auroc_LOF', float('nan')):.3f}")
        print(f"    Linear macro_f1={fm.get('linear_macro_f1', float('nan')):.3f}  "
              f"GOF={fm.get('linear_auroc_GOF', float('nan')):.3f}  "
              f"DN={fm.get('linear_auroc_DN', float('nan')):.3f}  "
              f"LOF={fm.get('linear_auroc_LOF', float('nan')):.3f}")

        # Free GPU memory
        del fold_model
        import torch
        torch.cuda.empty_cache()

    def agg(fold_list):
        if not fold_list:
            return {"error": "no folds"}
        all_keys = set().union(*[set(f.keys()) for f in fold_list])
        out = {}
        for k in all_keys:
            vals = [f[k] for f in fold_list if k in f and not np.isnan(f[k])]
            if vals:
                out[f"{k}_mean"] = float(np.mean(vals))
                out[f"{k}_std"]  = float(np.std(vals))
        out["n_folds"] = len(fold_list)
        return out

    return agg(fold_results)

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_dir",    default="../data")
    parser.add_argument("--out_dir",     default="../results/20260524_baseline_run/run_0")
    parser.add_argument("--model_name",  default="esm2_t33_650M_UR50D",
                        choices=["esm2_t33_650M_UR50D", "esm2_t36_3B_UR50D",
                                 "esm2_t30_150M_UR50D"])
    parser.add_argument("--lora_rank",   type=int,   default=8)
    parser.add_argument("--lora_alpha",  type=int,   default=16)
    parser.add_argument("--lora_dropout", type=float, default=0.1)
    parser.add_argument("--lr",          type=float, default=2e-4)
    parser.add_argument("--max_epochs",  type=int,   default=20)
    parser.add_argument("--patience",    type=int,   default=5)
    parser.add_argument("--batch_size",  type=int,   default=16,
                        help="Triplet batch size per gradient step")
    parser.add_argument("--grad_accum",  type=int,   default=4,
                        help="Gradient accumulation steps")
    parser.add_argument("--embed_batch", type=int,   default=8,
                        help="Sequences per forward pass for embedding extraction")
    parser.add_argument("--max_len",     type=int,   default=1022,
                        help="Max sequence length (ESM-2 supports up to 1022)")
    parser.add_argument("--margin",      type=float, default=1.0)
    parser.add_argument("--n_folds",     type=int,   default=5)
    parser.add_argument("--seed",        type=int,   default=0)
    parser.add_argument("--merged",      action="store_true",
                        help="Use merged dataset (19100 variants)")
    parser.add_argument("--gene_split_only", action="store_true",
                        help="Only run gene-split CV (faster, for debugging)")
    args = parser.parse_args()

    import torch
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")
    if device.type == "cuda":
        print(f"GPU: {torch.cuda.get_device_name(0)}")
        print(f"VRAM: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")

    print("\n=== Loading data ===")
    variants, labels, genes, sequences = load_data(args.data_dir, merged=args.merged)
    gene_pfam, pfam_map = load_pfam(args.data_dir, genes)

    le = LabelEncoder()
    le.fit(labels)
    print(f"Classes: {list(le.classes_)}")

    print("\n=== Building CV splits ===")
    gene_splits = gene_split_cv(genes, n_folds=args.n_folds, seed=args.seed)
    fam_splits  = family_split_cv(genes, gene_pfam, pfam_map,
                                   n_folds=args.n_folds, seed=args.seed)
    print(f"Gene-split: {len(gene_splits)} folds | Family-split: {len(fam_splits)} folds")

    print(f"\n=== Loading ESM-2: {args.model_name} ===")
    base_model, alphabet = load_esm2(args.model_name, device)

    model_tag = args.model_name.replace("esm2_", "").replace("_UR50D", "")
    dataset_tag = "merged" if args.merged else "geras"

    results = {
        "model": args.model_name,
        "dataset": dataset_tag,
        "lora_rank": args.lora_rank,
        "lora_alpha": args.lora_alpha,
        "lora_dropout": args.lora_dropout,
        "lr": args.lr,
        "max_epochs": args.max_epochs,
        "patience": args.patience,
        "batch_size": args.batch_size,
        "grad_accum": args.grad_accum,
        "margin": args.margin,
        "seed": args.seed,
        "n_folds": args.n_folds,
    }

    shared_kwargs = dict(
        le=le, base_model=base_model, alphabet=alphabet, device=device,
        lora_rank=args.lora_rank, lora_alpha=args.lora_alpha,
        lora_dropout=args.lora_dropout, lr=args.lr,
        max_epochs=args.max_epochs, patience=args.patience,
        batch_size=args.batch_size, grad_accum=args.grad_accum,
        margin=args.margin, embed_batch=args.embed_batch,
        max_len=args.max_len, seed=args.seed,
    )

    print("\n" + "=" * 60)
    print("GENE-SPLIT CV")
    print("=" * 60)
    gene_results = run_cv_lora(
        variants, labels, genes, gene_pfam, pfam_map, sequences,
        splits=gene_splits, split_name="gene-split", **shared_kwargs
    )
    results["gene_split"] = gene_results

    if not args.gene_split_only:
        print("\n" + "=" * 60)
        print("FAMILY-SPLIT CV")
        print("=" * 60)
        fam_results = run_cv_lora(
            variants, labels, genes, gene_pfam, pfam_map, sequences,
            splits=fam_splits, split_name="family-split", **shared_kwargs
        )
        results["family_split"] = fam_results
    else:
        results["family_split"] = {"skipped": True}

    # Summary
    print("\n" + "=" * 60)
    print("HEADLINE SUMMARY")
    print("=" * 60)
    baselines = {
        "MLP (result_7) family-split F1":         0.364,
        "Contrastive proj (result_9) family-split F1": 0.397,
    }
    for split_name, split_key in [("Gene-split", "gene_split"), ("Family-split", "family_split")]:
        r = results.get(split_key, {})
        if r.get("skipped"):
            continue
        knn_f1    = r.get("knn_macro_f1_mean",    float("nan"))
        linear_f1 = r.get("linear_macro_f1_mean", float("nan"))
        print(f"\n{split_name}:")
        print(f"  k-NN    macro_f1 = {knn_f1:.3f} ± {r.get('knn_macro_f1_std', float('nan')):.3f}"
              f"  GOF={r.get('knn_auroc_GOF_mean', float('nan')):.3f}"
              f"  DN={r.get('knn_auroc_DN_mean', float('nan')):.3f}"
              f"  LOF={r.get('knn_auroc_LOF_mean', float('nan')):.3f}")
        print(f"  Linear  macro_f1 = {linear_f1:.3f} ± {r.get('linear_macro_f1_std', float('nan')):.3f}"
              f"  GOF={r.get('linear_auroc_GOF_mean', float('nan')):.3f}"
              f"  DN={r.get('linear_auroc_DN_mean', float('nan')):.3f}"
              f"  LOF={r.get('linear_auroc_LOF_mean', float('nan')):.3f}")

    if not args.gene_split_only:
        fam_knn_f1 = results["family_split"].get("knn_macro_f1_mean", float("nan"))
        print(f"\nVs. baselines (family-split k-NN F1 = {fam_knn_f1:.3f}):")
        for name, val in baselines.items():
            delta = fam_knn_f1 - val
            symbol = "✓" if delta > 0.03 else ("~" if delta > 0 else "✗")
            print(f"  {symbol} vs {name} ({val:.3f}): Δ = {delta:+.3f}")

    os.makedirs(args.out_dir, exist_ok=True)
    out_path = os.path.join(
        args.out_dir,
        f"lora_contrastive_{model_tag}_{dataset_tag}_seed{args.seed}.json"
    )
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults written to {out_path}")


if __name__ == "__main__":
    main()
