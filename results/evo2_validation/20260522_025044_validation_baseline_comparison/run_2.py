import argparse
import json
import os
import time

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset

# ---------------------------------------------------------------------------
# Validation experiment: Glycolysis vs T-cell receptor signaling
# Completely different functional contrast from evo2_function (DNA repair vs TSG)
# Purpose: verify that the 0.78 AUROC baseline generalises
# ---------------------------------------------------------------------------

GLYCOLYSIS_PATHWAY = "hsa00010"
TCELL_PATHWAY = "hsa04660"

ENSEMBL_REST = "https://rest.ensembl.org"
PROMOTER_BP = 2000
EVO2_LAYER = "blocks.28.mlp.l3"


# ---------------------------------------------------------------------------
# KEGG
# ---------------------------------------------------------------------------

def fetch_kegg_genes(pathway_id, retries=3, delay=1.0):
    import urllib.request
    url = f"https://rest.kegg.jp/link/hsa/{pathway_id}"
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                text = resp.read().decode()
            return [line.split("\t")[1].strip() for line in text.strip().split("\n") if line]
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(delay)
            else:
                raise e


def kegg_id_to_symbol(kegg_id, retries=3, delay=0.5):
    import urllib.request
    url = f"https://rest.kegg.jp/list/{kegg_id}"
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=15) as resp:
                text = resp.read().decode().strip()
            if not text:
                return None
            parts = text.split("\t")
            if len(parts) < 2:
                return None
            return parts[1].split(",")[0].split(";")[0].strip().upper()
        except Exception:
            if attempt < retries - 1:
                time.sleep(delay)
    return None


def get_pathway_symbols(pathway_id):
    kegg_ids = fetch_kegg_genes(pathway_id)
    symbols = []
    for kid in kegg_ids:
        sym = kegg_id_to_symbol(kid)
        if sym:
            symbols.append(sym)
        time.sleep(0.2)
    return list(set(symbols))


# ---------------------------------------------------------------------------
# Ensembl
# ---------------------------------------------------------------------------

def symbol_to_ensembl(symbol, retries=3, delay=1.0):
    import urllib.request
    url = f"{ENSEMBL_REST}/lookup/symbol/homo_sapiens/{symbol}?content-type=application/json"
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                return json.loads(resp.read())
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(delay)
    return None


def fetch_promoter_sequence(chrom, tss, strand, bp=PROMOTER_BP, retries=3, delay=1.0):
    import urllib.request
    if strand == 1:
        start = max(1, tss - bp)
        end = tss
        strand_param = 1
    else:
        start = tss
        end = tss + bp
        strand_param = -1
    url = (f"{ENSEMBL_REST}/sequence/region/human/{chrom}:{start}..{end}:{strand_param}"
           f"?content-type=application/json")
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                return json.loads(resp.read())["seq"]
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(delay)
    return None


def fetch_cds_sequence(ensembl_gene_id, retries=3, delay=1.0):
    import urllib.request
    url = (f"{ENSEMBL_REST}/sequence/id/{ensembl_gene_id}"
           f"?content-type=application/json&type=cds&multiple_sequences=1")
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                data = json.loads(resp.read())
            if isinstance(data, list):
                return max(data, key=lambda x: len(x["seq"]))["seq"]
            return data["seq"]
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(delay)
    return None


def get_gene_sequence(symbol):
    info = symbol_to_ensembl(symbol)
    if info is None:
        return None
    chrom = info["seq_region_name"]
    strand = info["strand"]
    tss = info["start"] if strand == 1 else info["end"]
    ensembl_id = info["id"]
    promoter = fetch_promoter_sequence(chrom, tss, strand)
    time.sleep(0.3)
    cds = fetch_cds_sequence(ensembl_id)
    time.sleep(0.3)
    if promoter is None or cds is None:
        return None
    return (promoter + cds).upper()


# ---------------------------------------------------------------------------
# Evo2
# ---------------------------------------------------------------------------

def load_evo2(model_name="evo2_7b", device="cuda"):
    from evo2 import Evo2
    return Evo2(model_name)


def get_evo2_embeddings(sequences, model, device, max_len=8192, batch_size=4):
    embeddings = []
    for i in range(0, len(sequences), batch_size):
        batch = sequences[i: i + batch_size]
        for seq in batch:
            seq = seq[:max_len]
            input_ids = torch.tensor(
                model.tokenizer.tokenize(seq), dtype=torch.int
            ).unsqueeze(0).to(device)
            with torch.no_grad():
                _, emb_dict = model(input_ids, return_embeddings=True, layer_names=[EVO2_LAYER])
            emb = emb_dict[EVO2_LAYER].squeeze(0).mean(0).float().cpu().numpy()
            embeddings.append(emb)
    return np.stack(embeddings)


# ---------------------------------------------------------------------------
# MLP
# ---------------------------------------------------------------------------

class MLP(nn.Module):
    def __init__(self, input_dim, hidden_dims, dropout=0.2):
        super().__init__()
        layers = []
        prev = input_dim
        for h in hidden_dims:
            layers += [nn.Linear(prev, h), nn.ReLU(), nn.Dropout(dropout)]
            prev = h
        layers.append(nn.Linear(prev, 1))
        self.net = nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x).squeeze(-1)


def train_epoch(model, loader, optimizer, device):
    model.train()
    total_loss = 0.0
    for X_batch, y_batch in loader:
        X_batch, y_batch = X_batch.to(device), y_batch.to(device)
        optimizer.zero_grad()
        loss = nn.functional.binary_cross_entropy_with_logits(model(X_batch), y_batch)
        loss.backward()
        optimizer.step()
        total_loss += loss.item() * len(y_batch)
    return total_loss / len(loader.dataset)


def evaluate(model, loader, device):
    model.eval()
    all_logits, all_labels = [], []
    with torch.no_grad():
        for X_batch, y_batch in loader:
            all_logits.append(model(X_batch.to(device)).cpu().numpy())
            all_labels.append(y_batch.numpy())
    logits = np.concatenate(all_logits)
    labels = np.concatenate(all_labels)
    probs = 1 / (1 + np.exp(-logits))
    preds = (probs >= 0.5).astype(float)
    from sklearn.metrics import roc_auc_score, f1_score
    auroc = float(roc_auc_score(labels, probs)) if len(np.unique(labels)) > 1 else 0.5
    f1 = float(f1_score(labels, preds, zero_division=0))
    return {"accuracy": float(np.mean(preds == labels)), "auroc": auroc, "f1": f1}


# ---------------------------------------------------------------------------
# Dataset
# ---------------------------------------------------------------------------

def build_dataset(cache_dir, force_refetch=False):
    cache_path = os.path.join(cache_dir, "dataset.json")
    if os.path.exists(cache_path) and not force_refetch:
        with open(cache_path) as f:
            return json.load(f)

    # Fall back to baseline template cache if available
    script_dir = os.path.dirname(os.path.abspath(__file__))
    baseline_cache = os.path.join(script_dir, "run_0", "data", "dataset.json")
    if os.path.exists(baseline_cache) and not force_refetch:
        print(f"Using baseline cache from {baseline_cache}")
        os.makedirs(cache_dir, exist_ok=True)
        import shutil
        shutil.copy(baseline_cache, cache_path)
        with open(cache_path) as f:
            return json.load(f)

    os.makedirs(cache_dir, exist_ok=True)
    dataset = []

    # Glycolysis genes (label=0)
    print("Fetching glycolysis gene list from KEGG...")
    glycolysis_symbols = get_pathway_symbols(GLYCOLYSIS_PATHWAY)
    print(f"  Found {len(glycolysis_symbols)} glycolysis genes")
    for sym in glycolysis_symbols:
        print(f"  Fetching sequence for {sym}...")
        seq = get_gene_sequence(sym)
        if seq:
            dataset.append({"symbol": sym, "label": 0, "label_name": "glycolysis", "seq": seq})

    # T-cell receptor signaling genes (label=1)
    print("Fetching T-cell receptor signaling gene list from KEGG...")
    tcell_symbols = get_pathway_symbols(TCELL_PATHWAY)
    print(f"  Found {len(tcell_symbols)} T-cell signaling genes")
    for sym in tcell_symbols:
        print(f"  Fetching sequence for {sym}...")
        seq = get_gene_sequence(sym)
        if seq:
            dataset.append({"symbol": sym, "label": 1, "label_name": "tcell_signaling", "seq": seq})

    with open(cache_path, "w") as f:
        json.dump(dataset, f)

    n0 = sum(1 for d in dataset if d["label"] == 0)
    n1 = sum(1 for d in dataset if d["label"] == 1)
    print(f"Dataset: {n0} glycolysis, {n1} T-cell signaling genes")
    return dataset


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

def run(out_dir, seed, model_name="evo2_7b", hidden_dims=(256, 128),
        dropout=0.2, lr=1e-3, epochs=10, batch_size=16, train_frac=0.8):
    os.makedirs(out_dir, exist_ok=True)
    torch.manual_seed(seed)
    np.random.seed(seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    cache_dir = os.path.join(out_dir, "data")
    dataset = build_dataset(cache_dir)
    sequences = [d["seq"] for d in dataset]
    labels = np.array([d["label"] for d in dataset], dtype=np.float32)

    print(f"Total genes: {len(sequences)} | Glycolysis: {int((labels==0).sum())} | T-cell: {int((labels==1).sum())}")

    emb_cache = os.path.join(out_dir, "data", f"embeddings_{model_name}.npy")
    script_dir = os.path.dirname(os.path.abspath(__file__))
    baseline_emb_cache = os.path.join(script_dir, "run_0", "data", f"embeddings_{model_name}.npy")
    if os.path.exists(emb_cache):
        print("Loading cached embeddings...")
        embeddings = np.load(emb_cache)
    elif os.path.exists(baseline_emb_cache):
        print(f"Using baseline embeddings cache from {baseline_emb_cache}")
        import shutil
        shutil.copy(baseline_emb_cache, emb_cache)
        embeddings = np.load(emb_cache)
    else:
        print(f"Extracting Evo2 embeddings ({model_name})...")
        evo2 = load_evo2(model_name, str(device))
        embeddings = get_evo2_embeddings(sequences, evo2, device)
        np.save(emb_cache, embeddings)
        del evo2
        torch.cuda.empty_cache()

    input_dim = embeddings.shape[1]
    
    # Compute unsupervised centroid cosine similarity
    glycolysis_mask = labels == 0
    tcell_mask = labels == 1
    glycolysis_centroid = embeddings[glycolysis_mask].mean(axis=0)
    tcell_centroid = embeddings[tcell_mask].mean(axis=0)
    
    # Normalize centroids
    glycolysis_centroid_norm = glycolysis_centroid / np.linalg.norm(glycolysis_centroid)
    tcell_centroid_norm = tcell_centroid / np.linalg.norm(tcell_centroid)
    
    # Compute cosine similarity
    centroid_cosine_sim = float(np.dot(glycolysis_centroid_norm, tcell_centroid_norm))
    centroid_cosine_dist = 1.0 - centroid_cosine_sim
    
    print(f"Unsupervised analysis:")
    print(f"  Centroid cosine similarity: {centroid_cosine_sim:.4f}")
    print(f"  Centroid cosine distance: {centroid_cosine_dist:.4f}")
    
    # Compute within-class and between-class distances
    glycolysis_embs = embeddings[glycolysis_mask]
    tcell_embs = embeddings[tcell_mask]
    
    # Within-class average cosine similarity
    def avg_cosine_sim(embs):
        if len(embs) < 2:
            return 0.0
        embs_norm = embs / np.linalg.norm(embs, axis=1, keepdims=True)
        sim_matrix = embs_norm @ embs_norm.T
        # Exclude diagonal
        mask = ~np.eye(len(embs), dtype=bool)
        return float(sim_matrix[mask].mean())
    
    within_glycolysis_sim = avg_cosine_sim(glycolysis_embs)
    within_tcell_sim = avg_cosine_sim(tcell_embs)
    
    print(f"  Within-class cosine similarity (glycolysis): {within_glycolysis_sim:.4f}")
    print(f"  Within-class cosine similarity (T-cell): {within_tcell_sim:.4f}")
    rng = np.random.RandomState(seed)
    idx = rng.permutation(len(sequences))
    n_train = int(len(sequences) * train_frac)
    train_idx, test_idx = idx[:n_train], idx[n_train:]

    X_train = torch.tensor(embeddings[train_idx])
    y_train = torch.tensor(labels[train_idx])
    X_test = torch.tensor(embeddings[test_idx])
    y_test = torch.tensor(labels[test_idx])

    train_loader = DataLoader(TensorDataset(X_train, y_train), batch_size=batch_size, shuffle=True)
    test_loader = DataLoader(TensorDataset(X_test, y_test), batch_size=batch_size)

    model = MLP(input_dim, list(hidden_dims), dropout).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs)

    train_log, val_log = [], []
    for epoch in range(1, epochs + 1):
        train_loss = train_epoch(model, train_loader, optimizer, device)
        val_metrics = evaluate(model, test_loader, device)
        scheduler.step()
        train_log.append({"epoch": epoch, "train_loss": train_loss})
        val_log.append({"epoch": epoch, **val_metrics})
        if epoch % 5 == 0:
            print(f"  epoch {epoch}: loss={train_loss:.4f} auroc={val_metrics['auroc']:.3f}")

    final_info = {
        "final_train_loss": train_log[-1]["train_loss"],
        "final_val_accuracy": val_log[-1]["accuracy"],
        "final_val_auroc": val_log[-1]["auroc"],
        "final_val_f1": val_log[-1]["f1"],
        "best_val_auroc": float(max(v["auroc"] for v in val_log)),
        "best_val_accuracy": float(max(v["accuracy"] for v in val_log)),
        "n_train": int(n_train),
        "n_test": int(len(sequences) - n_train),
        "n_glycolysis": int((labels == 0).sum()),
        "n_tcell": int((labels == 1).sum()),
        "input_dim": int(input_dim),
        "centroid_cosine_similarity": centroid_cosine_sim,
        "centroid_cosine_distance": centroid_cosine_dist,
        "within_glycolysis_cosine_sim": within_glycolysis_sim,
        "within_tcell_cosine_sim": within_tcell_sim,
    }
    print(final_info)

    with open(os.path.join(out_dir, f"final_info_seed{seed}.json"), "w") as f:
        json.dump(final_info, f)

    return final_info, train_log, val_log


parser = argparse.ArgumentParser()
parser.add_argument("--out_dir", type=str, default="run_0")
args = parser.parse_args()

if __name__ == "__main__":
    seeds = [0, 1, 2, 3, 4]
    all_results = {}
    final_infos_list = []

    for seed in seeds:
        print(f"\n=== Seed {seed} ===")
        final_info, train_log, val_log = run(args.out_dir, seed)
        all_results[f"seed{seed}_final_info"] = final_info
        all_results[f"seed{seed}_train_log"] = train_log
        all_results[f"seed{seed}_val_log"] = val_log
        final_infos_list.append(final_info)

    numeric_keys = [k for k in final_infos_list[0] if isinstance(final_infos_list[0][k], (int, float))]
    final_infos = {
        "means": {k: float(np.mean([d[k] for d in final_infos_list])) for k in numeric_keys},
        "stderrs": {k: float(np.std([d[k] for d in final_infos_list]) / len(seeds)) for k in numeric_keys},
        "final_info_list": final_infos_list,
    }

    with open(os.path.join(args.out_dir, "final_info.json"), "w") as f:
        json.dump(final_infos, f)

    with open(os.path.join(args.out_dir, "all_results.npy"), "wb") as f:
        np.save(f, all_results)
