import argparse
import json
import os
import time

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset

# ---------------------------------------------------------------------------
# ESM-2 functional classification: DNA repair vs Tumor Suppressor
# Direct comparison to evo2_function — same genes, same MLP, different embeddings.
# Input: protein sequences (translated CDS) rather than DNA sequences.
# Model: ESM-2 (protein language model) vs Evo2 (DNA language model).
# ---------------------------------------------------------------------------

DNA_REPAIR_PATHWAYS = ["hsa03430", "hsa03420", "hsa03440", "hsa03410"]

TSG_SEED_SYMBOLS = [
    "TP53", "RB1", "APC", "PTEN", "VHL", "BRCA1", "BRCA2", "NF1", "NF2",
    "WT1", "CDKN2A", "CDKN1B", "CDKN1C", "MLH1", "MSH2", "MSH6", "PMS2",
    "STK11", "TSC1", "TSC2", "SMAD4", "SMAD2", "RET", "MEN1", "PTCH1",
    "SUFU", "BAP1", "PALB2", "ATM", "CHEK2", "CDH1", "RUNX1", "TET2",
    "DNMT3A", "ASXL1", "EZH2", "KDM6A", "ARID1A", "ARID1B", "PBRM1",
    "SETD2", "KDM5C", "FBXW7", "PPP2R1A", "PIK3R1", "INPP4B",
]

ENSEMBL_REST = "https://rest.ensembl.org"
ESM2_MODEL = "esm2_t33_650M_UR50D"  # 650M — strong representation, fits on A100


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


def get_dna_repair_symbols():
    seen_ids = set()
    symbols = []
    for pathway in DNA_REPAIR_PATHWAYS:
        try:
            kegg_ids = fetch_kegg_genes(pathway)
            for kid in kegg_ids:
                if kid in seen_ids:
                    continue
                seen_ids.add(kid)
                sym = kegg_id_to_symbol(kid)
                if sym:
                    symbols.append(sym.upper())
                time.sleep(0.2)
        except Exception as e:
            print(f"Warning: could not fetch {pathway}: {e}")
    return list(set(symbols))


# ---------------------------------------------------------------------------
# Ensembl — fetch protein sequence (translated CDS)
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


def fetch_protein_sequence(ensembl_gene_id, retries=3, delay=1.0):
    """Fetch the canonical protein sequence (longest translation) for a gene."""
    import urllib.request
    url = (
        f"{ENSEMBL_REST}/sequence/id/{ensembl_gene_id}"
        f"?content-type=application/json&type=protein&multiple_sequences=1"
    )
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                data = json.loads(resp.read())
            if isinstance(data, list):
                # Pick longest protein sequence
                return max(data, key=lambda x: len(x["seq"]))["seq"]
            return data["seq"]
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(delay)
    return None


def get_gene_protein_sequence(symbol):
    """Return canonical protein sequence for a gene symbol."""
    info = symbol_to_ensembl(symbol)
    if info is None:
        return None
    ensembl_id = info["id"]
    time.sleep(0.3)
    protein_seq = fetch_protein_sequence(ensembl_id)
    time.sleep(0.3)
    if protein_seq is None:
        return None
    return protein_seq.upper()


# ---------------------------------------------------------------------------
# ESM-2 embedding extraction
# ---------------------------------------------------------------------------

def get_esm2_embeddings(sequences, model_name=ESM2_MODEL, device="cuda", batch_size=8):
    import esm
    model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
    model = model.to(device).eval()
    batch_converter = alphabet.get_batch_converter()

    all_embeddings = []
    for i in range(0, len(sequences), batch_size):
        batch_seqs = sequences[i: i + batch_size]
        # Truncate long sequences to avoid OOM (ESM-2 handles up to ~1022 tokens)
        batch_seqs_trunc = [s[:1022] for s in batch_seqs]
        data = [(f"seq{j}", s) for j, s in enumerate(batch_seqs_trunc)]
        _, _, tokens = batch_converter(data)
        tokens = tokens.to(device)
        with torch.no_grad():
            results = model(tokens, repr_layers=[model.num_layers])
        reps = results["representations"][model.num_layers]
        for k, seq in enumerate(batch_seqs_trunc):
            # Mean pool over residues (excluding BOS/EOS tokens)
            emb = reps[k, 1: len(seq) + 1].mean(0).cpu().float().numpy()
            all_embeddings.append(emb)
        if i % 50 == 0:
            print(f"  Embedded {i + len(batch_seqs)}/{len(sequences)} sequences")
    return np.stack(all_embeddings)


# ---------------------------------------------------------------------------
# MLP classifier (identical to evo2_function for fair comparison)
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

    # DNA repair genes (label=0)
    print("Fetching DNA repair gene list from KEGG...")
    repair_symbols = get_dna_repair_symbols()
    print(f"  Found {len(repair_symbols)} DNA repair genes")
    for sym in repair_symbols:
        print(f"  Fetching protein sequence for {sym}...")
        seq = get_gene_protein_sequence(sym)
        if seq:
            dataset.append({"symbol": sym, "label": 0, "label_name": "DNA_repair", "seq": seq})

    # Tumor suppressor genes (label=1)
    print(f"Fetching sequences for {len(TSG_SEED_SYMBOLS)} tumor suppressor genes...")
    for sym in TSG_SEED_SYMBOLS:
        print(f"  Fetching protein sequence for {sym}...")
        seq = get_gene_protein_sequence(sym)
        if seq:
            dataset.append({"symbol": sym, "label": 1, "label_name": "tumor_suppressor", "seq": seq})

    with open(cache_path, "w") as f:
        json.dump(dataset, f)

    n0 = sum(1 for d in dataset if d["label"] == 0)
    n1 = sum(1 for d in dataset if d["label"] == 1)
    print(f"Dataset: {n0} DNA repair, {n1} tumor suppressor genes")
    return dataset


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

def run(out_dir, seed, model_name=ESM2_MODEL, hidden_dims=(256, 128),
        dropout=0.2, lr=1e-3, epochs=10, batch_size=16, train_frac=0.8):
    os.makedirs(out_dir, exist_ok=True)
    torch.manual_seed(seed)
    np.random.seed(seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    cache_dir = os.path.join(out_dir, "data")
    dataset = build_dataset(cache_dir)
    sequences = [d["seq"] for d in dataset]
    labels = np.array([d["label"] for d in dataset], dtype=np.float32)
    symbols = [d["symbol"] for d in dataset]

    print(f"Total genes: {len(sequences)} | DNA repair: {int((labels==0).sum())} | TSG: {int((labels==1).sum())}")

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
        print(f"Extracting ESM-2 embeddings ({model_name})...")
        embeddings = get_esm2_embeddings(sequences, model_name, str(device))
        np.save(emb_cache, embeddings)

    input_dim = embeddings.shape[1]
    print(f"Embedding shape: {embeddings.shape}")

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
            print(f"  epoch {epoch}: loss={train_loss:.4f} auroc={val_metrics['auroc']:.3f} acc={val_metrics['accuracy']:.3f}")

    # Unsupervised geometry metrics for direct comparison with evo2_function
    from sklearn.metrics import silhouette_score
    from sklearn.cluster import KMeans
    from sklearn.metrics import adjusted_rand_score
    kmeans = KMeans(n_clusters=2, n_init=10, random_state=seed)
    cluster_labels = kmeans.fit_predict(embeddings)
    ari = float(adjusted_rand_score(labels, cluster_labels))
    sil = float(silhouette_score(embeddings, labels, metric="cosine"))
    from scipy.spatial.distance import cdist
    repair_embs = embeddings[labels == 0]
    tsg_embs = embeddings[labels == 1]
    within_repair = float(np.mean(cdist(repair_embs, repair_embs, metric="cosine")[np.triu_indices(len(repair_embs), k=1)]))
    within_tsg = float(np.mean(cdist(tsg_embs, tsg_embs, metric="cosine")[np.triu_indices(len(tsg_embs), k=1)]))
    between = float(np.mean(cdist(repair_embs, tsg_embs, metric="cosine")))
    within_mean = (within_repair + within_tsg) / 2
    distance_ratio = within_mean / between

    final_info = {
        "final_train_loss": train_log[-1]["train_loss"],
        "final_val_accuracy": val_log[-1]["accuracy"],
        "final_val_auroc": val_log[-1]["auroc"],
        "final_val_f1": val_log[-1]["f1"],
        "best_val_auroc": float(max(v["auroc"] for v in val_log)),
        "best_val_accuracy": float(max(v["accuracy"] for v in val_log)),
        "kmeans_ari": ari,
        "silhouette_score": sil,
        "within_class_distance": within_mean,
        "between_class_distance": between,
        "distance_ratio": distance_ratio,
        "n_train": int(n_train),
        "n_test": int(len(sequences) - n_train),
        "n_dna_repair": int((labels == 0).sum()),
        "n_tsg": int((labels == 1).sum()),
        "input_dim": int(input_dim),
        "model": model_name,
    }
    print(final_info)

    with open(os.path.join(out_dir, f"final_info_seed{seed}.json"), "w") as f:
        json.dump(final_info, f)

    return final_info, train_log, val_log


parser = argparse.ArgumentParser()
parser.add_argument("--out_dir", type=str, default="run_0")
args = parser.parse_args()

if __name__ == "__main__":
    seeds = [0]
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
