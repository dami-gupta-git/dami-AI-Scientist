import argparse
import json
import os
import time

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset

# ---------------------------------------------------------------------------
# ESM-2 validation: glycolysis vs T-cell signaling
# Direct comparison to evo2_validation — same genes, same MLP, protein embeddings.
# Gene set reused from evo2_validation dataset (symbols only, protein seqs fetched fresh).
# ---------------------------------------------------------------------------

ENSEMBL_REST = "https://rest.ensembl.org"
ESM2_MODEL = "esm2_t33_650M_UR50D"

# Relative path to evo2_validation dataset (source of gene symbols + labels)
EVO2_VALIDATION_DATASET = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "evo2_validation", "run_0", "data", "dataset.json"
)


# ---------------------------------------------------------------------------
# Ensembl — fetch protein sequence
# ---------------------------------------------------------------------------

def symbol_to_ensembl(symbol, retries=3, delay=1.0):
    import urllib.request
    url = f"{ENSEMBL_REST}/lookup/symbol/homo_sapiens/{symbol}?content-type=application/json"
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                return json.loads(resp.read())
        except Exception:
            if attempt < retries - 1:
                time.sleep(delay)
    return None


def fetch_protein_sequence(ensembl_gene_id, retries=3, delay=1.0):
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
                return max(data, key=lambda x: len(x["seq"]))["seq"]
            return data["seq"]
        except Exception:
            if attempt < retries - 1:
                time.sleep(delay)
    return None


def get_protein_sequence(symbol):
    info = symbol_to_ensembl(symbol)
    if info is None:
        return None
    ensembl_id = info["id"]
    time.sleep(0.3)
    seq = fetch_protein_sequence(ensembl_id)
    time.sleep(0.3)
    return seq.upper() if seq else None


# ---------------------------------------------------------------------------
# ESM-2 embeddings
# ---------------------------------------------------------------------------

def get_esm2_embeddings(sequences, model_name=ESM2_MODEL, device="cuda", batch_size=8):
    import esm
    model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
    model = model.to(device).eval()
    batch_converter = alphabet.get_batch_converter()

    all_embeddings = []
    for i in range(0, len(sequences), batch_size):
        batch_seqs = sequences[i: i + batch_size]
        batch_seqs_trunc = [s[:1022] for s in batch_seqs]
        data = [(f"seq{j}", s) for j, s in enumerate(batch_seqs_trunc)]
        _, _, tokens = batch_converter(data)
        tokens = tokens.to(device)
        with torch.no_grad():
            results = model(tokens, repr_layers=[model.num_layers])
        reps = results["representations"][model.num_layers]
        for k, seq in enumerate(batch_seqs_trunc):
            emb = reps[k, 1: len(seq) + 1].mean(0).cpu().float().numpy()
            all_embeddings.append(emb)
        if i % 50 == 0:
            print(f"  Embedded {i + len(batch_seqs)}/{len(sequences)} sequences")
    return np.stack(all_embeddings)


# ---------------------------------------------------------------------------
# MLP classifier (identical to esm2_function / evo2_validation)
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
# Dataset — reuse symbols from evo2_validation, fetch protein seqs fresh
# ---------------------------------------------------------------------------

def build_dataset(cache_dir):
    cache_path = os.path.join(cache_dir, "dataset.json")
    if os.path.exists(cache_path):
        print("Loading cached dataset...")
        with open(cache_path) as f:
            return json.load(f)

    print(f"Loading gene list from evo2_validation: {EVO2_VALIDATION_DATASET}")
    with open(EVO2_VALIDATION_DATASET) as f:
        evo2_data = json.load(f)

    os.makedirs(cache_dir, exist_ok=True)
    dataset = []
    for entry in evo2_data:
        symbol = entry["symbol"]
        label = entry["label"]
        label_name = entry["label_name"]
        print(f"  Fetching protein sequence for {symbol}...")
        seq = get_protein_sequence(symbol)
        if seq and len(seq) >= 10:
            dataset.append({"symbol": symbol, "label": label, "label_name": label_name, "seq": seq})
        else:
            print(f"    Warning: no protein sequence for {symbol}")

    with open(cache_path, "w") as f:
        json.dump(dataset, f)

    counts = {l: sum(1 for d in dataset if d["label_name"] == l) for l in set(d["label_name"] for d in dataset)}
    print(f"Dataset: {counts}")
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

    label_names = sorted(set(d["label_name"] for d in dataset))
    print(f"Total genes: {len(sequences)}")
    for ln in label_names:
        print(f"  {ln}: {sum(1 for d in dataset if d['label_name'] == ln)}")

    emb_cache = os.path.join(cache_dir, f"embeddings_{model_name}.npy")
    if os.path.exists(emb_cache):
        print("Loading cached embeddings...")
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

    # Unsupervised geometry metrics
    from sklearn.metrics import silhouette_score, adjusted_rand_score
    from sklearn.cluster import KMeans
    from scipy.spatial.distance import cdist
    kmeans = KMeans(n_clusters=2, n_init=10, random_state=seed)
    cluster_labels = kmeans.fit_predict(embeddings)
    ari = float(adjusted_rand_score(labels, cluster_labels))
    sil = float(silhouette_score(embeddings, labels, metric="cosine"))
    emb0 = embeddings[labels == 0]
    emb1 = embeddings[labels == 1]
    within0 = float(np.mean(cdist(emb0, emb0, metric="cosine")[np.triu_indices(len(emb0), k=1)]))
    within1 = float(np.mean(cdist(emb1, emb1, metric="cosine")[np.triu_indices(len(emb1), k=1)]))
    between = float(np.mean(cdist(emb0, emb1, metric="cosine")))
    distance_ratio = ((within0 + within1) / 2) / between

    final_info = {
        "final_train_loss": train_log[-1]["train_loss"],
        "final_val_accuracy": val_log[-1]["accuracy"],
        "final_val_auroc": val_log[-1]["auroc"],
        "final_val_f1": val_log[-1]["f1"],
        "best_val_auroc": float(max(v["auroc"] for v in val_log)),
        "best_val_accuracy": float(max(v["accuracy"] for v in val_log)),
        "kmeans_ari": ari,
        "silhouette_score": sil,
        "within_class_distance": (within0 + within1) / 2,
        "between_class_distance": between,
        "distance_ratio": distance_ratio,
        "n_train": int(n_train),
        "n_test": int(len(sequences) - n_train),
        "n_class0": int((labels == 0).sum()),
        "n_class1": int((labels == 1).sum()),
        "input_dim": int(input_dim),
        "model": model_name,
        "class0": label_names[0],
        "class1": label_names[1],
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
