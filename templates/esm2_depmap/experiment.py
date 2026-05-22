"""
Does ESM-2 embedding distance correlate with CRISPR essentiality correlation in DepMap?

Pipeline:
1. Download DepMap 26Q1 CRISPRGeneEffect.csv (Chronos-corrected scores)
2. Sample N_GENES random genes that appear in both DepMap and Ensembl
3. Fetch canonical protein sequences from Ensembl
4. Extract ESM-2 650M embeddings (frozen)
5. Compute pairwise cosine distances in embedding space
6. Compute pairwise Pearson correlations of essentiality profiles across cell lines
7. Mantel test: is there a correlation between embedding distance and essentiality correlation?
8. Baseline: same test with random shuffled embeddings
"""

import argparse
import json
import os
import time
import urllib.request

import numpy as np

N_GENES = 2000
ESM2_MODEL = "esm2_t33_650M_UR50D"
ENSEMBL_REST = "https://rest.ensembl.org"
DEPMAP_URL = "https://depmap.org/portal/api/download/files?release=depmap_public_26q1"
DEPMAP_FILENAME = "CRISPRGeneEffect.csv"


# ---------------------------------------------------------------------------
# DepMap
# ---------------------------------------------------------------------------

def get_depmap_download_url(filename=DEPMAP_FILENAME):
    import csv, io
    with urllib.request.urlopen(DEPMAP_URL, timeout=30) as resp:
        text = resp.read().decode()
    reader = csv.reader(io.StringIO(text))
    for row in reader:
        if len(row) >= 4 and row[2].strip() == filename and "26Q1" in row[0]:
            return row[3].strip()
    raise ValueError(f"{filename} not found in DepMap 26Q1 release")


def download_depmap(cache_dir):
    os.makedirs(cache_dir, exist_ok=True)
    dest = os.path.join(cache_dir, DEPMAP_FILENAME)
    if os.path.exists(dest):
        print(f"Using cached DepMap file: {dest}")
        return dest
    print("Downloading DepMap CRISPRGeneEffect.csv (~500MB)...")
    url = get_depmap_download_url()
    urllib.request.urlretrieve(url, dest)
    print(f"Downloaded to {dest}")
    return dest


def load_depmap(cache_dir):
    """Load DepMap gene effect matrix. Returns (gene_symbols, cell_lines, matrix)."""
    import pandas as pd
    path = download_depmap(cache_dir)
    print("Loading DepMap matrix...")
    df = pd.read_csv(path, index_col=0)
    # Columns are "GENE (ENTREZ)" format — extract gene symbol
    df.columns = [c.split(" ")[0] for c in df.columns]
    # Rows are cell line IDs, columns are genes
    print(f"DepMap shape: {df.shape} ({df.shape[0]} cell lines, {df.shape[1]} genes)")
    return df


# ---------------------------------------------------------------------------
# Ensembl protein sequences
# ---------------------------------------------------------------------------

def fetch_protein_sequence(symbol, retries=3, delay=1.0):
    """Fetch canonical protein sequence for a gene symbol."""
    # Step 1: symbol -> ensembl ID
    url = f"{ENSEMBL_REST}/lookup/symbol/homo_sapiens/{symbol}?content-type=application/json"
    info = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                info = json.loads(resp.read())
            break
        except Exception:
            if attempt < retries - 1:
                time.sleep(delay)
    if info is None:
        return None

    ensembl_id = info["id"]
    time.sleep(0.3)

    # Step 2: ensembl ID -> protein sequence
    url2 = (f"{ENSEMBL_REST}/sequence/id/{ensembl_id}"
            f"?content-type=application/json&type=protein&multiple_sequences=1")
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url2, timeout=30) as resp:
                data = json.loads(resp.read())
            if isinstance(data, list):
                return max(data, key=lambda x: len(x["seq"]))["seq"].upper()
            return data["seq"].upper()
        except Exception:
            if attempt < retries - 1:
                time.sleep(delay)
    return None


# ---------------------------------------------------------------------------
# ESM-2 embeddings
# ---------------------------------------------------------------------------

def get_esm2_embeddings(sequences, model_name=ESM2_MODEL, device="cuda", batch_size=8):
    import esm
    import torch
    model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
    model = model.to(device).eval()
    batch_converter = alphabet.get_batch_converter()

    all_embeddings = []
    for i in range(0, len(sequences), batch_size):
        batch = sequences[i: i + batch_size]
        batch_trunc = [s[:1022] for s in batch]
        data = [(f"seq{j}", s) for j, s in enumerate(batch_trunc)]
        _, _, tokens = batch_converter(data)
        tokens = tokens.to(device)
        with torch.no_grad():
            results = model(tokens, repr_layers=[model.num_layers])
        reps = results["representations"][model.num_layers]
        for k, seq in enumerate(batch_trunc):
            emb = reps[k, 1: len(seq) + 1].mean(0).cpu().float().numpy()
            all_embeddings.append(emb)
        if (i // batch_size) % 10 == 0:
            print(f"  Embedded {min(i + batch_size, len(sequences))}/{len(sequences)}")
    return np.stack(all_embeddings)


# ---------------------------------------------------------------------------
# Mantel test
# ---------------------------------------------------------------------------

def mantel_test(dist_matrix_1, dist_matrix_2, n_permutations=1000, seed=0):
    """
    Mantel test: correlation between two distance matrices.
    Returns (r, p_value, null_distribution).
    """
    from scipy.stats import pearsonr, spearmanr
    rng = np.random.RandomState(seed)
    n = dist_matrix_1.shape[0]
    idx = np.triu_indices(n, k=1)
    v1 = dist_matrix_1[idx]
    v2 = dist_matrix_2[idx]

    r_obs, _ = spearmanr(v1, v2)

    null = []
    perm_idx = np.arange(n)
    for _ in range(n_permutations):
        rng.shuffle(perm_idx)
        v2_perm = dist_matrix_2[np.ix_(perm_idx, perm_idx)][idx]
        r_perm, _ = spearmanr(v1, v2_perm)
        null.append(r_perm)

    null = np.array(null)
    p_value = float(np.mean(np.abs(null) >= np.abs(r_obs)))
    return float(r_obs), p_value, null.tolist()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run(out_dir, seed=0, n_genes=N_GENES, model_name=ESM2_MODEL, n_permutations=1000):
    import torch
    os.makedirs(out_dir, exist_ok=True)
    np.random.seed(seed)
    device = str(torch.device("cuda" if torch.cuda.is_available() else "cpu"))

    cache_dir = os.path.join(out_dir, "data")
    os.makedirs(cache_dir, exist_ok=True)

    # --- Load DepMap ---
    depmap_df = load_depmap(cache_dir)
    depmap_genes = list(depmap_df.columns)
    print(f"DepMap genes available: {len(depmap_genes)}")

    # --- Sample genes ---
    dataset_cache = os.path.join(cache_dir, f"dataset_{n_genes}.json")
    if os.path.exists(dataset_cache):
        print("Loading cached gene dataset...")
        with open(dataset_cache) as f:
            dataset = json.load(f)
    else:
        rng = np.random.RandomState(seed)
        sampled = list(rng.choice(depmap_genes, size=min(n_genes * 2, len(depmap_genes)), replace=False))
        dataset = []
        print(f"Fetching protein sequences for up to {len(sampled)} genes...")
        for sym in sampled:
            if len(dataset) >= n_genes:
                break
            seq = fetch_protein_sequence(sym)
            if seq and len(seq) >= 10:
                dataset.append({"symbol": sym, "seq": seq})
                if len(dataset) % 100 == 0:
                    print(f"  Fetched {len(dataset)} sequences...")
        with open(dataset_cache, "w") as f:
            json.dump(dataset, f)
        print(f"Dataset: {len(dataset)} genes with protein sequences")

    symbols = [d["symbol"] for d in dataset]
    sequences = [d["seq"] for d in dataset]
    print(f"Using {len(symbols)} genes")

    # --- Extract embeddings ---
    emb_cache = os.path.join(cache_dir, f"embeddings_{model_name}_{len(symbols)}.npy")
    if os.path.exists(emb_cache):
        print("Loading cached embeddings...")
        embeddings = np.load(emb_cache)
    else:
        print(f"Extracting ESM-2 embeddings ({model_name})...")
        embeddings = get_esm2_embeddings(sequences, model_name, device)
        np.save(emb_cache, embeddings)
    print(f"Embeddings shape: {embeddings.shape}")

    # --- Compute pairwise ESM-2 cosine distances ---
    print("Computing pairwise ESM-2 cosine distances...")
    from scipy.spatial.distance import pdist, squareform
    emb_norm = embeddings / np.linalg.norm(embeddings, axis=1, keepdims=True)
    esm2_dist = squareform(pdist(emb_norm, metric="cosine"))

    # --- Compute pairwise essentiality correlations ---
    print("Computing pairwise essentiality correlations...")
    ess_matrix = depmap_df[symbols].values  # shape: (cell_lines, n_genes)
    # Remove cell lines with NaN for any of our genes
    valid_rows = ~np.any(np.isnan(ess_matrix), axis=1)
    ess_matrix = ess_matrix[valid_rows]
    print(f"  Using {valid_rows.sum()} cell lines (after NaN removal)")

    # Pearson correlation between essentiality profiles
    from numpy.linalg import norm
    ess_centered = ess_matrix - ess_matrix.mean(axis=0, keepdims=True)
    ess_norm = ess_centered / (norm(ess_centered, axis=0, keepdims=True) + 1e-8)
    ess_corr = (ess_norm.T @ ess_norm) / ess_matrix.shape[0]
    # Convert correlation to distance: dist = 1 - corr
    ess_dist = 1 - ess_corr
    np.fill_diagonal(ess_dist, 0)

    # --- Mantel test ---
    print(f"Running Mantel test ({n_permutations} permutations)...")
    r_obs, p_value, null_dist = mantel_test(esm2_dist, ess_dist, n_permutations=n_permutations, seed=seed)
    print(f"Mantel r = {r_obs:.4f}, p = {p_value:.4f}")

    # --- Summary statistics ---
    idx = np.triu_indices(len(symbols), k=1)
    esm2_dists_flat = esm2_dist[idx]
    ess_dists_flat = ess_dist[idx]

    # Top 1% closest ESM-2 pairs — are they more correlated?
    threshold = np.percentile(esm2_dists_flat, 1)
    top_pairs_mask = esm2_dists_flat <= threshold
    mean_ess_corr_top = float(np.mean(1 - ess_dists_flat[top_pairs_mask]))
    mean_ess_corr_random = float(np.mean(1 - ess_dists_flat))

    final_info = {
        "mantel_r": r_obs,
        "mantel_p_value": p_value,
        "n_genes": len(symbols),
        "n_cell_lines": int(valid_rows.sum()),
        "embedding_dim": int(embeddings.shape[1]),
        "mean_esm2_dist": float(np.mean(esm2_dists_flat)),
        "mean_ess_correlation": float(np.mean(1 - ess_dists_flat)),
        "mean_ess_corr_top1pct_esm2_pairs": mean_ess_corr_top,
        "mean_ess_corr_random_pairs": mean_ess_corr_random,
        "enrichment_top1pct": mean_ess_corr_top - mean_ess_corr_random,
        "null_mean_r": float(np.mean(null_dist)),
        "null_std_r": float(np.std(null_dist)),
    }
    print(final_info)

    with open(os.path.join(out_dir, f"final_info_seed{seed}.json"), "w") as f:
        json.dump(final_info, f, indent=2)

    # Save distance matrices for plotting
    np.save(os.path.join(out_dir, "esm2_dist.npy"), esm2_dist)
    np.save(os.path.join(out_dir, "ess_dist.npy"), ess_dist)
    np.save(os.path.join(out_dir, "null_dist.npy"), np.array(null_dist))
    with open(os.path.join(out_dir, "symbols.json"), "w") as f:
        json.dump(symbols, f)

    return final_info


parser = argparse.ArgumentParser()
parser.add_argument("--out_dir", type=str, default="run_0")
parser.add_argument("--n_genes", type=int, default=N_GENES)
parser.add_argument("--n_permutations", type=int, default=1000)
args = parser.parse_args()

if __name__ == "__main__":
    import numpy as np

    final_info = run(args.out_dir, seed=0, n_genes=args.n_genes, n_permutations=args.n_permutations)

    final_infos = {
        "means": {k: v for k, v in final_info.items() if isinstance(v, (int, float))},
        "stderrs": {k: 0.0 for k, v in final_info.items() if isinstance(v, (int, float))},
        "final_info_list": [final_info],
    }

    with open(os.path.join(args.out_dir, "final_info.json"), "w") as f:
        json.dump(final_infos, f, indent=2)

    import numpy as np
    all_results = {"seed0_final_info": final_info}
    with open(os.path.join(args.out_dir, "all_results.npy"), "wb") as f:
        np.save(f, all_results)
