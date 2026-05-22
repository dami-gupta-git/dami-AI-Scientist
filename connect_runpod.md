# RunPod SSH Setup

## Current Pod (May 22 2026)
- **IP:** 216.81.245.125
- **Port:** 10075
- **User:** root
- **GPU:** NVIDIA A100-SXM4-80GB, 80GB VRAM
- **Working directory:** /esm2val
- **Proxy:** h9apndcravlmsq-64410c90@ssh.runpod.io

## SSH Key

The key that works is `~/.ssh/id_runpod` (unencrypted ed25519).

Public key (must be in pod's `~/.ssh/authorized_keys`):
```
ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIKiI2zsJyXGzNwpNSNa7plSjAwxTKEIFAkhZLuBHj92j dgupta@Mac.hsd1.ma.comcast.net
```

**Important:** New RunPod pods often have a corrupt or wrong `authorized_keys`. Always do this first in the pod terminal:
```bash
echo "ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIKiI2zsJyXGzNwpNSNa7plSjAwxTKEIFAkhZLuBHj92j dgupta@Mac.hsd1.ma.comcast.net" > ~/.ssh/authorized_keys
chmod 600 ~/.ssh/authorized_keys
```

## How to Connect

### From your local terminal
```bash
ssh -i ~/.ssh/id_runpod root@216.81.245.125 -p 10075
```

### For Claude to run commands programmatically
```bash
# Step 1: write key to temp file (do this at start of each session)
cat > /tmp/runpod_key << 'EOF'
-----BEGIN OPENSSH PRIVATE KEY-----
b3BlbnNzaC1rZXktdjEAAAAABG5vbmUAAAAEbm9uZQAAAAAAAAABAAAAMwAAAAtzc2gtZW
QyNTUxOQAAACCoiNs7CclxszcKTUjWu6ZUowMMUyhCBQJIWS7gR4/dowAAAKj2aahM9mmo
TAAAAAtzc2gtZWQyNTUxOQAAACCoiNs7CclxszcKTUjWu6ZUowMMUyhCBQJIWS7gR4/dow
AAAECuj0NFcVrvrIYO+oGUdIR9AXX3X8VkRIf2CrtmkM6SfKiI2zsJyXGzNwpNSNa7plSj
AwxTKEIFAkhZLuBHj92jAAAAHmRndXB0YUBNYWMuaHNkMS5tYS5jb21jYXN0Lm5ldAECAw
QFBgc=
-----END OPENSSH PRIVATE KEY-----
EOF
chmod 600 /tmp/runpod_key

# Step 2: run commands
ssh -o StrictHostKeyChecking=no -i /tmp/runpod_key root@216.81.245.125 -p 10075 "your command here"
```

## New Pod Setup Checklist

When a new RunPod pod is created:

1. Open pod terminal in browser and fix authorized_keys (see above)
2. Install tmux: `apt-get update -q && apt-get install -y tmux`
3. Install Python deps: `pip install fair-esm scikit-learn scipy -q`
4. Clone repo: `git clone https://github.com/dami-gupta-git/dami-AI-Scientist.git /esm2val`
5. Checkout branch: `cd /esm2val && git checkout <branch>`
6. Upload any required cached data via scp

## Notes
- The proxy connection (`ssh.runpod.io`) also works but is sometimes flaky — prefer direct IP
- Always use `-o StrictHostKeyChecking=no` to avoid host key prompts
- The old pod (103.196.86.40:38349) is no longer active

## Repo Branches
| Branch | Purpose |
|---|---|
| `evo2-template` | Evo2 DNA embeddings — DNA repair vs TSG |
| `evo2-validation` | Evo2 DNA embeddings — glycolysis vs T-cell |
| `esm2-function` | ESM-2 protein embeddings — DNA repair vs TSG |
| `esm2-depmap` | ESM-2 vs DepMap Mantel test (2000 genes) |
| `esm2-validation` | ESM-2 protein embeddings — glycolysis vs T-cell |
