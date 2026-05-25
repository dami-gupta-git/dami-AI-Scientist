# RunPod Connection Reference

## SSH Key

Working key: `~/.ssh/id_runpod_2` (registered as `runpod_2` in RunPod account settings).

Public key:
```
ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIMuXdSvZt602kwA42h3d78lyGwRgK35z4TA26mG8qlhw runpod_2
```

## Connect

```bash
ssh -o StrictHostKeyChecking=no -i ~/.ssh/id_runpod_2 root@<POD_IP> -p <PORT>
```

## Current pod (May 25 2026)

```bash
ssh -o StrictHostKeyChecking=no -i ~/.ssh/id_runpod_2 root@216.81.245.143 -p 11019
```

A100 SXM4 80GB, 128 cores.

## New pod setup

RunPod injects `runpod_2` automatically on pod start. If it doesn't appear, run in the **pod web terminal**:

```bash
echo "ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIMuXdSvZt602kwA42h3d78lyGwRgK35z4TA26mG8qlhw runpod_2" > ~/.ssh/authorized_keys
chmod 600 ~/.ssh/authorized_keys
```

Then install deps:

```bash
apt-get update -y && apt-get install -y tmux git
pip install fair-esm scikit-learn scipy
git clone https://github.com/dami-gupta-git/dami-AI-Scientist.git /workspace/dami-AI-Scientist
cd /workspace/dami-AI-Scientist && git checkout esm2-mechanism
```


Git token for dami-gupta-git
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