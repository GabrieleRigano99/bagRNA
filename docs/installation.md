# Installation

## 1. Install Nextflow and Java 21

```bash
export NXF_VER=25.04.6
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

Nextflow needs Java 21. On Debian/Ubuntu:

```bash
sudo apt install openjdk-21-jre-headless
export JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64
```

Set both `NXF_VER` and `JAVA_HOME` in every shell before running the pipeline (add them to `~/.bashrc` to persist).

## 2. Install Docker

All processes run in Docker containers — no per-tool conda/module setup needed.

```bash
curl -fsSL https://get.docker.com | sh
sudo usermod -aG docker $USER   # log out/in after this
```

## 3. (Optional) GPU support

Required only for `--use_gpu` (Helixer, ANNEVO, TransDecoder2, InterProScan6 TMbed/SignalP6).

- NVIDIA driver installed on host
- [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html)

Verify with:

```bash
docker run --rm --gpus all nvidia/cuda:12.0.0-base-ubuntu22.04 nvidia-smi
```

## 4. Clone the pipeline

```bash
git clone https://github.com/GabrieleRigano99/bagRNA.git
cd bagRNA
```

## 5. Download reference databases (one-time)

```bash
nextflow run main.nf -entry SETUP --db_dir /path/to/databases
```

Downloads EggNOG, KEGG, RFAM, and other databases used by functional annotation. Re-run only if you need to refresh or relocate them.

Next: [Usage](usage.md).
