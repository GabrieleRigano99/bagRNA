# bagRNA

Nextflow DSL2 pipeline for end-to-end eukaryotic genome annotation. Integrates RNA-seq evidence, ab initio gene predictors, structural annotation, and functional annotation into one reproducible workflow. All steps run in Docker containers.

![bagRNA pipeline overview](https://raw.githubusercontent.com/GabrieleRigano99/bagRNA/main/bagRNA_pipeline.png)

## What it does

- **Structural annotation** — genome + RNA-seq evidence → gene models (~25 steps: alignment, assembly, Mikado integration, cleanup)
- **Functional annotation** — gene models → InterProScan6, eggNOG-mapper, KEGG (DeepKOALA), Rfam, SignalP, dbCAN, MEROPS, PHI-base, AntiSMASH, plus an HTML/PDF report

See [Usage](usage.md) for every way to run it, or [Pipeline steps](structural-annotation.md) for what each stage does.

## Requirements

- [Nextflow](https://www.nextflow.io/) `>=25.04.6,<26`
- Java 21
- Docker
- (optional) NVIDIA GPU + [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html) for `--use_gpu`

```bash
export NXF_VER=25.04.6
export JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64
```

These two variables must be set before every `nextflow run`.
