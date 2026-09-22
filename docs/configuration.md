# Configuration

## Process-level overrides go in `nextflow.config`

Never edit `.nf` files after a run has started or if you intend to `-resume` — Nextflow fingerprints process definitions, and any change invalidates the cache for all downstream processes. Use `withName:` blocks instead:

```groovy
process {
    withName: 'STRUCTURAL_ANNOTATION:SOME_PROCESS' {
        containerOptions = '-v /extra/mount:/extra/mount'
        cache = 'deep'   // for directory inputs staged with stageAs
    }
}
```

## Adding a new optional parameter

1. Add to `params {}` in `nextflow.config` with a `null` default.
2. Resolve in the workflow with a null-check, falling back to `Channel.value(file("${projectDir}/assets/NO_FILE"))`.

## Directory input caching

Nextflow fingerprints directory inputs by path, not content. Use `stageAs: 'fixed_name'` + `cache 'deep'` to force content-based hashing when the input is a directory.

## Docker symlinks

Nextflow stages inputs as symlinks. If the real file lives outside the working directory (e.g. a database on a separate mount), add `-v /path/to/mount:/path/to/mount` to `containerOptions` so Docker can follow the symlink chain. A few processes (`RUN_DBCAN`, `TRANSDECODER2_PFAM_SCAN`, `EGGNOG`) already have such a mount configured for this machine's layout — edit the `withName:` block if your databases live elsewhere.

## Resource sizing

`process_high` / `process_gpu` / `process_long` labels default to (physical RAM − 4 GB). Override with `--max_memory` if your Docker daemon has a lower memory limit than the host.

## Execution profiles

| Profile | Executor |
|---|---|
| `standard` (default) | local |
| `slurm` | SLURM cluster |
| `sge` | SGE cluster |

```bash
nextflow run main.nf [params] -profile slurm
```

## GPU labels

`withLabel: 'use_gpu'` adds `--gpus all` to `containerOptions`. Requires the NVIDIA Container Toolkit on the host.

## Pull timeout

`pullTimeout = '1 hour'` is set globally to handle large image pulls (IPS6 member-database images in particular).
