# ComRepka

**Comparative analysis of repetitive elements using conserved genomic flanks.**

ComRepka is a Nextflow pipeline for comparing repeat-associated genomic positions across assembled genomes. It identifies repeats in a source genome, builds flanking-sequence queries, aligns them to a target genome, and classifies the resulting correspondence. The scientific processing steps are implemented in Python; Nextflow manages dependencies, parallel execution, and partial runs.

The pipeline supports any assembled genomes for which suitable sequences and repeat annotations are available. It is not restricted to primates or to a particular repeat family.

> **Reproducibility:** The repository includes a 20 kb reference example with recorded outputs and checksums. Tool versions and environment specifications are documented in the [technical reference](docs/technical_reference.md#validation-and-reproducibility).

[Русская документация](README.ru.md) · [Detailed technical reference](docs/technical_reference.md)

## How it works

```text
Source genome (FASTA)
       │
       ▼
  RepeatMasker ── or ── existing repeat annotations (.out / BED6 / custom5)
       │
       ▼
Python: extract left/right genomic flanks and record source coordinates
       │
       ▼
minimap2: align flanks to a target genome
       │
       ▼
Python: classify alignment pairs → statuses.tsv
       │
       ├── Optional HAL-based validation
       └── Optional Jupyter analysis
```

One source genome can be compared against multiple targets. Every comparison is **directed** (`source → target`), and each target genome can reuse its minimap2 index across compatible comparisons. Individual stages can be entered using precomputed files.

## Requirements

- Linux x86_64 (the documented platform)
- Nextflow **25.10.2** or newer, plus Java (the recorded launcher used OpenJDK 17)
- Conda, with access to Bioconda and conda-forge
- RepeatMasker and a compatible repeat library **only** when starting with unannotated genomes

Nextflow can create the process environments defined in [`envs/`](envs/). The historical environment snapshots for Linux x86_64 are in [`envs/locked/linux-64/`](envs/locked/linux-64/); the main processing environment is specified in [`envs/core.yml`](envs/core.yml).

To recreate the recorded Nextflow launcher environment:

```bash
conda create --prefix "$PWD/.conda/nextflow" --file envs/locked/linux-64/nextflow.explicit.txt
conda activate "$PWD/.conda/nextflow"
```

If Nextflow is already available, use that installation instead. For reference-library installation and HPC-specific environments, see the [technical reference](docs/technical_reference.md#repeatmasker-and-reference-data).

## Quick start — no genome-scale computation

This example uses the **provided PAF alignments** and runs only the Python classifier. It does not run RepeatMasker or minimap2:

```bash
nextflow run main.nf -profile local \
  --genomes examples/genomes.tsv \
  --pairs examples/pairs.tsv \
  --from paf --until statuses \
  --outdir results/example-paf
```

The expected output path is:

```text
results/example-paf/pairs/genome_A__to__genome_B/statuses.tsv
```

For a reference of selected output fields, see [`examples/expected/classification_summary.tsv`](examples/expected/classification_summary.tsv). The example uses synthetic data and a hand-authored PAF; it is a workflow-interface example, not a biological validation dataset.

## Running your own comparison

Create two **tab-separated** manifests:

- `genomes.tsv` — genome identifiers, FASTA paths, and optionally precomputed repeat annotations or flanks. See [`examples/genomes.tsv`](examples/genomes.tsv).
- `pairs.tsv` — explicit source and target genome IDs, with optional precomputed PAF or status files. See [`examples/pairs.tsv`](examples/pairs.tsv).

Paths within a manifest are relative to that manifest's directory (or absolute). The genome IDs must be consistent across both manifests. Input genome assemblies use **FASTA**, not raw FASTQ reads.

A full run, including RepeatMasker when annotations are not already supplied, can be started with:

```bash
nextflow run main.nf -profile local \
  -params-file reference.yml \
  --genomes genomes.tsv --pairs pairs.tsv \
  --from genome --until statuses \
  --outdir results/my-study
```

`reference.yml` is your local configuration for the actual installed repeat database; it is deliberately not included in the repository. Preparing it is described in [RepeatMasker and reference data](docs/technical_reference.md#repeatmasker-and-reference-data). **Do not run the full command before configuring the reference library and computational resources.**

### Start from an intermediate stage

| `--from` | What you supply | What runs next |
|---|---|---|
| `genome` | Source FASTA, repeat annotation or RepeatMasker reference | Remaining selected stages |
| `annotation` | `.out`, BED6, or `custom5`; source FASTA | Flank extraction onward |
| `flanks` | Flank FASTA and flank TSV | minimap2 onward |
| `paf` | Pairwise PAF and source flank TSV | Classification onward |
| `statuses` | Pairwise status TSV and source flank TSV | Optional HAL validation |

Select the stopping point with `--until` (`annotation`, `flanks`, `paf`, `statuses`, or `validation`). Restart interrupted runs with `-resume` *only while the corresponding Nextflow cache and work directory are retained*. Starting from externally supplied results uses `--from`, not `-resume`.

[Complete parameter list and intermediate-format contracts](docs/technical_reference.md#full-and-partial-execution)

## Running on SLURM

A SLURM profile is included, with CPU, memory, time, partition and account parameters configurable for your cluster:

```bash
nextflow run main.nf -profile slurm \
  -params-file reference.yml \
  --genomes genomes.tsv --pairs pairs.tsv \
  --from genome --until statuses \
  --outdir results/my-study
```

Use a shared filesystem for the repository, data, Nextflow `work/`, and Conda cache. Configure the cluster-specific partition, account, filesystem and resource parameters before submitting jobs. See [SLURM and resources](docs/technical_reference.md#slurm-and-resources).

## Results

Final results are published under the chosen `--outdir`:

```text
results/my-study/
├── genomes/<genome_id>/
│   ├── repeats.bed
│   ├── flanks.fa
│   ├── flanks.tsv
│   └── targets.tsv
├── pairs/<source>__to__<target>/
│   ├── alignments.paf
│   └── statuses.tsv
└── provenance/
    └── ...
```

Only files produced by the selected run stages will be present. `statuses.tsv` contains the original **54-column classification output** with status calls, target positions, alignment evidence, and confidence-related fields. A high-confidence call is not independent proof of orthology. See [Outputs and interpretation](docs/technical_reference.md#outputs-and-interpretation) for field definitions and scientific caveats.

## Optional downstream analysis

The [`notebooks/`](notebooks/) directory contains three standalone notebooks:

- `pairwise_statistics.ipynb` — per-comparison status counts and optional repeat/HAL matching.
- `prepare_comparison_matrix.ipynb` — HDF5 representation of directed comparisons.
- `multigenome_statistics.ipynb` — numerical summaries of the comparison matrix.

They are **not** executed by Nextflow. Use [`envs/analysis.yml`](envs/analysis.yml), set manifest/result paths in the notebooks, and run them manually. Configure input and output paths in the notebooks for your dataset.

## Examples, documentation, and reproducibility

- [`examples/data/`](examples/data/) — small synthetic genomic sequences, annotations, flanks and PAF.
- [`examples/full/`](examples/full/) — an archived 20 kb technical positive control, with expected outputs and SHA-256 provenance. The expected outputs correspond to earlier tool versions.
- [Technical reference](docs/technical_reference.md) — installation details, RepeatMasker database setup, parameters, formats, full output schema, resource settings and validation limits.
- [Русская версия документации](README.ru.md) — подробное руководство на русском языке.

For reproducible scientific comparisons, record assembly accessions and sequence checksums, RepeatMasker database versions/partitions, software versions, all parameters, and the exact directed genome pairs. Large genomes and reference databases are not part of this repository.

## License and citation

**To be completed before public release:** the authors must choose a license and provide `LICENSE`, `CITATION.cff`, and release metadata. No license or DOI is claimed by this draft.
