# XYalign for Terra

Run XYalign chromosome analysis on indexed BAM files. Enable read remapping when you have selected the Y status for each sample.

<!-- workflow-badges:start -->
[![Docker Image CI](https://github.com/AoU-Multiomics-Analysis/SexChromosome/actions/workflows/docker-image.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/SexChromosome/actions/workflows/docker-image.yml)
[![R lint](https://github.com/AoU-Multiomics-Analysis/SexChromosome/actions/workflows/r-lint.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/SexChromosome/actions/workflows/r-lint.yml)
[![Update README workflow badges](https://github.com/AoU-Multiomics-Analysis/SexChromosome/actions/workflows/update-readme-badges.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/SexChromosome/actions/workflows/update-readme-badges.yml)
[![WDL validation](https://github.com/AoU-Multiomics-Analysis/SexChromosome/actions/workflows/wdl-validation.yml/badge.svg)](https://github.com/AoU-Multiomics-Analysis/SexChromosome/actions/workflows/wdl-validation.yml)
<!-- workflow-badges:end -->

This repository uses WDL 1.0. It targets Terra's managed Cromwell environment. XYalign source and scientific calculations are unchanged.

## Workflow

Import `workflows/xyalign.wdl` through Dockstore or the Terra workflow repository.

The workflow processes samples in parallel:

1. Check the localized BAM and reference files.
2. Run XYalign `CHROM_STATS` for chromosome depth and mapping quality.
3. Run `CHARACTERIZE_SEX_CHROMS` when `characterize` is true.
4. Prepare shared masked references when `run_remapping` is true.
5. Remap X/Y reads when `run_remapping` is true. Use each sample's explicit Y status.

Reference preparation can run in parallel with analysis. Remapping waits for that sample's analysis and the prepared references.

The default workflow performs chromosome statistics only. It does not select Y status automatically. Characterization, variant calling, and remapping are separate options.

## Container

GitHub Actions builds a Linux AMD64 image. The image uses `mambaorg/micromamba:2.10.0` and pinned conda-forge and bioconda packages. It installs unchanged XYalign 1.1.6 from commit `a68f484cd256dec37abb1166bb3a8a0fb994351f`. The build verifies the source archive checksum.

XYalign requires Python 2.7. The container supplies this runtime and compatible package versions. Your local computer does not need Python 2.

The Docker workflow builds and tests each image before publication. Pull requests build and test without publication. A successful run on `main` publishes:

- `ghcr.io/aou-multiomics-analysis/sexchromosome:main`
- `ghcr.io/aou-multiomics-analysis/sexchromosome:sha-<commit>`

Set the GHCR package visibility to public before Terra uses it. The default image tag is `main`. Set `XYalign.docker_image` to the published image digest for reproducible runs.

Automatic image builds run only when `scripts/xyalign_runner.py` changes. Use **Run workflow** after changing the Dockerfile or environment. WDL-only and documentation-only changes do not rebuild the image. No local Docker build is required.

## Inputs

Start with `inputs/xyalign.analysis.example.json`. Replace the example cloud paths with accessible files.

| Input | Type | Use |
| --- | --- | --- |
| `samples` | Array of `XYalignSample` | One entry per BAM. The array must not be empty. |
| `samples[].sample_id` | String | Output label. Use letters, numbers, underscores, dots, colons, or hyphens. Start with a letter, number, or underscore. |
| `samples[].bam` | File | Coordinate-sorted BAM. SAM and CRAM are not supported. |
| `samples[].bam_index` | File | BAM index. Declare it separately so Terra can localize it. |
| `samples[].y_status` | Optional String | `present` or `absent`. Required when remapping is enabled. |
| `reference_fasta` | File | Uncompressed reference FASTA used for the BAM alignment. |
| `reference_fasta_index` | Optional File | Matching FASTA index. The task creates an index if absent. |
| `reference_mask` | Optional File | BED intervals to mask in both remapping references. |
| `chromosomes` | Array of String | Default: `chr19`, `chrX`, `chrY`. Include X, Y, and at least one autosome. |
| `x_chromosome`, `y_chromosome` | String | Exact contig names. Defaults: `chrX`, `chrY`. |
| `characterize` | Boolean | Default: false. Add window metrics, plots, and chromosome comparisons. |
| `call_variants` | Boolean | Default: false. Enable Platypus calls during characterization. Requires `characterize: true`. |
| `run_remapping` | Boolean | Default: false. Remap reads with the supplied Y status. |
| `single_end` | Boolean | Default: false. Set true for single-end reads. |
| `window_size` | Int | Default: 5000 bases. Must be positive. |
| `num_permutations`, `num_bootstraps` | Int | Default: 10000 each. Must be positive. |
| `docker_image` | String | Container tag or digest accessible to Terra. |

BAM and FASTA contig names, lengths, and order must match. Chromosome names use the same permitted characters as sample IDs. The wrapper rejects unreadable inputs and unresolved cloud URIs before analysis.

The workflow preserves File types inside sample structs. Terra localizes the declared files. The tasks then pass local paths through named CLI arguments. A path listed inside a manifest is not an input dependency.

Use `inputs/xyalign.remap.example.json` for remapping. The Y status selects the reference. It is not a statement about gender or a complete karyotype. Supply a reference mask that matches your assembly and analysis design. No human PAR coordinates are assumed.

## Resources

| Step | CPUs | Memory | Working disk |
| --- | ---: | ---: | ---: |
| Analysis | 2 | 16 GB | 100 GB |
| Reference preparation | 4 | 16 GB | 200 GB |
| Remapping | 8 | 32 GB | 300 GB |

Set `analysis_cpu`, `analysis_memory_gb`, and `analysis_disk_gb` to change analysis resources. Use the corresponding `reference_*` and `remap_*` inputs for the other tasks. Increase disk space for large BAMs. Remapping can create FASTQ files and temporary BAMs. Each remapping task also localizes and extracts the reference archive.

`java_memory` defaults to `8g`. Keep this allocation below the remapping memory limit. `preemptible_attempts` defaults to 0.

## Outputs

Output arrays follow the input sample order. `sample_ids` gives the labels for that order.

| Output | Contents |
| --- | --- |
| `chromosome_depth` | Per-sample chromosome depth tables. |
| `chromosome_mapq` | Per-sample chromosome mapping-quality tables. |
| `analysis_reports` | Per-sample `.tar.gz` archives with upstream results, logs, and optional characterization files. |
| `analysis_logs` | Wrapper command logs. |
| `prepared_references` | Optional `.tar` archive with both masked FASTAs and indexes. |
| `reference_log` | Optional preparation command log. |
| `remapped_bams` | Optional BAM for each sample. XYalign merges remapped sex-chromosome reads with the other alignments. |
| `remapped_bam_indexes` | Optional BAM indexes. |
| `remapping_reports`, `remapping_logs` | Optional upstream logs and wrapper logs. |

Remapping arrays contain null entries when remapping is disabled. The statistics retain upstream table contents. Their BAM column is named `input.bam`, the fixed task-local alias. Use `sample_ids` to identify each table.

XYalign's characterization file with a `.csv` extension is tab-delimited. Archives preserve the upstream file names. Remapping does not add a duplicate-marking step. It remaps reads selected from X/Y alignments, rather than all original reads.

## Known upstream limits

The pipeline does not patch XYalign.

In `xyalign/ploidy.py`, the permutation p-value calculation uses `len(np.where(perms > observed)) / num_perms`. For a one-dimensional array, the numerator is always 1. With 10000 permutations, the reported value is 0.0001. Do not interpret these values as valid permutation p-values. This limit applies to optional characterization. See the [pinned upstream source](https://github.com/SexChrLab/XYalign/blob/a68f484cd256dec37abb1166bb3a8a0fb994351f/xyalign/ploidy.py#L63-L72).

The wrapper keeps the upstream default duplicate handling. It does not expose `--ignore_duplicates`, whose implementation appears opposite to its help text. It keeps the chromosome tests enabled because the upstream test-disable flags can leave result variables undefined.

## Validation

The complete workflow has **not been tested on Terra**. No cloud jobs were submitted. WDL syntax checks and container smoke tests do not establish Terra compatibility.

The repository includes:

- A `miniwdl check` syntax check.
- An AST regression check for workflow-scope file-writing functions, including nested workflow expressions.
- Unit tests for unresolved cloud paths and safe task-local aliases.
- Command-rendering tests for all tasks after simulated cloud-to-local localization.
- GitHub Actions smoke tests with synthetic data for analysis, characterization, masked references, and both Y-status remapping branches.

Run local checks with Python 3.11 and miniwdl installed:

```bash
miniwdl check workflows/xyalign.wdl
python tools/check_terra_wdl.py workflows/xyalign.wdl
python -m unittest discover -s tests -p 'test_*.py' -v
```

Use GitHub Actions for the real container smoke test. Its artifacts contain the synthetic inputs and outputs. Approval is required before any Terra test submission.

## Upstream

[XYalign source and documentation](https://github.com/SexChrLab/XYalign)

Webster TH et al. 2019. Identifying, understanding, and correcting technical artifacts on the sex chromosomes in next-generation sequencing data. *GigaScience* 8(7): giz074. [Publication](https://doi.org/10.1093/gigascience/giz074).
