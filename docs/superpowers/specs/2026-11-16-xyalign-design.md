# XYalign for Terra

Build a WDL 1.0 workflow in SexChromosome. Analyze indexed BAM files in a sample scatter. Run XYalign chromosome depth and mapping-quality statistics by default. Offer detailed chromosome characterization and read remapping as separate options. The user selected chromosome analysis with optional remapping.

## Design

Use typed sample inputs: sample ID, BAM, BAM index, and optional Y status. Preserve File types until each task renders its command. Require a reference FASTA and optional FASTA index. Accept chromosome names, analysis settings, resource settings, and an optional reference mask. Do not use a JSON argument wrapper.

Use three tasks: analyze one BAM, prepare masked references once, and remap one BAM. Prepare and remap only when enabled. Remapping requires an explicit Y status for each sample. Analysis must finish before remapping that sample. Return statistics, report archives, logs, and optional remapped BAMs with indexes. Preserve sample order in outputs.

Use a Python wrapper compatible with Python 2.7 and Python 3. Validate all inputs before XYalign runs. Use safe task-local names because upstream commands use shell strings. Reject unresolved cloud paths. Validate BAM sort order, indexes, contig names, lengths, and order, and reference masks. Copy or link readable files into writable task directories. Supply only named CLI arguments. Create all argument lists inside tasks.

Use XYalign 1.1.6 at source commit a68f484cd256dec37abb1166bb3a8a0fb994351f. Use a micromamba base and pinned conda-forge/bioconda packages. Keep the Python 2 environment inside the container. The wrapper must not change upstream scientific calculations. Document the upstream permutation p-value defect. Disable characterization by default. Disable variant calling unless requested.

## Approaches

A single upstream full-pipeline command couples inference to remapping. Separate tasks permit explicit remapping decisions and shared reference preparation. A new implementation of the scientific method would exceed the requested tool wrapper. Use separate tasks.

## Checks

Run local wrapper unit tests, WDL syntax checks, and an AST check that rejects file-writing functions at workflow scope. Test commands after a fake cloud-to-local mapping. Use synthetic BAM/reference data in GitHub Actions for real analysis, characterization, reference preparation, and both remapping branches. Build no local Docker image. Trigger image rebuilds automatically only when the contained wrapper changes. Use manual dispatch for environment changes. Publish images only on main after smoke tests pass. Open a draft PR. Do not merge it or submit cloud jobs.

The complete workflow has not been tested on Terra. GitHub Actions smoke tests do not establish Terra compatibility.
