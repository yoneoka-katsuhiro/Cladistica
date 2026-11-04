# Cladistica manual (English)

[Overview](../README.md) · [日本語](MANUAL.ja.md)

## Installation

The native ZIP packages include the Cladistica CLI and Python dependencies. Extract the complete archive and keep the executable beside `_internal`. In the examples below, replace `cladistica` with `./cladistica` on macOS or `.\cladistica.exe` in Windows PowerShell. Use one line for each command in PowerShell; Bash's `\` line continuation is not a PowerShell continuation.

For source use, install Python 3.10+, create a virtual environment, activate it, and run `python -m pip install .` in the repository. The `setup_mac.sh` and `run_cladistica.sh` scripts remain available for existing macOS users. Installed CLI commands write default results under the current working directory, not inside the Python installation.

Install only the external programs needed for the stages you run:

| Stage | Program | Location/option |
| --- | --- | --- |
| Alignment | [MUSCLE](https://drive5.com/muscle5/) | `--muscle /path/to/muscle` |
| Maximum likelihood | [IQ-TREE](https://iqtree.github.io/) | `--iqtree /path/to/iqtree3` |
| Bayesian inference | [MrBayes](https://nbisweden.github.io/MrBayes/) 3.2+ | `--mrbayes /path/to/mb` |

On Windows, quote paths containing spaces: `--mrbayes "C:\Tools\MrBayes\mb.exe"`. Each executable must also have its own required libraries/runtime installed. Cladistica does not install or redistribute these programs. `--skip-ml` and `--skip-bi` omit the corresponding inference stages. `--threads` applies to IQ-TREE; it does not configure MPI for MrBayes. The current CLI launches one MrBayes executable directly.

For NCBI searches or accession downloads, provide `--email your.email@example.com` or set `NCBI_EMAIL` in a `.env` file in the working directory. An optional API key belongs in `NCBI_API_KEY`. Keep `.env` private. Local FASTA workflows do not contact NCBI.

## First checks

```bash
cladistica --version
cladistica self-check
cladistica --help
cladistica tree --input examples --output preview --dry-run
```

The last command creates NEXUS, partition and statistics files without running external programs. The four synthetic taxa in `examples/tiny.fasta` are a functional example, not a research dataset.

## Survey and curate GenBank accessions

```bash
cladistica survey --genus Hymenasplenium --outgroup "Asplenium setoi" --markers rbcL atpA atpB --email your.email@example.com --output survey_run
```

`accession_all.csv` records candidates and automatic recommendations. Inspect the source records, taxon names, specimen/voucher identity, extracted region, length and QC before choosing records. By default, `survey` retrieves all matching search results; `--retmax 100` bounds a small trial. `run` has a different default limit of 500 records per query name.

Create `accession_selected.csv` with one row per specimen, a `taxon` or `sample_id`, and one column per requested marker. A marker cell may contain at most one accession. Blank marker cells represent missing loci. An explicit stable `sample_id` is recommended:

```csv
sample_id,taxon,voucher_or_isolate,rbcL,matK
specimen_001,Example species,V001,ACCESSION_1,ACCESSION_2
```

The accessions above are placeholders; do not submit them. Without `sample_id`, Cladistica uses the taxon and an explicitly supplied voucher. This identifier remains the same across marker columns. Repeated or normalization-colliding identifiers are rejected. Use separate stable IDs for distinct specimens of the same species. An ID cannot contain `|` or a control character.

```bash
cladistica resume --accession-selected input/accession_selected.csv --accession-all survey_run/accession_all.csv --markers rbcL atpA atpB --email your.email@example.com --output selected_run
```

`--accession-all` is optional and retains the survey in the final result. `run --genus ...` performs automatic selection and all following stages; manual review before inference is preferable when specimen identity is uncertain.

## Supply your own sequences

For unaligned data, use one file per marker (`rbcL.fasta`, `matK.fasta`, etc.). In a FASTA header, text before the first `|` is the sample ID. Use the same ID across loci. Duplicate headers, empty sequences and malformed FASTA are errors.

```bash
cladistica resume --fasta-dir input/markers --markers rbcL matK --output my_run
```

Add `--add-fasta-dir input/additional` to an accession or marker workflow to add your sequences. A matching sample and marker in the additional directory **replaces** the primary sequence. Confirm this intended override before analysis. `--markers` must list the actual loci; the default marker set is not automatically inferred from arbitrary filenames.

For an already aligned and concatenated matrix:

```bash
cladistica resume --concatenated-fasta input/matrix.fasta --partition-file input/partitions.tsv --output matrix_run
```

The partition file may be CSV/TSV with `marker,start,end`, IQ-TREE text (`DNA, rbcL = 1-100`), or NEXUS charset lines (`charset rbcL = 1-100;`). Coordinates are 1-based and inclusive. Without a partition file, the concatenated matrix is one partition. All sequences must have equal length. Review partition boundaries yourself before inference.

## Methods and interpretation

Marker extraction uses GenBank annotations and curated aliases. CDS quality checks account for translated coding regions; noncoding loci are not subjected to coding-sequence stop-codon checks. Introns and spacers require identifiable boundaries. Whole plastomes without suitable annotations are not automatically annotated or searched by similarity.

MUSCLE aligns each marker separately. Concatenation joins the alignments by sample ID and fills an absent marker with `?`. It does not establish that loci share a genealogy. Review paralogy, orthology, recombination, alignment quality and missing data for your study.

ML uses IQ-TREE `-m MFP+MERGE` with partitions and **standard nonparametric bootstrap** (`-b`, default 1000). It does not use ultrafast bootstrap. IQ-TREE 2 or 3 is recommended. The chosen models and raw output are retained under `_work/trees/`.

BI uses GTR+Gamma (`nst=6 rates=gamma`), Dirichlet(1,1,1,1) base-frequency priors, and, with multiple partitions, unlinked state frequencies, exchangeabilities and gamma shapes plus variable partition rates. Other priors retain MrBayes defaults. The default is two independent runs, four chains per run, 1,000,000 generations, sampling/checkpointing every 1,000 generations, diagnostics every 5,000 generations and 25% relative burn-in. The `tree` command exposes `--nruns`, `--nchains` and `--burninfrac`; the higher-level workflows use their defaults. Burn-in is applied explicitly to both `sump` and `sumt`. BI model settings are separate from ML model selection.

A completed run is not necessarily converged. Compare independent traces and likelihood stationarity, inspect ESS and MrBayes convergence statistics, and assess sensitivity to burn-in and run length. A 2,000-generation demonstration is far too short for substantive inference. Including outgroup accessions does not automatically root trees. Final interpretation and rooting are the investigator's responsibility.

## BI recovery and extension

Keep the complete `BI_resume/` folder. Its original basenames are part of the MrBayes restart mechanism:

| File | Purpose |
| --- | --- |
| `mrbayes_analysis.nex` | Original data, model, priors and MCMC controls |
| `mrbayes_analysis.nex.ckp` and `.ckp~` | Current and backup chain checkpoints |
| `mrbayes_analysis.nex.runN.p` | Parameter samples for each run |
| `mrbayes_analysis.nex.runN.t` | Tree samples for each run |
| `mrbayes_analysis.nex.mcmc` | MCMC diagnostics |
| `resume_manifest.json` | Versions, input/output hashes, state and continuation history |
| `run_mrbayes.log`, `history/` | Current and previous command logs/scripts |

With one run, MrBayes uses `.p` and `.t` instead of `.run1.p` and `.run1.t`.

Extend a completed result:

```bash
cladistica resume-bi --bundle output/my_run --output output/my_run_extended --ngen 4000000
```

`--ngen` is the **total target**, not the number of additional generations. It must exceed the checkpoint generation and be a multiple of 1,000. The destination must be a new directory. The command checks the original script hash, recorded file hashes, checkpoint/sample IDs, sample presence and MrBayes version, then copies the bundle and changes only the target generation and `append=yes`. Previous logs and scripts remain in `history/`. The source is unchanged. Use the same MrBayes installation when possible; the CLI requires the same version string.

If a high-level run stops before packaging:

```bash
cladistica resume-bi --bundle output/my_run/_work/trees/03_mrbayes --output output/recovered --ngen 4000000
```

If the standalone `tree` command stops, use its `03_mrbayes/` directory. Stop the source process before resuming. An active local run is rejected. A run whose manifest still says `running` after a hard kill can be recovered on its original host when its recorded process is absent; recovery is refused while that PID exists. A running-status bundle copied from another host must first be recovered on the original host. Ordinary Ctrl-C records an `interrupted` state and is portable once the required files are complete.

A crash before the first checkpoint, or during checkpoint replacement, can leave no usable current `.ckp`. The backup is retained for expert recovery but Cladistica does not silently substitute it. Missing or damaged checkpoint/tree files require a verified backup or a new analysis. Old 0.1.x flat results that lack `.ckp` and `.t` cannot be resumed from `.p` and the consensus tree alone.

If you change sequences, taxa, partitions, priors, chains or models, start a **new** analysis from the changed data. Do not edit the manifest to bypass an integrity error. Historical samples are meaningful only under the original analysis definition. See the [official MrBayes manual](https://nbisweden.github.io/MrBayes/manual.html) for its checkpoint and append semantics.

## Output and storage

High-level commands write dated directories under `output/` unless `--output` is given. Existing nonempty analysis directories are rejected rather than cleared. Important files include:

- `concatenated.fasta`, `partitions.txt`: final aligned matrix and partitions.
- `ML.tre`, `BI.tre`: ML and BI consensus trees when those stages run.
- `BI.nex`, `run1.p`, `run2.p`: convenient copies for inspection.
- `BI_resume/`: the authoritative portable BI continuation bundle.
- `_work/`: raw stages retained even after failure; includes alignment QC and full inference outputs.
- `summly.txt`, `run.log`: legacy-named summary and combined log.

`_work/` and `BI_resume/` duplicate BI files after a successful workflow. Budget disk space for both copies, and keep a backup of the full result. Deleting `_work/` after validating and backing up a complete `BI_resume/` is a manual storage choice; Cladistica never performs that cleanup automatically.

`package --run-dir RAW_RUN --output NEW_DIRECTORY --archive results.zip` collects a staged run and its BI bundle. Output cannot replace the source directory or its parent. Zip archives can be large because tree and parameter samples are included.

## Troubleshooting

| Message or symptom | Action |
| --- | --- |
| Program not found | Install that external tool and pass its full executable path. |
| Output directory is not empty | Choose a new output; use `resume-bi` for MCMC continuation. |
| Checkpoint/sample missing or changed | Restore the intact complete bundle from backup. |
| Source still active | Let it finish or stop it before continuing. |
| Version mismatch | Use the original MrBayes version. |
| Duplicate/empty sample ID | Correct FASTA headers or selected-table IDs before concatenating. |
| Tree exists but ESS is poor | Review diagnostics and extend only after confirming the model/data are suitable. |

Include the version, OS, command, relevant log and a minimal synthetic reproducer when reporting a bug. Do not publish API keys or private specimen data. See [AUDIT](AUDIT.md), [ROADMAP](ROADMAP.md) and [CITATIONS](../CITATIONS.md).
