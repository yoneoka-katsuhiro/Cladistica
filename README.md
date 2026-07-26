# Cladistica v0.1.1

Cladistica is a collection of modular chloroplast DNA pipelines. It can run
from GenBank discovery to maximum-likelihood (ML) and Bayesian-inference (BI)
trees, or resume from an accession table, marker FASTA files, or an aligned
concatenated FASTA.

## Setup

On macOS:

```bash
bash setup_mac.sh
```

Create `.env` from the example file and set the contact address required by
NCBI Entrez:

```bash
cp .env.example .env
```

```text
NCBI_EMAIL=your.email@example.com
```

Confirm the external programs before a long analysis:

```bash
command -v muscle
command -v iqtree3 || command -v iqtree2 || command -v iqtree
command -v mb
```

## Live progress

During `run`, `survey`, `resume`, and the individual analysis commands,
Cladistica displays a fixed 40-column by 5-row sequence stream above the
pipeline status:

```text
+----------------------------------------+
|Hymenasplenium hondoense                |
|        NC_035840.1                     |
|ATGTCACCACAAACAGAGACTAAAGCAAGTGTTGGATT|
|                 rbcL                   |
|                         Asplenium setoi|
+----------------------------------------+
Sequence stream [RUNNING]
[================================] 100% ModelFinder
[========================--------]  75% Bootstrap replicates - 750/1000 replicates
[========------------------------] ~25% ML tree search - candidate trees
[--------------------------------]   0% MrBayes MCMC
[--------------------------------]   0% BI summary and consensus
```

Taxon names, sample IDs, accessions, markers, and representative 40-base
windows from the sequences actually used by the workflow move continuously
from left to right. Whole alignments are never copied into animation memory.
The vocabulary changes immediately at each stage:

- accession survey: `NCBI GenBank`, `DDBJ`, `ENA`, `INSDC`, `Entrez`, queried taxa, markers, and accessions;
- download and extraction: FASTA, CDS/noncoding, feature extraction, QC, and downloaded sequence windows;
- alignment and concatenation: `MUSCLE`, marker/sample names, aligned bases, partitions, and missing data;
- ML model selection: `JC`, `K2P`, `HKY`, `TN`, `TIM`, `TVM`, `SYM`, `GTR`, `+F`, `+I`, `+G4`, `BIC`, and models reported by IQ-TREE;
- BI and packaging: `MrBayes`, MCMC runs/chains, convergence terms, and the final output filenames.

Cladistica currently queries NCBI Entrez directly. `DDBJ` and `ENA` identify
the other INSDC partners that exchange sequence records with GenBank; their
appearance does not mean that separate DDBJ or ENA API queries were made.
When Cladistica detects an exception, failed external command, or interruption,
the visible strings break into individual characters, fall through the five
rows, and accumulate at the bottom. This also runs when the user interrupts
Cladistica with Control-C. Cladistica terminates the active external child
process before presenting the collapse effect. A silent long-running IQ-TREE or MrBayes
process remains `RUNNING`; the animation does not guess that valid computation
has stalled.

Long combined labels are split into separate `ModelFinder`, `Bootstrap
replicates`, `ML tree search`, `MrBayes MCMC`, and `BI summary and consensus`
rows. Bootstrap and MCMC percentages come directly from replicate and
generation numbers in the external program output. A leading `~` marks
phase/ETA-based estimates where the final number of iterations is not known in
advance. When output is redirected to a file,
Cladistica automatically switches to plain status lines without terminal
control codes. Add `--no-progress` to disable the display.

Preview normal flow and detected-failure behavior without starting an analysis:

```bash
bash run_cladistica.sh demo
bash run_cladistica.sh demo --fail
```

## Start Here: Lightweight Full Workflow

Running the full workflow with all default markers can take time, especially
when many accessions are found and both ML and BI are enabled. If this is your
first run, start with the two-marker example below.

```bash
bash run_cladistica.sh run \
  --genus Hymenasplenium \
  --outgroup "Asplenium setoi" "Asplenium nidus" \
  --markers rbcL trnL-F \
  --bootstrap 1000 \
  --ngen 1000000
```

The default selection is one sample per taxon. Ranking prioritizes extracted
marker coverage, clean or CDS-QC-passing markers, publication evidence, and
total extracted sequence length.

For a broader production analysis, add more markers after confirming that the
survey, sequence extraction, alignment, and tree inference steps work as
expected for your taxonomic group.

## Recommended Usage Examples

The same examples are available in the terminal:

```bash
bash run_cladistica.sh examples
```

### Example 1: Use self-made FASTA, then run ML and BI together

For unaligned FASTA files separated by marker, name each file after the marker,
for example `rbcL.fasta`, `matK.fasta`, and `trnL-F.fasta`:

```bash
bash run_cladistica.sh resume \
  --fasta-dir input/my_fasta_by_marker \
  --markers rbcL matK trnL-F \
  --bootstrap 1000 \
  --ngen 1000000
```

Cladistica performs marker-wise MUSCLE alignment, concatenation, partition
generation, IQ-TREE model selection and ML, and two independent MrBayes runs.

If the FASTA is already aligned and concatenated, MUSCLE and concatenation are
skipped:

```bash
bash run_cladistica.sh resume \
  --concatenated-fasta input/concatenated.fasta \
  --partition-file input/partitions.txt \
  --bootstrap 1000 \
  --ngen 1000000
```

`--partition-file` accepts:

- CSV or TSV with `marker`, `start`, and `end` columns;
- IQ-TREE text such as `DNA, rbcL = 1-1200`;
- NEXUS charset lines such as `charset rbcL = 1-1200;`.

Without a partition file, the concatenated alignment is analyzed as one cpDNA
partition.

### Example 2: Create only `accession_all.csv` to inspect NCBI records

```bash
bash run_cladistica.sh survey \
  --genus Hymenasplenium \
  --outgroup "Asplenium setoi" "Asplenium nidus" \
  --markers rbcL trnL-F
```

`survey` pages through all matching NCBI search results by default. For a quick
bounded check, add `--retmax 500`.

This stops after writing:

```text
accession_all.csv
summly.txt
run.log
```

`accession_all.csv` contains accepted candidates, automatically recommended
representatives, unselected candidates, and rejected records with reasons. The
`selection_status` column is a recommendation for review; `survey` does not
continue to FASTA alignment or tree inference.

### Example 3: Manually select accessions, then resume to tree inference

First run `survey`. Open `accession_all.csv`, retain the desired sample rows,
and save them as `accession_selected.csv`.

The selected table must follow these rules:

- one row per sample;
- `taxon` or `sample_id` must be present;
- one accession at most in each marker cell;
- marker column names must match Cladistica marker names.

Then run:

```bash
bash run_cladistica.sh resume \
  --accession-selected input/accession_selected.csv \
  --accession-all output/20260726_1/accession_all.csv \
  --bootstrap 1000 \
  --ngen 1000000
```

Cladistica validates the table, downloads and extracts the selected marker
sequences, aligns and concatenates them, and runs ML and BI. Supplying
`--accession-all` is optional, but preserves the full survey beside the final
selected table.

### Example 4: Combine Cladistica FASTA with your own research sequences

Put the new sequences in a second marker-wise directory:

```text
input/my_sequences/
  rbcL.fasta
  matK.fasta
  trnL-F.fasta
```

Use the completed Cladistica output as the primary FASTA source:

```bash
bash run_cladistica.sh resume \
  --fasta-dir output/20260726_1 \
  --add-fasta-dir input/my_sequences \
  --markers rbcL matK trnL-F \
  --bootstrap 1000 \
  --ngen 1000000
```

The part of each FASTA header before the first `|` is the sample ID. Use the
same sample ID across markers. If the same sample ID occurs in both inputs for
the same marker, `--add-fasta-dir` takes precedence and the replacement is
recorded in `run.log`.

The marker FASTA files returned in the final output contain the combined
dataset, including both Cladistica and user sequences.

To keep a first combined run lightweight, limit the marker set in the same way:

```bash
bash run_cladistica.sh resume \
  --fasta-dir output/20260726_1 \
  --add-fasta-dir input/my_sequences \
  --markers rbcL trnL-F \
  --bootstrap 1000 \
  --ngen 1000000
```

## Output

Every high-level `run`, `resume`, or `survey` command creates one dated
directory. Repeated runs on the same day increment the suffix:

```text
output/
  20260726_1/
  20260726_2/
```

A complete tree analysis can contain:

```text
accession_all.csv
accession_selected.csv
<marker>.fasta
concatenated.fasta
partitions.txt
BI.nex
ML.tre
BI.tre
run1.p
run2.p
summly.txt
run.log
```

Files that do not apply to the chosen starting point are omitted. For example,
a self-made concatenated FASTA analysis has no accession CSV or marker FASTA.
Intermediate alignments, QC tables, and external-tool working files stay in a
temporary directory and are summarized into `summly.txt` and `run.log`.

Open both `run1.p` and `run2.p` in Tracer. Check trace stationarity and ESS for
each parameter before using the BI tree in a publication.

## Citation guidance

Analyses produced with Cladistica may use public sequence records and external
phylogenetic software. For manuscripts, cite the sequence accessions or datasets
you used and the external tools actually run by your analysis. See
`CITATIONS.md` for suggested references.

At minimum, check whether your workflow used:

- NCBI GenBank / E-utilities for accession discovery and sequence retrieval;
- MUSCLE for marker-wise multiple sequence alignment;
- IQ-TREE and ModelFinder for model selection and maximum-likelihood inference;
- MrBayes for Bayesian inference;
- Biopython for GenBank parsing and sequence handling.

## Low-level pipelines

The individual stages remain independently callable:

```bash
bash run_cladistica.sh accessions --help
bash run_cladistica.sh download --help
bash run_cladistica.sh align --help
bash run_cladistica.sh concat --help
bash run_cladistica.sh tree --help
```

Use `package` to flatten an older v0.1.0 or v0.1.1 working directory:

```bash
bash run_cladistica.sh package \
  --run-dir /path/to/old/run \
  --output output/20260726_1 \
  --archive Cladistica_results.zip
```
