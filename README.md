# Cladistica v0.1.0

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
Cladistica displays an animated growing flower followed by the pipeline status:

```text
    -- * --
      \|/
      / \
  Cladistica is growing
[OK] NCBI records and representative selection - 28 selected samples
[>>] MUSCLE marker alignments - matK: 24 sequences
[  ] Concatenated alignment
[  ] ModelFinder and ML tree
[  ] MrBayes two-run analysis
[  ] Final output
```

The active marker, bootstrap count, and MrBayes generation count are shown.
Completed stages remain marked `[OK]`. When output is redirected to a file,
Cladistica automatically switches to plain status lines without terminal
control codes. Add `--no-progress` to disable the display.

## Full Workflow

```bash
bash run_cladistica.sh run \
  --genus Hymenasplenium \
  --outgroup "Asplenium setoi" "Asplenium nidus" \
  --email "your.email@example.com" \
  --bootstrap 1000 \
  --ngen 1000000
```

The default selection is one sample per taxon. Ranking prioritizes extracted
marker coverage, clean or CDS-QC-passing markers, publication evidence, and
total extracted sequence length.

## Recommended Usage Examples

The same usage captions are available in the terminal:

```bash
bash run_cladistica.sh examples
```

### 1. Start From Self-Made FASTA, Then Run ML And BI

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

### 2. Create Only `accession_all.csv` To Inspect NCBI Data

```bash
bash run_cladistica.sh survey \
  --genus Hymenasplenium \
  --outgroup "Asplenium setoi" "Asplenium nidus" \
  --email "your.email@example.com"
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

### 3. Manually Select Accessions, Then Resume To Tree Inference

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
  --email "your.email@example.com" \
  --bootstrap 1000 \
  --ngen 1000000
```

Cladistica validates the table, downloads and extracts the selected marker
sequences, aligns and concatenates them, and runs ML and BI. Supplying
`--accession-all` is optional, but preserves the full survey beside the final
selected table.

### 4. Combine Cladistica FASTA With Your Own Research Sequences

Put the new sequences in a second marker-wise directory:

```text
input/my_sequences/
  rbcL.fasta
  matK.fasta
  trnL-F.fasta
```

Use the completed Cladistica output directory as the primary FASTA source. If
you are using a low-level `accessions` output directory instead, point
`--fasta-dir` to its `fasta_by_marker/` directory.

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

## Citation Guidance

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

Use `package` to flatten an existing working directory:

```bash
bash run_cladistica.sh package \
  --run-dir /path/to/old/run \
  --output output/20260726_1 \
  --archive Cladistica_results.zip
```
