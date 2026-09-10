# Cladistica v0.1.3

Cladistica is a plant molecular phylogenetic pipeline for surveying sequences
in NCBI GenBank, selecting accessions, combining them with your own sequences,
and performing maximum-likelihood (ML) and Bayesian-inference (BI) analyses.

## Main workflow

Cladistica is designed around two common use cases:

1. Search NCBI GenBank and create `accession_all.csv` for manual review.
2. Read a curated `accession_selected.csv`, combine the selected GenBank
   sequences with your own sequences, and continue through alignment,
   concatenation, ML, and BI analyses.

## Setup

On macOS:

```bash
bash setup_mac.sh
```

Create `.env` from the example file and enter the contact email required by
NCBI Entrez:

```bash
cp .env.example .env
```

```text
NCBI_EMAIL=your.email@example.com
```

Confirm that the external programs are available before starting a long run:

```bash
command -v muscle
command -v iqtree3 || command -v iqtree2 || command -v iqtree
command -v mb
```

## Marker selection

In Cladistica v0.1.3, `--markers` can specify markers beyond predefined loci
such as `rbcL`. When sequences are retrieved automatically from GenBank, each
record must contain annotations that allow Cladistica to identify the requested
marker.

Cladistica checks the following GenBank feature qualifiers:

- `gene`
- `standard_name`
- `label`
- `locus_tag`
- `product`
- `note`
- `function`

The title or description can also identify a marker when the GenBank record is
a short, marker-specific sequence. Intergenic regions and introns require
identifiable annotated boundaries. Cladistica does not identify unannotated
regions by sequence similarity.

Specify the markers used in your study explicitly. For example:

```bash
--markers rbcL atpA atpB
```

An additional marker such as `GapCp` can be requested in the same way:

```bash
--markers GapCp
```

## 1. Create an accession survey from NCBI

Use `survey` to search GenBank and create an accession table without running a
phylogenetic analysis:

```bash
bash run_cladistica.sh survey \
  --genus Hymenasplenium \
  --outgroup "Asplenium setoi" "Asplenium nidus" \
  --markers rbcL atpA atpB \
  --email "your.email@example.com"
```

The results are written to a dated directory such as `output/20260910_1/`:

```text
accession_all.csv
summly.txt
run.log
```

`accession_all.csv` lists the GenBank candidates and Cladistica's automatic
selection recommendation. Review this table manually, retain the samples and
accessions required for your study, and save the selected rows as
`accession_selected.csv`.

For a short test search, add `--retmax 100`. Without `--retmax`, `survey`
retrieves all matching search results.

## 2. Combine selected accessions with your own sequences

Prepare the two inputs:

- A manually curated `accession_selected.csv` made from the survey results.
- One unaligned FASTA file per marker containing your own sequences.

Place your FASTA files in one directory and name each file after its marker:

```text
input/my_sequences/
  rbcL.fasta
  atpA.fasta
  atpB.fasta
```

Then run:

```bash
bash run_cladistica.sh resume \
  --accession-selected input/accession_selected.csv \
  --accession-all output/20260910_1/accession_all.csv \
  --add-fasta-dir input/my_sequences \
  --markers rbcL atpA atpB \
  --email "your.email@example.com" \
  --bootstrap 1000 \
  --ngen 1000000
```

Cladistica downloads and extracts the selected GenBank sequences, combines
them with your sequences, aligns each marker with MUSCLE, concatenates the
alignments, and runs IQ-TREE and MrBayes.

`--accession-all` is optional, but including it preserves the complete survey
table with the final results. In each FASTA header, the text before the first
`|` is treated as the sample ID. Use the same sample ID across markers. If the
same sample and marker occur in both inputs, the sequence in
`--add-fasta-dir` takes precedence.

The selected accession table must contain:

- one row per sample;
- a `taxon` or `sample_id` value for every row;
- no more than one accession in each marker cell;
- marker columns matching the names supplied to `--markers`.

## Other starting points

To analyze only your own unaligned marker FASTA files:

```bash
bash run_cladistica.sh resume \
  --fasta-dir input/my_sequences \
  --markers rbcL atpA atpB \
  --bootstrap 1000 \
  --ngen 1000000
```

To start from an aligned, concatenated FASTA file:

```bash
bash run_cladistica.sh resume \
  --concatenated-fasta input/concatenated.fasta \
  --partition-file input/partitions.txt \
  --bootstrap 1000 \
  --ngen 1000000
```

The partition file may be a CSV or TSV containing `marker`, `start`, and `end`
columns, an IQ-TREE partition file, or a NEXUS file containing charset lines.

## Output

Each `survey`, `resume`, or `run` command creates a dated output directory.
Repeated runs on the same day receive a sequential suffix:

```text
output/
  20260910_1/
  20260910_2/
```

A complete analysis can contain:

```text
accession_all.csv
accession_selected.csv
<marker>.fasta
concatenated.fasta
partitions.txt
ML.tre
BI.tre
BI.nex
run1.p
run2.p
summly.txt
run.log
```

Files that do not apply to the selected starting point are omitted.

Open both `run1.p` and `run2.p` in Tracer and evaluate trace stationarity and
effective sample size (ESS) before using the BI tree in a publication.

## Citation guidance

Analyses produced with Cladistica may use public sequence records and external
phylogenetic software. Cite the sequence accessions or datasets used in your
study and the external tools actually run by your analysis. See
`CITATIONS.md` for suggested references.

Depending on the workflow, these may include:

- NCBI GenBank and E-utilities for accession discovery and retrieval;
- Biopython for GenBank parsing and sequence handling;
- MUSCLE for multiple sequence alignment;
- IQ-TREE and ModelFinder for model selection and ML inference;
- MrBayes for BI.

## Individual commands

Each stage can also be run separately:

```bash
bash run_cladistica.sh accessions --help
bash run_cladistica.sh download --help
bash run_cladistica.sh align --help
bash run_cladistica.sh concat --help
bash run_cladistica.sh tree --help
```
