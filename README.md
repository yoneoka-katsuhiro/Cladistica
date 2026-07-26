# Cladistica v0.1.1

Cladistica builds chloroplast DNA datasets and phylogenetic trees for
systematics and biodiversity research. It can start from
NCBI GenBank searches, curated accession tables, marker-wise FASTA files, or an
already concatenated alignment.

## What It Does

| Step | Main outputs |
| --- | --- |
| Survey GenBank records | `accession_all.csv`, `summly.txt`, `run.log` |
| Select representative accessions | `accession_selected.csv` |
| Download and extract markers | marker-wise FASTA files |
| Align and concatenate markers | `concatenated.fasta`, `partitions.txt`, `BI.nex` |
| Infer trees | `ML.tre`, `BI.tre`, `run1.p`, `run2.p` |

Cladistica queries NCBI Entrez directly. Tree workflows use MUSCLE, IQ-TREE /
ModelFinder, and MrBayes when those stages are requested.

## Requirements

- macOS or Linux
- Python 3.10 or newer
- A contact email for NCBI Entrez
- MUSCLE, IQ-TREE, and MrBayes for full tree workflows

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

## Recommended First Run

The full workflow can become slow when all default markers are used for a genus
with many accessions. Start with two chloroplast markers, then expand the marker
set after confirming that the survey, extraction, alignment, and tree steps work
for your group.

```bash
bash run_cladistica.sh run \
  --genus Hymenasplenium \
  --outgroup "Asplenium setoi" "Asplenium nidus" \
  --markers rbcL trnL-F \
  --email "your.email@example.com" \
  --bootstrap 1000 \
  --ngen 1000000
```

If `NCBI_EMAIL` is set in `.env`, the `--email` option can be omitted.

## Common Workflows

### Survey NCBI Records Only

Use this when you want to inspect all candidate records before deciding which
accessions should be used.

```bash
bash run_cladistica.sh survey \
  --genus Hymenasplenium \
  --outgroup "Asplenium setoi" "Asplenium nidus" \
  --markers rbcL trnL-F \
  --email "your.email@example.com"
```

This stops after writing `accession_all.csv`, `summly.txt`, and `run.log`.

### Resume From A Manually Curated Accession Table

First run `survey`, edit `accession_all.csv`, and save the selected rows as
`accession_selected.csv`. Then resume from that table:

```bash
bash run_cladistica.sh resume \
  --accession-selected input/accession_selected.csv \
  --accession-all output/20260726_1/accession_all.csv \
  --email "your.email@example.com" \
  --bootstrap 1000 \
  --ngen 1000000
```

The selected table should contain one row per sample, no more than one accession
per marker cell, and marker column names that match Cladistica marker names.

### Start From Your Own Marker FASTA Files

Use this when the FASTA files were prepared outside Cladistica. Keep one file
per marker, for example `rbcL.fasta`, `matK.fasta`, and `trnL-F.fasta`.

```bash
bash run_cladistica.sh resume \
  --fasta-dir input/my_fasta_by_marker \
  --markers rbcL matK trnL-F \
  --bootstrap 1000 \
  --ngen 1000000
```

Cladistica aligns each marker with MUSCLE, concatenates the alignments, creates
partitions, and runs ML and BI unless those stages are skipped.

### Start From A Concatenated Alignment

Use this when alignment and concatenation are already complete.

```bash
bash run_cladistica.sh resume \
  --concatenated-fasta input/concatenated.fasta \
  --partition-file input/partitions.txt \
  --bootstrap 1000 \
  --ngen 1000000
```

Without `--partition-file`, the concatenated alignment is analyzed as one cpDNA
partition.

### Combine Cladistica FASTA With New Research Sequences

Use this when you want to combine public GenBank-derived sequences with newly
generated sequences from your own study.

```bash
bash run_cladistica.sh resume \
  --fasta-dir output/20260726_1 \
  --add-fasta-dir input/my_sequences \
  --markers rbcL matK trnL-F \
  --bootstrap 1000 \
  --ngen 1000000
```

The part of each FASTA header before the first `|` is treated as the sample ID.
If the same sample ID appears in both inputs for the same marker, the sequence
from `--add-fasta-dir` is used.

Run `bash run_cladistica.sh --help` or `bash run_cladistica.sh examples` for
additional command-line guidance.

## Output

Each high-level `run`, `resume`, or `survey` command creates a dated output
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

Open both `run1.p` and `run2.p` in Tracer. Check trace stationarity and ESS for
each parameter before using the BI tree in a publication.

## Citations

Analyses produced with Cladistica may use public sequence records and external
phylogenetic software. For manuscripts, cite the sequence accessions or datasets
you used and the external tools actually run by your analysis. See
`CITATIONS.md` for suggested references.

## Low-level Commands

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
