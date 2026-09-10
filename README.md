# Cladistica v0.1.3

Cladistica is a collection of modular plant molecular phylogenetic pipelines. It can run
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

## Flexible marker recognition

`--markers` is not restricted to a fixed list. Cladistica combines a curated
registry of common plant loci with annotation-driven resolution:

```bash
bash run_cladistica.sh survey \
  --genus Hymenasplenium \
  --markers rbcL-accD accD trnK "trnK intron" "5.8S ribosomal RNA" ITS2 PHYC
```

- Common aliases resolve to stable names and biological types. For example,
  `5.8S ribosomal RNA`, `5.8S rRNA`, and `5.8S rDNA` resolve to `5.8S`.
- Unregistered names such as `PHYC` are accepted and matched against GenBank
  `gene`, `product`, `standard_name`, `label`, and `note` qualifiers.
- A two-locus name such as `rbcL-accD`, `atpB-rbcL`, or
  `trnL(UAA)-trnF(GAA)` is interpreted as a flanking-feature interval. On an
  annotated long genome, Cladistica locates both features and extracts the
  shortest intervening region, including intervals that cross a circular
  sequence origin.
- An intron request such as `trnK intron`, `trnL intron`, or
  `rps16 intron 1` is resolved from an annotated `intron` feature or from the
  gap between parts of a spliced CDS or RNA feature. Intron numbers follow
  transcriptional order, including on the reverse strand. When several
  introns match, specify a number; no longest-intron guess is made.
  This is kept separate from
  the mature `trnK` tRNA so that a 72 bp tRNA is not aligned against a roughly
  2.4 kb intron-containing region.
- The GenBank feature type determines whether the sequence is treated as CDS,
  rRNA, tRNA, or noncoding. Frame, stop-codon, and codon QC is applied only to
  a CDS. A partial Sanger CDS does not need to contain its start or terminal
  stop codon.

Short marker-specific GenBank records can still be accepted as whole records.
For long genomes, Cladistica deliberately refuses a requested marker when the
necessary annotation cannot be identified; it never silently substitutes the
whole plastome. Such cases are recorded with an extraction-failure reason in
`accession_all.csv` and `run.log`. Because public records vary in annotation
quality and naming, inspect these files before publication-scale analysis.

The default marker panel is unchanged. Pass `--markers` explicitly to request
additional loci. For supplied FASTA, use canonical names such as `5.8S.fasta`
and `trnK_intron.fasta`; spaces and punctuation in output filenames are
converted to underscores. Accession CSV columns also accept curated aliases
such as `5.8S ribosomal RNA`.

Recognition uses annotation, not a homology search or orthology inference.
Unannotated subregions, paralogs, and conflicting annotations still require
curation. An rRNA precursor listing 5.8S among several regions is not itself
a valid 5.8S feature. Nuclear and organellar markers should be combined only
when the study design supports that choice. See `RELEASE_NOTES.md` for the
verified examples and revision tests.

## Live progress

During `run`, `survey`, `resume`, and the individual analysis commands,
Cladistica displays a single progress-bar line for the active stage:

```text
[================================] 100% ModelFinder
[========================--------]  75% Bootstrap replicates - 750/1000 replicates
[========------------------------] ~25% ML tree search - candidate trees
[================================] 100% MrBayes MCMC - 1,000,000/1,000,000 generations
```

Long combined labels are split into separate `ModelFinder`, `Bootstrap
replicates`, `ML tree search`, `MrBayes MCMC`, and `BI summary and consensus`
stages. Bootstrap and MCMC percentages come directly from replicate and
generation numbers in the external program output. A leading `~` marks a
phase-based estimate where the final number of iterations is not known in
advance. There is no animation or background drawing thread. On an interactive
TTY, Cladistica reads the current terminal width for every redraw, erases the
active line with ANSI erase-line, and keeps output one column short of the
right edge to prevent automatic wrapping. `SIGWINCH` triggers an immediate
resize redraw on macOS and Linux when signal handling is available. Redirected
output uses the existing plain status lines without ANSI control codes. Add
`--no-progress` to disable progress output.

## Full workflow

Start with a small marker panel before expanding to a heavier genus-wide
analysis.

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

## Partial-workflow captions

The same captions are available in the terminal:

```bash
bash run_cladistica.sh examples
```

### Caption 1: Use self-made FASTA, then run ML and BI together

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

### Caption 2: Inspect all data registered in NCBI

```bash
bash run_cladistica.sh survey \
  --genus Hymenasplenium \
  --outgroup "Asplenium setoi" "Asplenium nidus"
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

### Caption 3: Manually select accessions, then resume

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

### Caption 4: Combine Cladistica FASTA with new research sequences

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

Use `package` to flatten an older v0.1.0, v0.1.1, or v0.1.2 working directory:

```bash
bash run_cladistica.sh package \
  --run-dir /path/to/old/run \
  --output output/20260726_1 \
  --archive Cladistica_results.zip
```
