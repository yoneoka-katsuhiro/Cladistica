# Cladistica v0.1.3 revision

Revision date: 2026-09-10. The package version is 0.1.3, as requested.

## Marker recognition

The former fixed-list rejection is replaced by a curated alias registry plus
annotation-driven resolution. Additional identifiers can be requested without
editing the Python code. Gene/product annotations, intergenic flanks, and
annotated or spliced-feature introns provide the extraction boundaries.

Reported names `rbcL-accD`, `accD`, `trnK`, and `5.8S ribosomal RNA` are
supported. `trnK` denotes the tRNA; use `trnK intron` for its intron. Broad rRNA
precursor annotations are not substituted for a specific 5.8S region.

Additional fixes made during verification:

- Wraparound extraction requires an explicitly circular record.
- Intron numbering follows strand direction; ambiguous unnumbered introns
  are rejected. Trans-spliced features are not used to infer genomic introns.
- Short identifiers and similar gene names are matched at token boundaries.
- tRNA anticodons are not mistaken for intergenic marker names.
- Accession CSV columns accept curated aliases; marker FASTA filenames are
  sanitized, and path-like marker identifiers are rejected.
- Existing progress bar and terminal resize handling are retained.
- Fable5 fixes are retained: MrBayes autoclose, partition rate variation, and
  IQ-TREE v1 thread-flag compatibility. Standard bootstrap `-b` is retained.

## Verification

`python -m unittest discover -s tests -v`: 69 tests passed, including the
existing progress/resize tests and new extraction regression tests.

Live NCBI extraction was repeated with the final implementation:

| Requested marker | Accession | Extracted bp | Method |
| --- | --- | ---: | --- |
| rbcL-accD | PQ187915.1 | 548 | Flanking features |
| accD | PQ187915.1 | 807 | CDS feature; coding QC passed |
| trnK | PZ094978.1 | 72 | tRNA feature |
| 5.8S ribosomal RNA | JX475167.1 | 160 | rRNA feature |
| trnK intron | EU161480.1 | 2422 | Intron feature |

A separate execution smoke test used six synthetic samples, two loci
(`accD`, `rbcL-accD`), two standard bootstrap replicates, and 10,000 MrBayes
generations. Actual MUSCLE, IQ-TREE 3, and MrBayes binaries completed the
workflow. Output included both trees, model summary, and both parameter traces
(`run1.p`, `run2.p`). This checks execution, not biological convergence.

## Limits

Marker recognition does not establish orthology and does not infer missing
annotations by sequence similarity. GenBank records without reliable region
boundaries can be rejected with reasons in the accession table/log. Review
paralogs, conflicting annotations, genomic compartment, and alignment before
using results for research. The default marker panel has not been expanded
automatically; request extra markers with `--markers`.

Annotation conventions:
[NCBI feature table](https://www.ncbi.nlm.nih.gov/genbank/feature_table/) and
[INSDC feature table](https://www.insdc.org/submitting-standards/feature-table/).
