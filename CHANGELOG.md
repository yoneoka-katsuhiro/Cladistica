# Changelog

## 0.2.0 — 2026-11-05

### Fixed

- Preserve MrBayes checkpoint, backup, tree samples, parameter samples, diagnostics, scripts and logs. Earlier high-level workflows deleted their temporary raw analysis directory after copying only a subset of results.
- Keep inference workspaces inside the chosen output directory, including when interrupted; package the complete `BI_resume` directory and protect existing nonempty outputs from replacement.
- Add `resume-bi` to validate and copy a stopped run before `append=yes`; reject altered inputs/samples, missing checkpoint/trees, mismatched run IDs, active local source processes, incompatible MrBayes versions and invalid total generation targets.
- Apply the selected burn-in fraction explicitly to MrBayes parameter and tree summaries.
- Prevent malformed/empty FASTA, duplicate headers and colliding alignment sample IDs from silently changing or omitting samples.
- Preserve specimen identity across loci when selected accession rows omit `sample_id`. Earlier fallback IDs incorporated the per-marker accession and split one specimen across concatenated rows. Reject normalized ID collisions and unsafe delimiters.
- Store default CLI outputs in the working directory, including installed and frozen distributions.

### Added

- English/Japanese overview and manuals, actual synthetic BI output illustration, audit report and prioritized feature proposals.
- Native macOS arm64 and Windows x64 build workflows, frozen CLI smoke checks, dependency inventory, license collection and archive checksums.
- Offline regression tests and opt-in real MrBayes completed/interrupted continuation tests.

### Compatibility

- Old results lacking `.ckp` and `.t` cannot be converted into resumable chains. Start a new analysis or recover complete original files from backup.
- `resume` continues the data-processing workflow from sequences; use `resume-bi` for MCMC continuation. `resume-bi --ngen` is the total target, not an additional count.
- New BI generation targets must be positive multiples of 1,000. Raw `_work` is now retained and requires additional disk space.

## 0.1.3

See [earlier release notes](RELEASE_NOTES.md) for the marker-resolution and workflow changes preceding this audit.
