# Audit and validation — Cladistica 0.2.0

Audit date: 2026-11-04–05. Scope: the public Cladistica repository, with priority on preserving and extending Bayesian analyses. The long-running research analysis outside this repository was inspected as a read-only reference and was not modified.

## Confirmed defects and corrections

| Finding | Consequence | Correction and regression evidence |
| --- | --- | --- |
| `TemporaryDirectory` plus selective result copying dropped `.ckp` and `.t` | Completed/interrupted high-level BI could not be continued | Durable `_work`, complete `BI_resume`, actual MrBayes extension and interruption tests |
| No MCMC continuation command or input/sample integrity check | Users had to restart or edit NEXUS manually | `resume-bi`, new destination, manifest hashes, sample/checkpoint identity and version checks |
| Packaging reset could target the source or another existing result | Analysis files could be erased | Reject source/ancestor replacement and nonempty destination |
| Burn-in set for MCMC diagnostics but not explicitly for summaries | A nondefault requested burn-in could differ from `sump`/`sumt` | Explicit relative burn-in on both summary commands |
| Duplicate FASTA/header or alignment sample IDs silently overwrote/skipped data | Changed specimen content without a clear error | Explicit validation; nine malformed/duplicate/empty input regression tests plus three within-source alignment identity tests |
| Taxon-only selected rows derived IDs from each accession | One specimen became separate rows across loci | Taxon/voucher-based stable IDs, collision checks; seven identity regression tests |
| Installed/frozen CLI default paths could point into the package | Output could be unwritable or placed in an unexpected directory | Use the current working directory |

Changes follow [ai-biocode-kata](https://github.com/dritoshi/ai-biocode-kata/tree/main/chapters), especially design (§1), testing (§8), debugging (§9), deliverables (§10), documentation (§18), and the final checklist. The BI lifecycle is isolated in one module; regression tests target data-loss/identity failures rather than mirroring helper implementation. Existing biological algorithms were not rewritten to claim improved scientific validity.

Additional release review found and fixed three reproducible failures: an arbitrary existing `--archive` path could be truncated; a parent killed with SIGKILL could leave a live MrBayes child whose samples were still writable; and an unclosed read-only Excel workbook caused Windows frozen self-check cleanup to fail with WinError 32. Regression cases now cover archive preflight and concurrent creation, a real orphan child, stopped-parent/child checkpoint recovery, and populated/empty/failed Excel reads with retained workbook references. Archive publication is complete-file, atomic and exclusive; unsupported hard-link filesystems fail without replacing user data.

## Verification

- Local verification: all 113 tests passed with real MrBayes tests enabled (2026-11-05); Python compilation also passed.
- Offline tests cover accession selection, marker extraction, reverse-strand/circular boundaries, FASTA identity, concatenation, packaging, progress, input validation and BI bundle checks.
- Real MrBayes 3.2.7a on macOS arm64: a synthetic four-taxon matrix completed at 2,000 generations, continued to 4,000, and retained generations 0/1000/2000/3000/4000 without duplicated sample rows. The source bundle hashes remained unchanged.
- Real interruption test: a long target was interrupted after checkpoint production; a new run resumed from the retained checkpoint and reached its new target while retaining the earlier sample series.
- The generated README image is from an actual synthetic run and extension. `docs/images/demo-provenance.json` records its source/settings.
- Local macOS arm64 PyInstaller package built successfully; frozen CLI `--version`, GenBank/Excel `self-check`, `tree --dry-run`, actual BI execution and 2,000→4,000 continuation passed. The local development bundle uses Python 3.13.9; release CI is fixed to Python 3.12.
- Native build workflow tests source on macOS/Windows and checks the frozen CLI's version and NEXUS generation. The release page/Actions run records determine which platform artifacts actually completed; a configured workflow alone is not a passed platform test.

Run locally:

```bash
python -m unittest discover -s tests -v
CLADISTICA_TEST_MRBAYES=1 python -m unittest discover -s tests -p test_bayes.py -v
python -m compileall -q cladistica tests packaging
```

The environment-variable example is Bash syntax. In PowerShell set `$env:CLADISTICA_TEST_MRBAYES='1'` before running Python. Integration tests require an installed `mb` and are skipped otherwise.

## Limits that remain

- Short synthetic runs test software behavior, not statistical adequacy. No empirical phylogeny or biological conclusion was revalidated by this software audit.
- ML/BI completion is not convergence. ESS, split frequencies, stationarity, model adequacy, orthology and alignment quality still require review.
- Full native Windows inference with independently installed MUSCLE/IQ-TREE/MrBayes has not been established by local macOS tests. External binaries are not bundled.
- Crash recovery requires a valid checkpoint and matching samples. A crash during checkpoint replacement may require expert use of the preserved backup; the program does not auto-promote `.ckp~`.
- A hard-killed run with `running` status can be recovered only on its original host after both recorded parent and child PIDs are absent. Missing/invalid PID metadata and PID reuse conservatively block automatic recovery; a dead launcher alone does not authorize copying samples.
- Old results that already lost checkpoint/tree samples are not recoverable from traces alone. No software change can reconstruct the absent chain state.
- Retained `_work` plus `BI_resume` increases disk usage. Automatic cleanup and diagnostic dashboards remain proposals in [ROADMAP](ROADMAP.md).

## 日本語要約

今回の主要な修正は、BIのcheckpointと樹形サンプルが一時ディレクトリの消去とともに失われる問題です。必要ファイル一式を永続化し、元結果を保持したまま新しい出力先へ継ぎ足すコマンドを追加しました。実MrBayesで正常終了後と中断後の再開を確認しています。

併せて、FASTAの重複・空配列の見落とし、マーカーごとにsample IDが分裂する問題、成果物の上書き、burn-in適用範囲を修正しました。短い人工データ試験はソフトウェアの挙動確認であり、収束や生物学的な妥当性を保証するものではありません。英日マニュアルには解析方法、再開条件、失われた旧ファイルの復元不能性、保存容量、外部依存を明記しました。
