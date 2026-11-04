# Feature ideas / 今後の機能案

These are proposals, not implemented capabilities. They build on the existing CLI and do not change the interpretation of current results.

| Priority | Proposal | User benefit | Validation needed |
| --- | --- | --- | --- |
| High | BI diagnostics report with per-run traces, ESS, PSRF and split-frequency trends | See why a run needs extension and review all diagnostics together | Compare against independent Tracer/MrBayes outputs; test short and nonstationary chains |
| High | `inspect-bi` with checkpoint generation, size, integrity and recovery readiness | Check a bundle before moving it to another machine | Missing, corrupt and partially written files; stopped/running processes on each OS |
| High | Explicit reproducible seeds and a dataset/model manifest for ML and BI | Repeat a defined analysis and audit software settings | Known fixtures, version differences, stochastic behavior documented |
| Medium | Alignment viewer/report with marker occupancy, flagged records and removed columns | Review biologically questionable regions before inference | Codon-aware fixtures, ambiguous states and sample identity tests |
| Medium | Explicit MAFFT backend in addition to MUSCLE | Choose an alignment method suited to the dataset | Backend-specific options/version detection and expected alignment fixtures |
| Medium | Safer manual checkpoint-backup recovery with a preview | Recover from a crash during `.ckp` replacement without editing files by hand | Confirm backup/sample compatibility and recovery on real MrBayes versions |
| Medium | Separate model/prior configuration for BI, with provenance | Make methodological decisions explicit without editing generated NEXUS | Parameter validation, partition tests and refusal to append under changed settings |
| Later | Workflow manifest, storage estimates and verified archive cleanup | Reduce duplicated storage while proving a portable bundle is complete | Fault injection and checksum validation before any deletion |
| Later | Guided desktop interface over the same tested CLI | Help occasional users select inputs and interpret errors | Real macOS/Windows usability tests, cancellation and recovery behavior |

## 日本語

最優先の案は、BI診断を1画面にまとめるレポートと、再開bundleの状態確認コマンドです。ESSやPSRFなどの値だけで機械的に「収束済み」と断定せず、独立runの比較やトレースを研究者が確認できる構成が有用です。

次に、seedとモデル設定の明示、マーカー欠損・整列QCの可視化、MAFFTの選択肢が考えられます。保存容量を削減する自動削除やcheckpointバックアップへの自動切替は、復旧性を損なわない検証が必要なため現版には追加していません。GUIを設ける場合も、検証済みCLIを共通の実行基盤にする方針です。
