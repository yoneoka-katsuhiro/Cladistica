# Cladistica 日本語マニュアル

[概要](../README-ja.md) · [English](MANUAL.en.md)

## 導入

配布ZIPにはCladistica本体とPython依存が含まれます。全体を展開し、実行ファイルの隣に`_internal`を残してください。以下の例で先頭の`cladistica`は、macOS配布版では`./cladistica`、Windows PowerShellでは`.\cladistica.exe`に置き換えます。コマンドは1行ずつ実行してください。

ソース版にはPython 3.10以降が必要です。仮想環境を作成・有効化した後、リポジトリ内で`python -m pip install .`を実行します。既存のmacOS向け`setup_mac.sh`、`run_cladistica.sh`も使用できます。インストール済みCLIは現在の作業フォルダに既定の出力を作り、Pythonのインストール先へは保存しません。

外部ソフトは必要な工程に応じて別途導入してください。

| 工程 | ソフト | パス指定 |
| --- | --- | --- |
| 整列 | [MUSCLE](https://drive5.com/muscle5/) | `--muscle /path/to/muscle` |
| 最尤法 | [IQ-TREE](https://iqtree.github.io/) | `--iqtree /path/to/iqtree3` |
| ベイズ法 | [MrBayes](https://nbisweden.github.io/MrBayes/) 3.2以降 | `--mrbayes /path/to/mb` |

Windowsで空白を含むパスは`--mrbayes "C:\Tools\MrBayes\mb.exe"`のように引用符で囲みます。各ソフトの実行に必要なライブラリも、そのソフトの説明に従って導入してください。Cladisticaは外部ソフトを同梱・自動導入しません。`--skip-ml`または`--skip-bi`で該当解析を省略できます。`--threads`はIQ-TREEの設定です。MrBayesのMPIプロセス数を設定する機能ではなく、現在のCLIはMrBayes実行ファイルを直接起動します。

GenBank検索・取得には`--email your.email@example.com`を指定するか、作業フォルダの`.env`に`NCBI_EMAIL`を設定します。APIキーを使う場合は`NCBI_API_KEY`に保存してください。`.env`は公開しないでください。ローカルFASTAだけの解析はNCBIへ接続しません。

## 最初の動作確認

```bash
cladistica --version
cladistica self-check
cladistica --help
cladistica tree --input examples --output preview --dry-run
```

最後のコマンドは外部ソフトを実行せず、NEXUS、partition、統計ファイルを作ります。`examples/tiny.fasta`は4分類群・32塩基の人工配列で、研究データではありません。

## GenBank候補の収集と選択

```bash
cladistica survey --genus Hymenasplenium --outgroup "Asplenium setoi" --markers rbcL atpA atpB --email your.email@example.com --output survey_run
```

`accession_all.csv`に候補と自動推薦が記録されます。元のGenBankレコード、分類群名、標本・voucher、抽出領域、配列長、QCを確認してから採用してください。`survey`は既定で検索結果を全件取得します。試行時には`--retmax 100`で制限できます。全自動の`run`では既定値が異なり、検索名ごとに500件です。

`accession_selected.csv`は1標本1行とし、`taxon`または`sample_id`と、各マーカーの列を設けます。1セルに入るアクセッションは1つまでです。空欄は欠損マーカーを示します。明示的な`sample_id`を推奨します。

```csv
sample_id,taxon,voucher_or_isolate,rbcL,matK
specimen_001,Example species,V001,ACCESSION_1,ACCESSION_2
```

上記アクセッションは書式を示す仮名です。実行には使わないでください。`sample_id`を省略した場合、Cladisticaはtaxonと明示されたvoucherからIDを作ります。マーカーごとのアクセッション番号はIDに含めず、同じ行の配列を同じ標本として扱います。重複IDや正規化後に衝突するIDは拒否します。同種の別標本には異なるIDを付けてください。IDに`|`や制御文字は使えません。

```bash
cladistica resume --accession-selected input/accession_selected.csv --accession-all survey_run/accession_all.csv --markers rbcL atpA atpB --email your.email@example.com --output selected_run
```

`--accession-all`は任意で、指定すると元の候補一覧も成果物に残ります。`run --genus ...`は自動選択から系統樹推定まで連続実行します。標本の対応が不確かな場合は、推定前に手動選択を確認してください。

## 自分の配列を追加する

未整列データは`rbcL.fasta`、`matK.fasta`など、マーカーごとに1ファイルを作ります。FASTAヘッダーの最初の`|`より前がsample IDです。全マーカーで同じ標本に同じIDを使ってください。重複ヘッダー、空配列、不正なFASTAはエラーになります。

```bash
cladistica resume --fasta-dir input/markers --markers rbcL matK --output my_run
```

アクセッションまたはマーカーFASTAからの解析に`--add-fasta-dir input/additional`を付けると、自分の配列を統合できます。同じsample ID・同じマーカーが両方にある場合は、**追加ディレクトリ側の配列が優先**されます。意図した置換か確認してください。`--markers`には使用するマーカーを明示します。任意のファイル名からすべてのマーカーを自動推定する仕様ではありません。

整列・連結済み行列から始める場合は次を使います。

```bash
cladistica resume --concatenated-fasta input/matrix.fasta --partition-file input/partitions.tsv --output matrix_run
```

partitionは`marker,start,end`列を持つCSV/TSV、IQ-TREE形式（`DNA, rbcL = 1-100`）、NEXUS charset（`charset rbcL = 1-100;`）を使えます。座標は1始まり・終点を含みます。省略時は全行列を1区画として扱います。配列長は全標本で揃える必要があります。解析前に領域境界も確認してください。

## 解析方法と解釈

マーカー抽出はGenBank注釈と別名対応表に基づきます。CDSでは翻訳に基づくQCを行い、非コード領域には終止コドン検査を適用しません。イントロン・遺伝子間領域には識別可能な境界が必要です。注釈のないプラストームへの自動注釈や類似性検索は行いません。

MUSCLEはマーカーごとに整列します。連結処理はsample IDを照合し、欠損マーカーを`?`で埋めます。この処理は、異なる遺伝子が同じ進化史を持つことを保証しません。paralogy、orthology、組換え、整列の妥当性、欠損割合を研究目的に応じて評価してください。

MLはIQ-TREEのpartition付き`-m MFP+MERGE`と、**通常のノンパラメトリックbootstrap**（`-b`、既定1,000回）です。ultrafast bootstrapではありません。IQ-TREE 2または3を推奨します。モデル選択結果と生出力は`_work/trees/`に残ります。

BIはGTR+Gamma（`nst=6 rates=gamma`）とDirichlet(1,1,1,1)の塩基頻度事前分布を使います。複数区画の場合は塩基頻度・交換率・Gamma形状をunlinkし、区画間速度を可変にします。その他の事前分布にはMrBayes既定値を使います。既定値は独立2 run、各4 chain、100万世代、1,000世代ごとのサンプリング・checkpoint、5,000世代ごとの診断、25%の相対burn-inです。`tree`では`--nruns`、`--nchains`、`--burninfrac`を指定できます。上位ワークフローは既定値を使います。burn-inは`sump`と`sumt`にも明示的に適用します。MLのモデル選択結果をBIへ自動転用しません。

計算終了は収束を意味しません。独立run間の一致、尤度トレースの定常性、ESS、MrBayes診断値、burn-in・世代数への感度を確認してください。2,000世代の例題は動作確認用です。outgroupを取得対象に加えても出力樹は自動rootingされません。樹形の解釈とrootingは利用者が確認してください。

## BIの再開と延長

`BI_resume/`はフォルダ全体を保管してください。元のbasenameはMrBayesの再開処理に必要です。

| ファイル | 役割 |
| --- | --- |
| `mrbayes_analysis.nex` | データ、モデル、事前分布、MCMC設定 |
| `.nex.ckp`と`.nex.ckp~` | 現在とバックアップのcheckpoint |
| `.nex.runN.p` | runごとのパラメータサンプル |
| `.nex.runN.t` | runごとの樹形サンプル |
| `.nex.mcmc` | MCMC診断 |
| `resume_manifest.json` | バージョン、ハッシュ、状態、継続履歴 |
| `run_mrbayes.log`、`history/` | 現在・過去のログと実行スクリプト |

1 runの場合、MrBayesは`.run1.p/.run1.t`ではなく`.p/.t`という名前を使います。

完了した解析を延長する例:

```bash
cladistica resume-bi --bundle output/my_run --output output/my_run_extended --ngen 4000000
```

`--ngen`は**追加世代数ではなく合計の目標世代数**です。checkpointの世代より大きく、1,000の倍数にしてください。出力先は未作成のディレクトリにします。Cladisticaは入力スクリプトと記録済みファイルのハッシュ、checkpointとサンプルのID、必要ファイル、MrBayesバージョンを確認してコピーします。その後、目標世代数と`append=yes`だけを変更します。元のbundleは変更せず、旧ログ・スクリプトは`history/`にも残します。可能なら元と同じMrBayes実体を使用してください。CLIは少なくともバージョン文字列の一致を要求します。

上位ワークフローが成果物作成前に中断した場合:

```bash
cladistica resume-bi --bundle output/my_run/_work/trees/03_mrbayes --output output/recovered --ngen 4000000
```

単独の`tree`コマンドなら、その出力内の`03_mrbayes/`を指定します。元の解析を停止してから再開してください。実行中のローカルrunは拒否します。強制終了でmanifestが`running`のまま残った場合、元のホストで記録PIDが存在しなければ再開できます。そのPIDが存在する間は拒否します。別ホストへコピーした`running`状態のbundleは、まず元ホストで復旧してください。通常のCtrl-Cでは`interrupted`が記録され、必要ファイルが揃えば移動後も再開できます。

最初のcheckpoint作成前やcheckpointの書換え途中に停止すると、正常な現在の`.ckp`がない場合があります。バックアップは専門的な復旧のために保存しますが、Cladisticaは自動で差し替えません。checkpointや樹形サンプルを失った場合は、確認済みのバックアップか新規解析が必要です。旧0.1.xの`.p`や合意樹だけから、失われた`.ckp`と`.t`を復元することはできません。

配列、分類群、区画、事前分布、chain、モデルを変更した場合は、変更後のデータから**新規解析**してください。整合性エラーを回避するためにmanifestを書き換えないでください。元のサンプルを再利用できるのは元と同じ解析定義の場合です。checkpointとappendの仕様は[MrBayes公式マニュアル](https://nbisweden.github.io/MrBayes/manual.html)を参照してください。

## 出力と保存容量

上位コマンドは、`--output`がなければ作業フォルダの`output/日付_番号/`へ保存します。空でない既存の解析ディレクトリは消去せず拒否します。

- `concatenated.fasta`、`partitions.txt`: 整列・連結行列と区画。
- `ML.tre`、`BI.tre`: 実行した方法の系統樹。
- `BI.nex`、`run1.p`、`run2.p`: 閲覧しやすいコピー。
- `BI_resume/`: BI再開に使う一式。
- `_work/`: QC、工程別の入力・出力、解析ログ。失敗時も保持。
- `summly.txt`、`run.log`: 互換性のため従来名を残した要約と統合ログ。

正常終了後は`_work/`と`BI_resume/`の両方にBIファイルが存在します。両方の容量を見積もり、結果一式をバックアップしてください。完全な`BI_resume/`を検証・バックアップした後に`_work/`を削除するかどうかは利用者が判断します。自動削除は行いません。

`package --run-dir RAW_RUN --output NEW_DIRECTORY --archive results.zip`は工程別出力とBI一式をまとめます。元のディレクトリやその親への上書きは拒否します。樹形・パラメータサンプルを含むため、ZIPは大きくなることがあります。

## 困ったとき

| 状況 | 対処 |
| --- | --- |
| 外部ソフトが見つからない | 対象ソフトを導入し、実行ファイルのフルパスを指定します。 |
| 出力先が空でない | 新しい出力先を指定します。MCMC継続には`resume-bi`を使います。 |
| checkpoint/sample欠損・改変 | 完全なbundleをバックアップから復元します。 |
| 元のrunが実行中 | 完了を待つか停止してから再開します。 |
| MrBayesバージョン不一致 | 元と同じバージョンを使います。 |
| sample IDが重複・空 | FASTAまたは選択表を修正し、標本とマーカーの対応を確認します。 |
| 樹はできたがESSが低い | 診断とモデル・データの妥当性を確認してから延長を検討します。 |

不具合報告にはバージョン、OS、コマンド、関連ログ、最小限の人工データを添えてください。APIキーや非公開標本データは公開しないでください。[監査](AUDIT.md)、[今後の機能案](ROADMAP.md)、[引用](../CITATIONS.md)も参照できます。
