# Citation and methods reporting

Suggested software citation:

> Yoneoka, K. (2026). *Cladistica: plant marker phylogenetics with resumable Bayesian analyses* (Version 0.2.0) [Computer software]. https://github.com/yoneoka-katsuhiro/Cladistica

[CITATION.cff](CITATION.cff) supplies machine-readable metadata. No DOI has been assigned to Cladistica.

Also cite the data and programs actually used:

- **GenBank:** report accession **versions**, retrieval date, voucher/specimen assignments, and publications associated with reused sequences. Retain the selected accession table; citing Cladistica does not replace attribution to sequence producers.
- **MUSCLE:** cite the installed major version. For v5, use [Edgar (2022), *Muscle5*](https://doi.org/10.1038/s41467-022-34630-w). See the [developer's citation guidance](https://drive5.com/muscle/) for older versions. Cladistica currently uses one alignment; it does not run ensemble analyses.
- **IQ-TREE:** use the [official citations for your installed version and selected features](https://iqtree.github.io/iqtree2/). IQ-TREE 2 and 3 have different general citations. For ModelFinder, cite [Kalyaanamoorthy et al. (2017)](https://doi.org/10.1038/nmeth.4285). Cladistica uses standard bootstrap (`-b`), so do not describe its default analysis as UFBoot.
- **MrBayes 3.2:** [Ronquist et al. (2012)](https://doi.org/10.1093/sysbio/sys029). The [official manual](https://nbisweden.github.io/MrBayes/manual.html) lists additional citations for parallel execution, BEAGLE and other features.
- **Biopython:** [Cock et al. (2009)](https://doi.org/10.1093/bioinformatics/btp163).

Record software versions, marker boundaries, partition definitions, missing-data handling, models/priors, bootstrap replicates, MCMC settings, seeds from native logs, total generations, continuation history, burn-in, and convergence diagnostics. Retain `_work/`, `BI_resume/`, and the original input data. Describe manual sequence or specimen decisions separately from automatic processing.

## 日本語

Cladisticaのバージョンとともに、実際に使用した外部ソフトウェアと配列の原著を引用してください。アクセッションの版、取得日、標本との対応、解析条件、手動修正、継ぎ足し履歴、収束診断を記録します。ソフトウェアの実行成功だけを、解析の収束や系統関係の妥当性の根拠にしないでください。
