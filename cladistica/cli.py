from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

from . import __version__
from .align import align_by_marker
from .concat import concatenate_alignments
from .config import DEFAULT_MARKERS
from .download import download_fasta_from_accession_table
from .package import package_run_outputs
from .progress import PipelineProgress
from .trees import run_tree_analyses
from .workflow import (
    default_run_dir,
    load_env_file,
    run_accession_pipeline,
    run_accession_survey,
    run_from_fasta,
    run_from_selected_accessions,
    run_full_workflow,
)


PROJECT_DIR = Path(__file__).resolve().parents[1]


PIPELINE_DESCRIPTIONS = {
    "accessions": (
        "GenBank からターゲット Genus と outgroup の cpDNA 候補レコードを取得し、marker 領域を抽出します。"
        "CDS marker と noncoding marker を分け、CDS は partial Sanger fragment を許容しつつ internal stop/pseudo/codon frame の QC を記録します。"
        "デフォルトでは不確実な同定名を除外し、抽出済み marker coverage、QC、出版済み/受理済みらしさ、抽出配列長の順で 1 taxon = 1 sample を選びます。"
    ),
    "download": (
        "確認済み accession_selected.csv を読み、marker ごとに FASTA を取得します。"
        "complete plastome など長い accession は全長を流さず、注釈 feature から marker 領域を抽出できる場合だけ出力します。"
    ),
    "align": (
        "GenBank 由来 FASTA と手元 FASTA を marker ごとに統合し、MUSCLE で alignment を作ります。"
        "同じ sample が重複した場合は手元 FASTA を優先し、terminal gap は欠損データ `?` として扱います。"
    ),
    "concat": (
        "marker ごとの処理済み alignment を連結し、欠けた region は `?` で埋めます。"
        "partition 座標を作成し、最終的に partitions.txt と BI.nex へ反映して下流の ML/BI に渡します。"
    ),
    "tree": (
        "連結済み alignment から IQ-TREE の ModelFinder/ML bootstrap と MrBayes 用 NEXUS を実行します。"
        "MrBayes の自動診断はスクリーニングなので、論文化前には Tracer 等で ESS と trace を手動確認する前提です。"
    ),
    "package": (
        "解析ディレクトリから accession、marker FASTA、連結済み alignment、partition、BI/ML tree、"
        "Tracer 用 run1.p/run2.p、summly、log だけを単一フォルダへ集約します。"
    ),
    "survey": "GenBank の全候補を accession_all.csv にまとめ、選抜や系統解析を行わずに停止します。",
    "resume": (
        "手作業の accession_selected.csv、マーカー別 FASTA、または連結済み FASTA から解析を再開します。"
        "追加 FASTA を統合して alignment、concatenation、ML、BI を一括実行できます。"
    ),
}

USAGE_CAPTIONS = [
    (
        "Example 1: 自作FASTAからML・BIを一括解析",
        "マーカー別FASTAなら `resume --fasta-dir DIR`、すでにalignment・concatenate済みなら "
        "`resume --concatenated-fasta FILE` を使います。",
    ),
    (
        "Example 2: NCBI登録データの全容だけを確認",
        "まずは `survey --genus GENUS --outgroup ... --markers rbcL trnL-F` の2領域から始めると軽量です。"
        "accession_all.csv、summly.txt、run.logだけを返します。",
    ),
    (
        "Example 3: 手選別したaccession_selected.csvから再開",
        "`resume --accession-selected accession_selected.csv --accession-all accession_all.csv` で、"
        "配列取得からML・BIまでを一括実行します。",
    ),
    (
        "Example 4: Cladistica配列と自作配列を統合",
        "`resume --fasta-dir CLADISTICA_OUTPUT --add-fasta-dir MY_FASTA` でsample IDを照合して統合し、"
        "alignment、concatenation、ML、BIまで進めます。同一sample IDは追加側を優先します。",
    ),
]

INFERENCE_PROGRESS_STAGES = [
    ("model", "ModelFinder"),
    ("bootstrap", "Bootstrap replicates"),
    ("ml", "ML tree search"),
    ("bi", "MrBayes MCMC"),
    ("bi_summary", "BI summary and consensus"),
]


def split_markers(values: list[str] | None) -> list[str] | None:
    if not values:
        return None
    markers: list[str] = []
    for value in values:
        markers.extend(item.strip() for item in value.split(",") if item.strip())
    return markers or None


def ncbi_email(args: argparse.Namespace) -> str:
    load_env_file(PROJECT_DIR / ".env")
    value = (getattr(args, "email", "") or os.environ.get("NCBI_EMAIL", "")).strip()
    if not value and sys.stdin.isatty():
        value = input("NCBI email: ").strip()
    if "@" not in value:
        raise ValueError("A valid NCBI email is required. Use --email or NCBI_EMAIL.")
    return value


def add_marker_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--markers",
        nargs="+",
        default=None,
        help="Marker names. Default: " + ", ".join(marker.name for marker in DEFAULT_MARKERS),
    )


def add_accession_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--genus", required=True, help="Target genus, e.g. Hymenasplenium.")
    parser.add_argument("--outgroup", "--outgroups", nargs="*", default=[], help="Outgroup species names.")
    parser.add_argument("--email", default="", help="NCBI contact email. Defaults to NCBI_EMAIL.")
    parser.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY", ""), help="Optional NCBI API key.")
    parser.add_argument("--retmax", type=int, default=500, help="Maximum GenBank IDs per query name.")
    parser.add_argument("--max-samples-per-taxon", type=int, default=1, help="Default: 1, enforcing 1 taxon = 1 sample.")
    parser.add_argument("--trusted-voucher-keywords", nargs="*", default=[], help="Voucher keywords that may be kept with --include-trusted-extra.")
    parser.add_argument("--include-trusted-extra", action="store_true", help="Keep trusted voucher matches in addition to the top sample per taxon.")
    add_marker_argument(parser)


def add_inference_arguments(parser: argparse.ArgumentParser, *, include_program_paths: bool = True) -> None:
    parser.add_argument("--threads", default="AUTO", help="IQ-TREE thread setting. Default: AUTO.")
    parser.add_argument("--bootstrap", type=int, default=1000, help="IQ-TREE non-parametric bootstrap replicates.")
    parser.add_argument("--ngen", type=int, default=1_000_000, help="MrBayes generations.")
    parser.add_argument("--skip-ml", action="store_true", help="Skip IQ-TREE ML.")
    parser.add_argument("--skip-bi", action="store_true", help="Skip MrBayes BI.")
    if include_program_paths:
        parser.add_argument("--muscle", default="", help="MUSCLE command name/path.")
        parser.add_argument("--iqtree", default="", help="IQ-TREE command name/path.")
        parser.add_argument("--mrbayes", default="", help="MrBayes command name/path.")


def add_progress_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--no-progress", action="store_true", help="Disable the animated progress display.")


def progress_display(args: argparse.Namespace, stages: list[tuple[str, str]]) -> PipelineProgress:
    return PipelineProgress(stages, enabled=not getattr(args, "no_progress", False))


def cmd_describe(args: argparse.Namespace) -> int:
    print(f"Cladistica v{__version__}")
    print()
    for name, text in PIPELINE_DESCRIPTIONS.items():
        print(f"{name}: {text}")
    return 0


def cmd_examples(args: argparse.Namespace) -> int:
    print(f"Cladistica v{__version__} usage captions")
    print()
    for title, text in USAGE_CAPTIONS:
        print(title)
        print(text)
        print()
    return 0


def cmd_demo(args: argparse.Namespace) -> int:
    stages = [
        ("survey", "NCBI survey"),
        ("download", "Sequence download"),
        ("align", "MUSCLE marker alignments"),
        ("model", "ModelFinder"),
        ("bootstrap", "Bootstrap replicates"),
        ("ml", "ML tree search"),
        ("bi", "MrBayes BI"),
        ("bi_summary", "BI summary and consensus"),
    ]
    pause = max(args.duration / (len(stages) * 5 + 1), 0.08)
    progress = PipelineProgress(stages)
    progress.feed(
        "Hymenasplenium hondoense",
        "Hymenasplenium obliquissimum",
        "Asplenium setoi",
        "rbcL",
        "matK",
        "NC_035840.1",
    )
    progress.feed_sequence(
        "H. hondoense | rbcL",
        "ATGTCACCACAAACAGAGACTAAAGCAAGTGTTGGATTCAAAGCTGGTGTTAAAGATTATAAATTG",
    )
    with progress:
        time.sleep(pause)
        for key, _ in stages:
            progress.start(key)
            for percent in (0, 25, 50, 75):
                if key == "bootstrap":
                    detail = f"{percent * 10}/1000 replicates"
                elif key == "bi":
                    detail = f"{percent * 10_000:,}/1,000,000 generations"
                elif key == "model":
                    detail = f"{max(1, percent // 5)} candidates"
                else:
                    detail = ""
                progress.set_progress(
                    key,
                    percent,
                    detail,
                    approximate=key in {"model", "ml", "bi_summary"},
                )
                time.sleep(pause)
                if args.fail and key == "align" and percent == 50:
                    progress.fail(key, "process stopped")
                    time.sleep(max(pause * 2, 1.0))
                    break
            if args.fail and key == "align":
                break
            progress.succeed(key, "done")
            time.sleep(pause)
    return 0


def cmd_survey(args: argparse.Namespace) -> int:
    with progress_display(
        args,
        [("accessions", "NCBI survey and marker extraction"), ("package", "Final output")],
    ) as progress:
        result = run_accession_survey(
            project_dir=PROJECT_DIR,
            genus=args.genus,
            outgroups=args.outgroup,
            output_dir=args.output.expanduser().resolve() if args.output else None,
            email=ncbi_email(args),
            api_key=args.api_key,
            markers=split_markers(args.markers),
            retmax=args.retmax,
            progress=progress,
        )
    print(f"Output: {result.output_dir}")
    print(f"Accession rows: {result.records_written}")
    return 0


def cmd_resume(args: argparse.Namespace) -> int:
    output_dir = args.output.expanduser().resolve() if args.output else None
    markers = split_markers(args.markers)
    additional = args.add_fasta_dir.expanduser().resolve() if args.add_fasta_dir else None
    common = {
        "project_dir": PROJECT_DIR,
        "output_dir": output_dir,
        "markers": markers,
        "muscle_command": args.muscle,
        "iqtree_command": args.iqtree,
        "mrbayes_command": args.mrbayes,
        "threads": args.threads,
        "bootstrap": args.bootstrap,
        "ngen": args.ngen,
        "skip_ml": args.skip_ml,
        "skip_bi": args.skip_bi,
    }
    if args.accession_selected:
        stages = [
            ("download", "Selected GenBank sequences"),
            ("align", "MUSCLE marker alignments"),
            ("concat", "Concatenated alignment"),
            *INFERENCE_PROGRESS_STAGES,
            ("package", "Final output"),
        ]
    elif args.fasta_dir:
        stages = [
            ("align", "MUSCLE marker alignments"),
            ("concat", "Concatenated alignment"),
            *INFERENCE_PROGRESS_STAGES,
            ("package", "Final output"),
        ]
    else:
        stages = [
            *INFERENCE_PROGRESS_STAGES,
            ("package", "Final output"),
        ]
    with progress_display(args, stages) as progress:
        common["progress"] = progress
        if args.accession_selected:
            if args.partition_file:
                raise ValueError("--partition-file is only valid with --concatenated-fasta.")
            result = run_from_selected_accessions(
                accession_selected=args.accession_selected.expanduser().resolve(),
                accession_all=args.accession_all.expanduser().resolve() if args.accession_all else None,
                additional_fasta_dir=additional,
                email=ncbi_email(args),
                api_key=args.api_key,
                **common,
            )
        else:
            if args.accession_all:
                raise ValueError("--accession-all is only valid with --accession-selected.")
            result = run_from_fasta(
                fasta_dir=args.fasta_dir.expanduser().resolve() if args.fasta_dir else None,
                additional_fasta_dir=additional,
                concatenated_fasta=args.concatenated_fasta.expanduser().resolve() if args.concatenated_fasta else None,
                partition_file=args.partition_file.expanduser().resolve() if args.partition_file else None,
                **common,
            )
    print(f"Output: {result.output_dir}")
    print(f"Samples analyzed: {result.records_written}")
    if result.warnings:
        print(f"Warnings: {len(result.warnings)}", file=sys.stderr)
        for warning in result.warnings:
            print(f"- {warning}", file=sys.stderr)
    return 0


def cmd_accessions(args: argparse.Namespace) -> int:
    markers = split_markers(args.markers)
    output_dir = args.output.expanduser().resolve()
    with progress_display(args, [("accessions", "NCBI records and representative selection")]) as progress:
        result = run_accession_pipeline(
            project_dir=PROJECT_DIR,
            genus=args.genus,
            outgroups=args.outgroup,
            output_dir=output_dir,
            email=ncbi_email(args),
            api_key=args.api_key,
            markers=markers,
            retmax=args.retmax,
            max_samples_per_taxon=args.max_samples_per_taxon,
            trusted_voucher_keywords=args.trusted_voucher_keywords,
            include_trusted_extra=args.include_trusted_extra,
            progress=progress,
        )
    print(f"Output: {result.output_dir}")
    print(f"All accessions: {output_dir / 'accession_all.csv'}")
    print(f"Selected accessions: {output_dir / 'accession_selected.csv'}")
    print(f"Selected FASTA: {output_dir / 'fasta_by_marker'}")
    print(f"Selected samples: {result.records_written}")
    return 0


def cmd_download(args: argparse.Namespace) -> int:
    with progress_display(args, [("download", "Selected GenBank sequences")]) as progress:
        result = download_fasta_from_accession_table(
            accession_table=args.accession_table.expanduser().resolve(),
            output_dir=args.output.expanduser().resolve(),
            email=ncbi_email(args),
            api_key=args.api_key,
            markers=split_markers(args.markers),
            progress=progress,
        )
    print(f"Output: {result.output_dir}")
    print(f"Sequences written: {result.records_written}")
    if result.warnings:
        print(f"Warnings: {len(result.warnings)}", file=sys.stderr)
        return 3
    return 0


def cmd_align(args: argparse.Namespace) -> int:
    with progress_display(args, [("align", "MUSCLE marker alignments")]) as progress:
        result = align_by_marker(
            genbank_dir=args.genbank_dir.expanduser().resolve(),
            user_dir=args.user_dir.expanduser().resolve(),
            output_dir=args.output.expanduser().resolve(),
            markers=split_markers(args.markers),
            muscle_command=args.muscle,
            dry_run=args.dry_run,
            progress=progress,
        )
    print(f"Output: {result.output_dir}")
    print(f"Marker alignments written: {result.records_written}")
    return 0


def cmd_concat(args: argparse.Namespace) -> int:
    with progress_display(args, [("concat", "Concatenated alignment")]) as progress:
        result = concatenate_alignments(
            input_dir=args.input.expanduser().resolve(),
            output_dir=args.output.expanduser().resolve(),
            markers=split_markers(args.markers),
            strict=args.strict,
            progress=progress,
        )
    print(f"Output: {result.output_dir}")
    print(f"Samples written: {result.records_written}")
    return 0


def cmd_tree(args: argparse.Namespace) -> int:
    with progress_display(
        args,
        INFERENCE_PROGRESS_STAGES,
    ) as progress:
        result = run_tree_analyses(
            input_dir=args.input.expanduser().resolve(),
            output_dir=args.output.expanduser().resolve(),
            threads=args.threads,
            bootstrap=args.bootstrap,
            skip_ml=args.skip_ml,
            skip_bi=args.skip_bi,
            iqtree_command=args.iqtree,
            mrbayes_command=args.mrbayes,
            ngen=args.ngen,
            nruns=args.nruns,
            nchains=args.nchains,
            burninfrac=args.burninfrac,
            dry_run=args.dry_run,
            progress=progress,
        )
    print(f"Output: {result.output_dir}")
    print(f"Analyses completed: {result.records_written}")
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    stages = [
        ("accessions", "NCBI records and representative selection"),
        ("align", "MUSCLE marker alignments"),
        ("concat", "Concatenated alignment"),
        *INFERENCE_PROGRESS_STAGES,
        ("package", "Final output"),
    ]
    with progress_display(args, stages) as progress:
        result = run_full_workflow(
            project_dir=PROJECT_DIR,
            genus=args.genus,
            outgroups=args.outgroup,
            output_dir=args.output.expanduser().resolve() if args.output else None,
            email=ncbi_email(args),
            api_key=args.api_key,
            markers=split_markers(args.markers),
            retmax=args.retmax,
            max_samples_per_taxon=args.max_samples_per_taxon,
            trusted_voucher_keywords=args.trusted_voucher_keywords,
            include_trusted_extra=args.include_trusted_extra,
            my_fasta_dir=args.my_fasta_dir.expanduser().resolve() if args.my_fasta_dir else None,
            skip_ml=args.skip_ml,
            skip_bi=args.skip_bi,
            muscle_command=args.muscle,
            iqtree_command=args.iqtree,
            mrbayes_command=args.mrbayes,
            threads=args.threads,
            bootstrap=args.bootstrap,
            ngen=args.ngen,
            progress=progress,
        )
    print(f"Output: {result.output_dir}")
    print(f"Selected samples: {result.records_written}")
    return 0


def cmd_package(args: argparse.Namespace) -> int:
    archive_path = args.archive.expanduser().resolve() if args.archive else None
    output_dir = args.output.expanduser().resolve() if args.output else default_run_dir(PROJECT_DIR)
    result = package_run_outputs(
        run_dir=args.run_dir.expanduser().resolve(),
        output_dir=output_dir,
        archive_path=archive_path,
    )
    print(f"Output: {result.output_dir}")
    print(f"Files written: {result.records_written}")
    if archive_path:
        print(f"Archive: {archive_path}")
    if result.warnings:
        print(f"Warnings: {len(result.warnings)}", file=sys.stderr)
        for warning in result.warnings:
            print(f"- {warning}", file=sys.stderr)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="cladistica", description="Cladistica: cpDNA accession-to-tree workflows.")
    parser.add_argument("--version", action="version", version=f"Cladistica {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    describe = subparsers.add_parser("describe", help="Explain the default behavior of each pipeline.")
    describe.set_defaults(func=cmd_describe)

    examples = subparsers.add_parser("examples", help="Show four recommended usage examples.")
    examples.set_defaults(func=cmd_examples)

    demo = subparsers.add_parser("demo", help="Preview the sequence-stream terminal animation.")
    demo.add_argument(
        "--duration",
        type=float,
        default=15.0,
        help="Approximate demo duration in seconds. Default: 15.",
    )
    demo.add_argument(
        "--fail",
        action="store_true",
        help="Demonstrate character collapse after a detected analysis failure.",
    )
    demo.set_defaults(func=cmd_demo)

    survey = subparsers.add_parser("survey", help="Stop after writing accession_all.csv.")
    survey.add_argument("--genus", required=True, help="Target genus, e.g. Hymenasplenium.")
    survey.add_argument("--outgroup", "--outgroups", nargs="*", default=[], help="Outgroup species names.")
    survey.add_argument("--email", default="", help="NCBI contact email. Defaults to NCBI_EMAIL.")
    survey.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY", ""), help="Optional NCBI API key.")
    survey.add_argument(
        "--retmax",
        type=int,
        default=0,
        help="Maximum unique GenBank IDs per query name. Default: 0 (no limit).",
    )
    survey.add_argument("--output", type=Path, help="Output directory. Default: output/YYYYMMDD_N.")
    add_marker_argument(survey)
    add_progress_argument(survey)
    survey.set_defaults(func=cmd_survey)

    resume = subparsers.add_parser("resume", help="Resume from a selected accession table or FASTA.")
    source = resume.add_mutually_exclusive_group(required=True)
    source.add_argument("--accession-selected", type=Path, help="Hand-edited accession_selected.csv.")
    source.add_argument("--fasta-dir", type=Path, help="Directory containing unaligned marker FASTA files.")
    source.add_argument("--concatenated-fasta", type=Path, help="Aligned concatenated FASTA.")
    resume.add_argument("--accession-all", type=Path, help="Optional accession_all.csv retained with selected accessions.")
    resume.add_argument("--add-fasta-dir", type=Path, help="Additional marker FASTA; duplicate sample IDs override the primary input.")
    resume.add_argument("--partition-file", type=Path, help="Optional partition CSV/TSV, IQ-TREE text, or NEXUS charset file.")
    resume.add_argument("--output", type=Path, help="Output directory. Default: output/YYYYMMDD_N.")
    resume.add_argument("--email", default="", help="NCBI contact email; required only with --accession-selected.")
    resume.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY", ""), help="Optional NCBI API key.")
    add_marker_argument(resume)
    add_inference_arguments(resume)
    add_progress_argument(resume)
    resume.set_defaults(func=cmd_resume)

    accessions = subparsers.add_parser("accessions", help="Build a selected cpDNA accession table from GenBank.")
    add_accession_arguments(accessions)
    accessions.add_argument("--output", type=Path, default=Path("output/accessions"), help="Output directory for accession CSV and marker FASTA.")
    add_progress_argument(accessions)
    accessions.set_defaults(func=cmd_accessions)

    download = subparsers.add_parser("download", help="Download marker FASTA files from an accession table.")
    download.add_argument("--accession-table", type=Path, required=True, help="Path to accession_selected.csv.")
    download.add_argument("--output", type=Path, default=Path("output/fasta"), help="Output directory.")
    download.add_argument("--email", default="", help="NCBI contact email. Defaults to NCBI_EMAIL.")
    download.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY", ""), help="Optional NCBI API key.")
    add_marker_argument(download)
    add_progress_argument(download)
    download.set_defaults(func=cmd_download)

    align = subparsers.add_parser("align", help="Align marker FASTA files with MUSCLE.")
    align.add_argument("--genbank-dir", type=Path, required=True, help="Directory containing GenBank FASTA files by marker.")
    align.add_argument("--user-dir", type=Path, default=Path("input/my_fasta_by_marker"), help="Directory containing user FASTA files by marker.")
    align.add_argument("--output", type=Path, default=Path("output/alignments"), help="Output directory.")
    align.add_argument("--muscle", default="", help="MUSCLE command name/path.")
    align.add_argument("--dry-run", action="store_true", help="Validate inputs without running MUSCLE.")
    add_marker_argument(align)
    add_progress_argument(align)
    align.set_defaults(func=cmd_align)

    concat = subparsers.add_parser("concat", help="Concatenate aligned marker FASTA files.")
    concat.add_argument("--input", type=Path, required=True, help="Directory containing aligned marker FASTA files.")
    concat.add_argument("--output", type=Path, default=Path("output/concatenated"), help="Output directory.")
    concat.add_argument("--strict", action="store_true", help="Fail when any requested marker is missing.")
    add_marker_argument(concat)
    add_progress_argument(concat)
    concat.set_defaults(func=cmd_concat)

    tree = subparsers.add_parser("tree", help="Run IQ-TREE and MrBayes from a concatenated alignment.")
    tree.add_argument("--input", type=Path, required=True, help="Directory containing concatenated FASTA and optional partition_coordinates.tsv.")
    tree.add_argument("--output", type=Path, default=Path("output/trees"), help="Output directory.")
    tree.add_argument("--threads", default="AUTO", help="IQ-TREE thread setting. Default: AUTO.")
    tree.add_argument("--bootstrap", type=int, default=1000, help="IQ-TREE standard non-parametric bootstrap replicates.")
    tree.add_argument("--skip-ml", action="store_true", help="Skip IQ-TREE ML.")
    tree.add_argument("--skip-bi", action="store_true", help="Skip MrBayes BI.")
    tree.add_argument("--iqtree", default="", help="IQ-TREE command name/path.")
    tree.add_argument("--mrbayes", default="", help="MrBayes command name/path.")
    tree.add_argument("--ngen", type=int, default=1_000_000, help="MrBayes generations.")
    tree.add_argument("--nruns", type=int, default=2, help="MrBayes nruns.")
    tree.add_argument("--nchains", type=int, default=4, help="MrBayes nchains.")
    tree.add_argument("--burninfrac", type=float, default=0.25, help="MrBayes burn-in fraction.")
    tree.add_argument("--dry-run", action="store_true", help="Write scripts/reports without running external tools.")
    add_progress_argument(tree)
    tree.set_defaults(func=cmd_tree)

    run = subparsers.add_parser("run", help="Run accession, FASTA, alignment, concatenation, and tree stages.")
    add_accession_arguments(run)
    run.add_argument("--output", type=Path, help="Run output directory. Default: output/YYYYMMDD_N.")
    run.add_argument("--my-fasta-dir", type=Path, help="Optional user FASTA directory by marker.")
    add_inference_arguments(run)
    add_progress_argument(run)
    run.set_defaults(func=cmd_run)

    package = subparsers.add_parser("package", help="Collect required results in one flat output directory.")
    package.add_argument("--run-dir", type=Path, required=True, help="Existing Cladistica run directory.")
    package.add_argument("--output", type=Path, help="Output directory. Default: output/YYYYMMDD_N.")
    package.add_argument("--archive", type=Path, help="Optional ZIP path.")
    package.set_defaults(func=cmd_package)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except (KeyboardInterrupt, EOFError):
        print("\nInterrupted.", file=sys.stderr)
        return 130
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
