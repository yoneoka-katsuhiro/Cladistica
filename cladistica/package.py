from __future__ import annotations

import shutil
import zipfile
from pathlib import Path

from . import __version__
from .io import read_table, write_csv
from .models import PipelineResult


def reset_directory(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)
    for child in path.iterdir():
        if child.is_dir():
            shutil.rmtree(child)
        else:
            child.unlink()


def copy_file(source: Path, destination: Path, warnings: list[str]) -> bool:
    if not source.exists():
        warnings.append(f"missing: {source}")
        return False
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    return True


def copy_glob_flat(source_dir: Path, pattern: str, output_dir: Path) -> int:
    if not source_dir.exists():
        return 0
    copied = 0
    output_dir.mkdir(parents=True, exist_ok=True)
    for source in sorted(source_dir.glob(pattern)):
        if source.is_file():
            shutil.copy2(source, output_dir / source.name)
            copied += 1
    return copied


def write_zip(source_dir: Path, archive_path: Path) -> None:
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(source_dir.rglob("*")):
            if path.is_file():
                archive.write(path, path.relative_to(source_dir))


def append_file_section(lines: list[str], title: str, path: Path) -> bool:
    if not path.exists():
        return False
    lines.extend(["", f"## {title}", ""])
    lines.extend(path.read_text(encoding="utf-8", errors="replace").splitlines())
    return True


def first_existing(paths: list[Path]) -> Path | None:
    return next((path for path in paths if path.exists()), None)


def copy_table_as_csv(source: Path, destination: Path, warnings: list[str]) -> bool:
    if source.suffix.lower() == ".csv":
        return copy_file(source, destination, warnings)
    try:
        rows, fieldnames = read_table(source)
        write_csv(destination, rows, fieldnames)
        return True
    except Exception as exc:
        warnings.append(f"could not convert {source} to CSV: {exc}")
        return False


def write_combined_log(work_dir: Path, output_dir: Path) -> None:
    lines: list[str] = ["Cladistica run log", ""]
    sections = [
        ("Accession selection", [work_dir / "accessions" / "run.log"]),
        ("GenBank query", [work_dir / "accessions" / "genbank_query_log.tsv", work_dir / "01_accessions" / "genbank_query_log.tsv"]),
        ("Rejected GenBank records", [work_dir / "01_accessions" / "rejected_records.tsv"]),
        ("Alignment QC", [work_dir / "alignments" / "qc" / "all_markers_qc.tsv", work_dir / "03_alignments" / "qc" / "all_markers_qc.tsv"]),
        ("Alignment rejected records", [work_dir / "alignments" / "qc" / "rejected_before_alignment.tsv", work_dir / "03_alignments" / "qc" / "rejected_before_alignment.tsv"]),
        ("Alignment duplicate records", [work_dir / "alignments" / "qc" / "duplicate_records.tsv", work_dir / "03_alignments" / "qc" / "duplicate_records.tsv"]),
        ("IQ-TREE", [work_dir / "trees" / "run_iqtree_ml.log", work_dir / "05_trees" / "run_iqtree_ml.log"]),
        ("MrBayes", [work_dir / "trees" / "run_mrbayes.log", work_dir / "05_trees" / "run_mrbayes.log"]),
    ]
    wrote = False
    for title, paths in sections:
        path = first_existing(paths)
        if path:
            wrote = append_file_section(lines, title, path) or wrote
    if not wrote:
        lines.append("No log files were produced.")
    (output_dir / "run.log").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_flat_summary(work_dir: Path, output_dir: Path) -> None:
    lines = [f"Cladistica v{__version__} analysis summary", ""]
    for title, paths in [
        ("Workflow", [work_dir / "workflow_summly.txt", work_dir / "workflow_summary.txt"]),
        ("Accession selection", [work_dir / "accessions" / "summly.txt", work_dir / "01_accessions" / "summary.txt"]),
        ("Alignment", [work_dir / "alignments" / "alignment_report.txt", work_dir / "03_alignments" / "alignment_report.txt"]),
        ("Concatenation", [work_dir / "concatenated" / "analysis_report.txt", work_dir / "04_concatenated" / "analysis_report.txt"]),
        ("Tree analysis", [work_dir / "trees" / "05_reports" / "analysis_summary.txt", work_dir / "05_trees" / "05_reports" / "analysis_summary.txt"]),
        (
            "Model selection",
            [
                work_dir / "trees" / "05_reports" / "model_selection_summary.txt",
                work_dir / "05_trees" / "05_reports" / "model_selection_summary.txt",
            ],
        ),
    ]:
        path = first_existing(paths)
        if path:
            append_file_section(lines, title, path)
    key_files = {
        "accession_all.csv": "every retrieved candidate or rejected GenBank record.",
        "accession_selected.csv": "samples selected for analysis.",
        "concatenated.fasta": "final aligned and concatenated matrix.",
        "partitions.txt": "IQ-TREE partition definitions.",
        "BI.nex": "executable MrBayes analysis file, including charset definitions.",
        "ML.tre": "IQ-TREE maximum-likelihood tree.",
        "BI.tre": "MrBayes consensus tree.",
        "run1.p": "MrBayes run 1 parameter trace for Tracer.",
        "run2.p": "MrBayes run 2 parameter trace for Tracer.",
        "run.log": "combined query, QC, IQ-TREE, and MrBayes log.",
    }
    lines.extend(["", "## Key files", ""])
    marker_fastas = sorted(
        path.name
        for path in output_dir.glob("*.fasta")
        if path.name != "concatenated.fasta"
    )
    if marker_fastas:
        lines.append(f"{', '.join(marker_fastas)}: combined unaligned sequences by marker.")
    for name, description in key_files.items():
        if (output_dir / name).exists() or name == "run.log":
            lines.append(f"{name}: {description}")
    (output_dir / "summly.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def package_run_outputs(
    *,
    run_dir: Path,
    output_dir: Path,
    archive_path: Path | None = None,
    reset_output: bool = True,
    expect_accessions: bool = True,
    expect_marker_fasta: bool = True,
    expect_concatenated: bool = True,
    expect_ml: bool = True,
    expect_bi: bool = True,
) -> PipelineResult:
    """Create one flat, reader-facing result directory from temporary stage outputs."""
    warnings: list[str] = []
    if reset_output:
        reset_directory(output_dir)
    else:
        output_dir.mkdir(parents=True, exist_ok=True)

    accession_dir = first_existing([run_dir / "accessions", run_dir / "01_accessions"])
    alignment_dir = first_existing([run_dir / "alignments", run_dir / "03_alignments"])
    concat_dir = first_existing([run_dir / "concatenated", run_dir / "04_concatenated"])
    tree_dir = first_existing([run_dir / "trees", run_dir / "05_trees"])

    if accession_dir:
        all_accessions = first_existing([accession_dir / "accession_all.csv", accession_dir / "all_candidate_accessions.tsv"])
        selected_accessions = first_existing([accession_dir / "accession_selected.csv", accession_dir / "accession_table.tsv"])
        if all_accessions:
            copy_table_as_csv(all_accessions, output_dir / "accession_all.csv", warnings)
        elif expect_accessions:
            warnings.append("missing: accession_all.csv")
        if selected_accessions:
            copy_table_as_csv(selected_accessions, output_dir / "accession_selected.csv", warnings)
        elif expect_accessions:
            warnings.append("missing: accession_selected.csv")
    elif expect_accessions:
        warnings.append("missing: accession stage output")

    fasta_candidates: list[Path] = []
    if alignment_dir:
        fasta_candidates.append(alignment_dir / "intermediate" / "combined_by_marker")
    if accession_dir:
        fasta_candidates.extend(
            [
                accession_dir / "fasta_by_marker",
                accession_dir / "selected_fasta_by_marker",
            ]
        )
    fasta_source = first_existing(fasta_candidates)

    fasta_count = copy_glob_flat(fasta_source, "*.fasta", output_dir) if fasta_source else 0
    if fasta_count == 0 and expect_marker_fasta:
        warnings.append("missing: marker FASTA files")

    if concat_dir:
        copy_file(concat_dir / "concatenated_cpDNA.fasta", output_dir / "concatenated.fasta", warnings)
    elif expect_concatenated:
        warnings.append("missing: concatenated stage output")

    if tree_dir:
        required_tree_files = {"partitions.txt"}
        if expect_ml:
            required_tree_files.add("ML.tre")
        if expect_bi:
            required_tree_files.update({"BI.nex", "BI.tre", "run1.p", "run2.p"})
        tree_sources = {
            "partitions.txt": first_existing(
                [tree_dir / "00_inputs" / "iqtree_partitions.txt", concat_dir / "partition_coordinates.nexus_block.txt"] if concat_dir
                else [tree_dir / "00_inputs" / "iqtree_partitions.txt"]
            ),
            "BI.nex": first_existing(
                [tree_dir / "03_mrbayes" / "mrbayes_analysis.nex", tree_dir / "00_inputs" / "mrbayes_analysis.nex"]
            ),
            "ML.tre": first_existing(
                [
                    tree_dir / "02_iqtree_ml" / "cladistica_ml.treefile",
                    tree_dir / "04_trees" / "cladistica_ml.treefile",
                    tree_dir / "02_iqtree_ml" / "phylostudio_ml.treefile",
                    tree_dir / "04_trees" / "phylostudio_ml.treefile",
                ]
            ),
            "BI.tre": first_existing(
                [
                    tree_dir / "03_mrbayes" / "mrbayes_analysis.nex.con.tre",
                    tree_dir / "04_trees" / "mrbayes_analysis.nex.con.tre",
                ]
            ),
            "run1.p": tree_dir / "03_mrbayes" / "mrbayes_analysis.nex.run1.p",
            "run2.p": tree_dir / "03_mrbayes" / "mrbayes_analysis.nex.run2.p",
        }
        for destination_name, source in tree_sources.items():
            if source:
                copy_file(source, output_dir / destination_name, warnings)
            elif destination_name in required_tree_files and (expect_ml or expect_bi):
                warnings.append(f"missing: {destination_name}")
    elif expect_ml or expect_bi:
        warnings.append("missing: tree stage output")

    write_flat_summary(run_dir, output_dir)
    write_combined_log(run_dir, output_dir)

    if archive_path:
        write_zip(output_dir, archive_path)

    files_written = sum(1 for path in output_dir.iterdir() if path.is_file())
    return PipelineResult(str(output_dir), records_written=files_written, warnings=warnings)
