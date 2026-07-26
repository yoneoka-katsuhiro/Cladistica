from __future__ import annotations

import os
import re
import shutil
import tempfile
from datetime import date
from pathlib import Path

from .accessions import build_accession_table_from_genbank_records
from .align import align_by_marker
from .concat import concatenate_alignments
from .config import marker_map
from .download import download_fasta_from_accession_table
from .io import read_fasta, read_table, write_tsv
from .models import PipelineResult
from .ncbi import query_genbank_seqrecords
from .package import package_run_outputs, reset_directory, write_combined_log
from .progress import PipelineProgress
from .trees import run_tree_analyses


def load_env_file(path: Path) -> None:
    if not path.exists():
        return
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", key):
            continue
        value = value.strip()
        if len(value) >= 2 and value[0] == value[-1] and value[0] in {"'", '"'}:
            value = value[1:-1]
        os.environ.setdefault(key, value)


def default_run_dir(root: Path, genus: str = "") -> Path:
    output_root = root / "output"
    run_date = date.today().strftime("%Y%m%d")
    index = 1
    while True:
        candidate = output_root / f"{run_date}_{index}"
        if not candidate.exists():
            return candidate
        index += 1


def require_files(paths: list[Path]) -> None:
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise RuntimeError("Required analysis output was not created: " + ", ".join(missing))


def inference_outputs(tree_dir: Path, *, skip_ml: bool, skip_bi: bool) -> list[Path]:
    required: list[Path] = []
    if not skip_ml:
        required.append(tree_dir / "02_iqtree_ml" / "cladistica_ml.treefile")
    if not skip_bi:
        required.extend(
            [
                tree_dir / "03_mrbayes" / "mrbayes_analysis.nex",
                tree_dir / "03_mrbayes" / "mrbayes_analysis.nex.con.tre",
                tree_dir / "03_mrbayes" / "mrbayes_analysis.nex.run1.p",
                tree_dir / "03_mrbayes" / "mrbayes_analysis.nex.run2.p",
            ]
        )
    return required


def run_marker_fasta_stages(
    *,
    primary_fasta_dir: Path,
    additional_fasta_dir: Path,
    work_dir: Path,
    markers: list[str],
    muscle_command: str,
    iqtree_command: str,
    mrbayes_command: str,
    threads: str,
    bootstrap: int,
    ngen: int,
    skip_ml: bool,
    skip_bi: bool,
    progress: PipelineProgress | None = None,
) -> PipelineResult:
    align_dir = work_dir / "alignments"
    concat_dir = work_dir / "concatenated"
    tree_dir = work_dir / "trees"
    align_by_marker(
        genbank_dir=primary_fasta_dir,
        user_dir=additional_fasta_dir,
        output_dir=align_dir,
        markers=markers,
        muscle_command=muscle_command,
        progress=progress,
    )
    concatenated = concatenate_alignments(
        input_dir=align_dir / "aligned_by_marker",
        output_dir=concat_dir,
        markers=markers,
        progress=progress,
    )
    run_tree_analyses(
        input_dir=concat_dir,
        output_dir=tree_dir,
        threads=threads,
        bootstrap=bootstrap,
        skip_ml=skip_ml,
        skip_bi=skip_bi,
        iqtree_command=iqtree_command,
        mrbayes_command=mrbayes_command,
        ngen=ngen,
        nruns=2,
        progress=progress,
    )
    require_files(
        [
            concat_dir / "concatenated_cpDNA.fasta",
            *inference_outputs(tree_dir, skip_ml=skip_ml, skip_bi=skip_bi),
        ]
    )
    return concatenated


def parse_partition_file(source: Path, alignment_length: int) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    if source.suffix.lower() in {".csv", ".tsv"}:
        table_rows, fieldnames = read_table(source)
        required = {"marker", "start", "end"}
        if not required.issubset(fieldnames):
            raise ValueError("Partition CSV/TSV needs marker, start, and end columns.")
        for row in table_rows:
            rows.append({"marker": row["marker"], "start": row["start"], "end": row["end"]})
    else:
        pattern = re.compile(
            r"^\s*(?:charset\s+|DNA\s*,\s*)?([^=,]+?)\s*=\s*(\d+)\s*-\s*(\d+)\s*;?\s*$",
            re.IGNORECASE,
        )
        for raw in source.read_text(encoding="utf-8-sig").splitlines():
            line = raw.strip()
            lowered = line.lower()
            if not line or lowered.startswith(("#", "[", "begin ", "end;")):
                continue
            match = pattern.match(line)
            if match:
                rows.append({"marker": match.group(1).strip(), "start": match.group(2), "end": match.group(3)})
    if not rows:
        raise ValueError(f"No supported partition definitions were found in {source}")

    normalized: list[dict[str, object]] = []
    seen_markers: set[str] = set()
    for row in rows:
        marker = str(row["marker"]).strip()
        start = int(str(row["start"]))
        end = int(str(row["end"]))
        if not marker or marker in seen_markers:
            raise ValueError(f"Partition names must be non-empty and unique: {marker!r}")
        if start < 1 or end < start or end > alignment_length:
            raise ValueError(f"Invalid partition coordinates for {marker}: {start}-{end}")
        seen_markers.add(marker)
        normalized.append({"marker": marker, "start": start, "end": end, "length": end - start + 1})
    for index, left in enumerate(normalized):
        for right in normalized[index + 1 :]:
            if int(left["start"]) <= int(right["end"]) and int(right["start"]) <= int(left["end"]):
                raise ValueError(
                    f"Overlapping partitions: {left['marker']} and {right['marker']}"
                )
    return normalized


def run_accession_survey(
    *,
    project_dir: Path,
    genus: str,
    outgroups: list[str],
    output_dir: Path | None,
    email: str,
    api_key: str = "",
    markers: list[str] | None = None,
    retmax: int = 500,
    progress: PipelineProgress | None = None,
) -> PipelineResult:
    run_dir = output_dir or default_run_dir(project_dir)
    marker_names = list(marker_map(markers))
    with tempfile.TemporaryDirectory(prefix="cladistica_") as temporary:
        work_dir = Path(temporary)
        accessions_dir = work_dir / "accessions"
        run_accession_pipeline(
            project_dir=project_dir,
            genus=genus,
            outgroups=outgroups,
            output_dir=accessions_dir,
            email=email,
            api_key=api_key,
            markers=marker_names,
            retmax=retmax,
            progress=progress,
        )
        rows, _ = read_table(accessions_dir / "accession_all.csv")
        if progress:
            progress.succeed("accessions", f"{len(rows)} accession rows reviewed")
        reset_directory(run_dir)
        shutil.copy2(accessions_dir / "accession_all.csv", run_dir / "accession_all.csv")
        survey_summary = [
            "Cladistica accession survey",
            "",
            f"Genus: {genus}",
            f"Outgroups: {'; '.join(outgroups) if outgroups else '(none)'}",
            f"Markers: {', '.join(marker_names)}",
            f"Rows in accession_all.csv: {len(rows)}",
            "",
            "selection_status is an automatic recommendation only.",
            "Edit or copy chosen rows into accession_selected.csv before resuming analysis.",
        ]
        (run_dir / "summly.txt").write_text("\n".join(survey_summary) + "\n", encoding="utf-8")
        if progress:
            progress.start("package", "Writing accession_all.csv, summly.txt, and run.log")
        write_combined_log(work_dir, run_dir)
        if progress:
            progress.succeed("package", f"{len(rows)} accession rows")
    return PipelineResult(str(run_dir), records_written=len(rows))


def run_accession_pipeline(
    *,
    project_dir: Path,
    genus: str,
    outgroups: list[str],
    output_dir: Path,
    email: str,
    api_key: str = "",
    markers: list[str] | None = None,
    retmax: int = 500,
    max_samples_per_taxon: int = 1,
    trusted_voucher_keywords: list[str] | None = None,
    include_trusted_extra: bool = False,
    progress: PipelineProgress | None = None,
) -> PipelineResult:
    load_env_file(project_dir / ".env")
    marker_names = list(marker_map(markers))
    if progress:
        progress.start("accessions")
        progress.set_progress("accessions", 5, "NCBI Entrez query", approximate=True)
    records = query_genbank_seqrecords(
        genus=genus,
        outgroups=outgroups,
        email=email,
        api_key=api_key,
        markers=marker_names,
        retmax=retmax,
        log_path=output_dir / "genbank_query_log.tsv",
    )
    if progress:
        progress.set_progress(
            "accessions",
            45,
            f"{len(records)} GenBank records",
            approximate=True,
        )
        progress.set_progress(
            "accessions",
            55,
            "marker extraction and ranking",
            approximate=True,
        )
    result = build_accession_table_from_genbank_records(
        records=records,
        output_dir=output_dir,
        markers=marker_names,
        max_samples_per_taxon=max_samples_per_taxon,
        trusted_voucher_keywords=trusted_voucher_keywords,
        include_trusted_extra=include_trusted_extra,
        parameters={
            "genus": genus,
            "outgroups": ";".join(outgroups),
            "retmax": retmax,
            "selection_priority": "extracted_marker_coverage > clean_or_coding_qc_passing_markers > publication_evidence > extracted_total_length",
        },
    )
    if progress:
        progress.succeed("accessions", f"{result.records_written} selected samples")
    return result


def run_full_workflow(
    *,
    project_dir: Path,
    genus: str,
    outgroups: list[str],
    output_dir: Path | None,
    email: str,
    api_key: str = "",
    markers: list[str] | None = None,
    retmax: int = 500,
    max_samples_per_taxon: int = 1,
    trusted_voucher_keywords: list[str] | None = None,
    include_trusted_extra: bool = False,
    my_fasta_dir: Path | None = None,
    skip_ml: bool = False,
    skip_bi: bool = False,
    muscle_command: str = "",
    iqtree_command: str = "",
    mrbayes_command: str = "",
    threads: str = "AUTO",
    bootstrap: int = 1000,
    ngen: int = 1_000_000,
    progress: PipelineProgress | None = None,
) -> PipelineResult:
    run_dir = output_dir or default_run_dir(project_dir, genus)
    marker_names = list(marker_map(markers))
    with tempfile.TemporaryDirectory(prefix="cladistica_") as temporary:
        work_dir = Path(temporary)
        accessions_dir = work_dir / "accessions"

        accession_result = run_accession_pipeline(
            project_dir=project_dir,
            genus=genus,
            outgroups=outgroups,
            output_dir=accessions_dir,
            email=email,
            api_key=api_key,
            markers=marker_names,
            retmax=retmax,
            max_samples_per_taxon=max_samples_per_taxon,
            trusted_voucher_keywords=trusted_voucher_keywords,
            include_trusted_extra=include_trusted_extra,
            progress=progress,
        )
        user_dir = my_fasta_dir or (project_dir / "input" / "my_fasta_by_marker")
        run_marker_fasta_stages(
            primary_fasta_dir=accessions_dir / "fasta_by_marker",
            additional_fasta_dir=user_dir,
            work_dir=work_dir,
            markers=marker_names,
            muscle_command=muscle_command,
            iqtree_command=iqtree_command,
            mrbayes_command=mrbayes_command,
            threads=threads,
            bootstrap=bootstrap,
            ngen=ngen,
            skip_ml=skip_ml,
            skip_bi=skip_bi,
            progress=progress,
        )
        require_files([accessions_dir / "accession_all.csv", accessions_dir / "accession_selected.csv"])

        workflow_summary = [
            "Cladistica full workflow",
            "",
            f"Genus: {genus}",
            f"Outgroups: {'; '.join(outgroups) if outgroups else '(none)'}",
            f"Markers: {', '.join(marker_names)}",
            f"Maximum samples per taxon: {max_samples_per_taxon}",
            f"IQ-TREE bootstrap: {bootstrap}",
            f"MrBayes generations: {ngen}",
            f"Skip ML: {skip_ml}",
            f"Skip BI: {skip_bi}",
        ]
        (work_dir / "workflow_summly.txt").write_text("\n".join(workflow_summary) + "\n", encoding="utf-8")
        if progress:
            progress.start("package", "Collecting final files")
        packaged = package_run_outputs(
            run_dir=work_dir,
            output_dir=run_dir,
            expect_ml=not skip_ml,
            expect_bi=not skip_bi,
        )
        if progress:
            progress.succeed("package", f"{packaged.records_written} files")

    return PipelineResult(
        str(run_dir),
        records_written=accession_result.records_written,
        warnings=packaged.warnings,
    )


def run_from_selected_accessions(
    *,
    project_dir: Path,
    accession_selected: Path,
    accession_all: Path | None,
    additional_fasta_dir: Path | None,
    output_dir: Path | None,
    email: str,
    api_key: str = "",
    markers: list[str] | None = None,
    muscle_command: str = "",
    iqtree_command: str = "",
    mrbayes_command: str = "",
    threads: str = "AUTO",
    bootstrap: int = 1000,
    ngen: int = 1_000_000,
    skip_ml: bool = False,
    skip_bi: bool = False,
    progress: PipelineProgress | None = None,
) -> PipelineResult:
    run_dir = output_dir or default_run_dir(project_dir)
    marker_names = list(marker_map(markers))
    with tempfile.TemporaryDirectory(prefix="cladistica_") as temporary:
        work_dir = Path(temporary)
        accessions_dir = work_dir / "accessions"
        fasta_dir = accessions_dir / "fasta_by_marker"
        accessions_dir.mkdir(parents=True)
        shutil.copy2(accession_selected, accessions_dir / "accession_selected.csv")
        if accession_all:
            shutil.copy2(accession_all, accessions_dir / "accession_all.csv")

        download_result = download_fasta_from_accession_table(
            accession_table=accession_selected,
            output_dir=fasta_dir,
            email=email,
            api_key=api_key,
            markers=marker_names,
            progress=progress,
        )
        shutil.copy2(fasta_dir / "summly.txt", accessions_dir / "summly.txt")
        shutil.copy2(fasta_dir / "run.log", accessions_dir / "run.log")
        user_dir = additional_fasta_dir or (work_dir / "empty_fasta")
        user_dir.mkdir(parents=True, exist_ok=True)
        concatenated = run_marker_fasta_stages(
            primary_fasta_dir=fasta_dir,
            additional_fasta_dir=user_dir,
            work_dir=work_dir,
            markers=marker_names,
            muscle_command=muscle_command,
            iqtree_command=iqtree_command,
            mrbayes_command=mrbayes_command,
            threads=threads,
            bootstrap=bootstrap,
            ngen=ngen,
            skip_ml=skip_ml,
            skip_bi=skip_bi,
            progress=progress,
        )
        workflow_summary = [
            "Cladistica resumed from accession_selected.csv",
            "",
            f"Selected accession table: {accession_selected}",
            f"Accession survey table: {accession_all or '(not supplied)'}",
            f"Additional FASTA: {additional_fasta_dir or '(none)'}",
            f"Markers: {', '.join(marker_names)}",
            f"IQ-TREE bootstrap: {bootstrap}",
            f"MrBayes generations: {ngen}",
        ]
        (work_dir / "workflow_summly.txt").write_text("\n".join(workflow_summary) + "\n", encoding="utf-8")
        if progress:
            progress.start("package", "Collecting final files")
        packaged = package_run_outputs(
            run_dir=work_dir,
            output_dir=run_dir,
            expect_accessions=accession_all is not None,
            expect_ml=not skip_ml,
            expect_bi=not skip_bi,
        )
        if progress:
            progress.succeed("package", f"{packaged.records_written} files")
    return PipelineResult(
        str(run_dir),
        records_written=concatenated.records_written,
        warnings=[*download_result.warnings, *packaged.warnings],
    )


def run_from_fasta(
    *,
    project_dir: Path,
    fasta_dir: Path | None,
    additional_fasta_dir: Path | None,
    concatenated_fasta: Path | None,
    partition_file: Path | None,
    output_dir: Path | None,
    markers: list[str] | None = None,
    muscle_command: str = "",
    iqtree_command: str = "",
    mrbayes_command: str = "",
    threads: str = "AUTO",
    bootstrap: int = 1000,
    ngen: int = 1_000_000,
    skip_ml: bool = False,
    skip_bi: bool = False,
    progress: PipelineProgress | None = None,
) -> PipelineResult:
    run_dir = output_dir or default_run_dir(project_dir)
    marker_names = list(marker_map(markers))
    if bool(fasta_dir) == bool(concatenated_fasta):
        raise ValueError("Provide exactly one of fasta_dir or concatenated_fasta.")
    if partition_file and not concatenated_fasta:
        raise ValueError("partition_file can only be used with concatenated_fasta.")
    if additional_fasta_dir and not fasta_dir:
        raise ValueError("additional_fasta_dir can only be used with marker FASTA input.")

    with tempfile.TemporaryDirectory(prefix="cladistica_") as temporary:
        work_dir = Path(temporary)
        concat_dir = work_dir / "concatenated"
        tree_dir = work_dir / "trees"
        expect_marker_fasta = bool(fasta_dir)

        if fasta_dir:
            user_dir = additional_fasta_dir or (work_dir / "empty_fasta")
            user_dir.mkdir(parents=True, exist_ok=True)
            concatenated = run_marker_fasta_stages(
                primary_fasta_dir=fasta_dir,
                additional_fasta_dir=user_dir,
                work_dir=work_dir,
                markers=marker_names,
                muscle_command=muscle_command,
                iqtree_command=iqtree_command,
                mrbayes_command=mrbayes_command,
                threads=threads,
                bootstrap=bootstrap,
                ngen=ngen,
                skip_ml=skip_ml,
                skip_bi=skip_bi,
                progress=progress,
            )
            sample_count = concatenated.records_written
            source_description = f"Marker FASTA directory: {fasta_dir}"
        else:
            concat_dir.mkdir(parents=True)
            destination = concat_dir / "concatenated_cpDNA.fasta"
            shutil.copy2(concatenated_fasta, destination)
            records = read_fasta(destination)
            if not records:
                raise ValueError(f"No FASTA records found: {concatenated_fasta}")
            lengths = {len(sequence) for sequence in records.values()}
            if len(lengths) != 1:
                raise ValueError("The concatenated FASTA contains inconsistent sequence lengths.")
            if partition_file:
                rows = parse_partition_file(partition_file, next(iter(lengths)))
                write_tsv(
                    concat_dir / "partition_coordinates.tsv",
                    rows,
                    ["marker", "start", "end", "length"],
                )
            run_tree_analyses(
                input_dir=concat_dir,
                output_dir=tree_dir,
                threads=threads,
                bootstrap=bootstrap,
                skip_ml=skip_ml,
                skip_bi=skip_bi,
                iqtree_command=iqtree_command,
                mrbayes_command=mrbayes_command,
                ngen=ngen,
                nruns=2,
                progress=progress,
            )
            require_files(inference_outputs(tree_dir, skip_ml=skip_ml, skip_bi=skip_bi))
            sample_count = len(records)
            source_description = f"Concatenated FASTA: {concatenated_fasta}"

        workflow_summary = [
            "Cladistica resumed from FASTA",
            "",
            source_description,
            f"Additional FASTA: {additional_fasta_dir or '(none)'}",
            f"Partition file: {partition_file or '(automatic single partition or marker-derived partitions)'}",
            f"IQ-TREE bootstrap: {bootstrap}",
            f"MrBayes generations: {ngen}",
        ]
        (work_dir / "workflow_summly.txt").write_text("\n".join(workflow_summary) + "\n", encoding="utf-8")
        if progress:
            progress.start("package", "Collecting final files")
        packaged = package_run_outputs(
            run_dir=work_dir,
            output_dir=run_dir,
            expect_accessions=False,
            expect_marker_fasta=expect_marker_fasta,
            expect_ml=not skip_ml,
            expect_bi=not skip_bi,
        )
        if progress:
            progress.succeed("package", f"{packaged.records_written} files")
    return PipelineResult(str(run_dir), records_written=sample_count, warnings=packaged.warnings)
