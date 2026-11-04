"""Durable MrBayes runs and conservative checkpoint continuation.

Keep original output basenames: MrBayes associates checkpoint, tree, and
parameter samples by that basename and by the ID stored in those files.
"""
from __future__ import annotations

import hashlib
import csv
import json
import os
from pathlib import Path
import re
import shutil
import socket
import subprocess
from datetime import datetime, timezone
from collections.abc import Callable

from . import __version__
from .models import PipelineResult


ANALYSIS = "mrbayes_analysis.nex"
MANIFEST = "resume_manifest.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def save_manifest(directory: Path, state: dict) -> None:
    state["updated_at"] = datetime.now(timezone.utc).isoformat()
    temporary = directory / (MANIFEST + ".tmp")
    temporary.write_text(json.dumps(state, indent=2) + "\n", encoding="utf-8")
    temporary.replace(directory / MANIFEST)


def mrbayes_version(command: str) -> str:
    result = subprocess.run(
        [command], input="version;\nquit;\n", capture_output=True,
        text=True, timeout=30, check=False,
    )
    match = re.search(r"MrBayes\s+(3\.\d+\.\d+[a-z]?)", result.stdout + result.stderr)
    if result.returncode or not match:
        raise RuntimeError("Could not identify MrBayes 3.x; check the executable and its runtime dependencies.")
    return match.group(1)


def process_is_running(pid: int) -> bool:
    if os.name == "nt":
        # os.kill(pid, 0) terminates a process on Windows; query tasklist instead.
        result = subprocess.run(["tasklist", "/FI", f"PID eq {pid}", "/FO", "CSV", "/NH"],
                                capture_output=True, text=True, check=True, timeout=10)
        return any(len(row) > 1 and row[1] == str(pid) for row in csv.reader(result.stdout.splitlines()))
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    return True


def resolve_command(command: str) -> str:
    found = shutil.which(command)
    if not found:
        raise ValueError(f"MrBayes executable not found: {command}")
    return str(Path(found).resolve())


def checkpoint_generation(directory: Path) -> int:
    checkpoint = directory / (ANALYSIS + ".ckp")
    if not checkpoint.is_file() or not checkpoint.stat().st_size:
        raise ValueError("Checkpoint .ckp is missing or empty. A new analysis is required; .p files alone cannot resume BI.")
    with checkpoint.open(encoding="utf-8") as handle:
        header = handle.read(4096)
    match = re.search(r"\[generation:\s*(\d+)\]", header, re.I)
    if not match:
        raise ValueError("Cannot read the checkpoint generation; keep the original .ckp file.")
    return int(match.group(1))


def sample_paths(directory: Path, nruns: int) -> list[Path]:
    suffixes = [f"run{i}.{kind}" for i in range(1, nruns + 1) for kind in ("p", "t")]
    if nruns == 1:
        suffixes = ["p", "t"]
    return [directory / f"{ANALYSIS}.{suffix}" for suffix in suffixes]


def last_parameter_generation(path: Path) -> int:
    last = -1
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            columns = line.split()
            if columns and columns[0].isdigit():
                last = int(columns[0])
    return last


def verify_bundle(directory: Path, state: dict) -> int:
    if state.get("schema") != 1:
        raise ValueError("Unsupported BI resume manifest. Use a bundle created by Cladistica 0.2 or later.")
    if sha256(directory / ANALYSIS) != state["analysis_sha256"]:
        raise ValueError("BI input/model/settings have changed. Refusing to append incompatible samples.")
    # Finished/interrupted runs have a complete inventory; abruptly stopped runs
    # still retain the input fingerprint and are checked by MrBayes itself.
    for name, expected in state.get("files", {}).items():
        if Path(name).name != name:
            raise ValueError("Invalid filename in BI manifest.")
        path = directory / name
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f"BI bundle file is missing or changed: {name}")
    generation = checkpoint_generation(directory)
    checkpoint_text = (directory / (ANALYSIS + ".ckp")).read_text(encoding="utf-8")
    identity = re.search(r"\[ID:\s*([^\]]+)\]", checkpoint_text)
    for path in sample_paths(directory, int(state["nruns"])):
        if not path.is_file() or not path.stat().st_size:
            raise ValueError(f"Required BI samples are missing: {path.name}")
        with path.open(encoding="utf-8") as handle:
            header = handle.read(4096)
        sample_id = re.search(r"\[ID:\s*([^\]]+)\]", header)
        if not identity or not sample_id or identity.group(1).strip() != sample_id.group(1).strip():
            raise ValueError(f"Checkpoint/sample IDs do not match: {path.name}")
        if path.suffix == ".p" and last_parameter_generation(path) < generation:
            raise ValueError(f"Parameter samples end before checkpoint: {path.name}")
    return generation


def execute_bayes(
    directory: Path, command: str, state: dict,
    runner: Callable, on_line: Callable[[str], None] | None = None,
) -> None:
    state.update(status="running", pid=os.getpid(), host=socket.gethostname(), files={})
    save_manifest(directory, state)
    try:
        runner([command, ANALYSIS], directory, directory / "run_mrbayes.log", on_line=on_line)
        for path in sample_paths(directory, int(state["nruns"])):
            if not path.is_file() or not path.stat().st_size:
                raise RuntimeError(f"MrBayes did not produce {path.name}; inspect run_mrbayes.log.")
            if path.suffix == ".p" and last_parameter_generation(path) < state["ngen"]:
                raise RuntimeError("MrBayes did not reach the target generation; samples have been preserved.")
        if not (directory / (ANALYSIS + ".con.tre")).is_file():
            raise RuntimeError("MrBayes did not produce a consensus tree; samples have been preserved.")
        checkpoint_generation(directory)
        state["status"] = "complete"
    except BaseException as error:
        state.update(status="interrupted", error=str(error) or type(error).__name__)
        raise
    finally:
        state["files"] = {
            path.name: sha256(path) for path in sorted(directory.iterdir())
            if path.is_file() and path.name.startswith(ANALYSIS)
        }
        save_manifest(directory, state)


def start_bayes(directory: Path, analysis: Path, command: str, *, ngen: int, nruns: int,
                runner: Callable, on_line: Callable[[str], None] | None = None) -> None:
    if directory.exists() and any(directory.iterdir()):
        raise ValueError("BI output already exists. Use resume-bi or a new output directory to protect previous samples.")
    command = resolve_command(command)
    version = mrbayes_version(command)
    directory.mkdir(parents=True, exist_ok=True)
    shutil.copy2(analysis, directory / ANALYSIS)
    state = dict(schema=1, cladistica_version=__version__, mrbayes_version=version,
                 command=command, ngen=ngen, nruns=nruns,
                 analysis_sha256=sha256(directory / ANALYSIS), history=[])
    execute_bayes(directory, command, state, runner, on_line)


def resume_bayes(*, bundle: Path, output_dir: Path, ngen: int, command: str = "mb") -> PipelineResult:
    """Copy a stopped run, validate it, and append to the new total generation."""
    from .trees import run_checked

    bundle = bundle.resolve()
    if (bundle / "BI_resume").is_dir():
        bundle = bundle / "BI_resume"
    manifest = bundle / MANIFEST
    if not manifest.is_file():
        raise ValueError("No resume_manifest.json found. Select BI_resume or _work/trees/03_mrbayes from Cladistica 0.2+.")
    state = json.loads(manifest.read_text(encoding="utf-8"))
    if state.get("status") == "running":
        if state.get("host") != socket.gethostname():
            raise ValueError("Source run has an unfinished status on another host. Recover it on the original host first.")
        if process_is_running(int(state["pid"])):
            raise ValueError("The source BI run is still active. Stop it before copying/resuming its samples.")
    generation = verify_bundle(bundle, state)
    if ngen <= generation or ngen % 1000:
        raise ValueError(f"--ngen is the total target, must exceed checkpoint generation {generation}, and must be a multiple of 1000.")
    command = resolve_command(command)
    if mrbayes_version(command) != state["mrbayes_version"]:
        raise ValueError("Use the same MrBayes version as the original run for checkpoint continuation.")
    output_dir = output_dir.resolve()
    if bundle == output_dir or bundle in output_dir.parents:
        raise ValueError("Resume output must be outside the source bundle.")
    output_dir.mkdir(parents=True, exist_ok=False)
    destination = output_dir / "BI_resume"
    shutil.copytree(bundle, destination)
    # Preserve earlier command log and analysis script before changing only the
    # target/append controls. All data, models, priors and chain settings persist.
    history = destination / "history"
    history.mkdir(exist_ok=True)
    index = len(state.get("history", [])) + 1
    shutil.copy2(destination / ANALYSIS, history / f"analysis_{index}.nex")
    if (destination / "run_mrbayes.log").exists():
        shutil.copy2(destination / "run_mrbayes.log", history / f"run_{index}.log")
    text = (destination / ANALYSIS).read_text(encoding="utf-8")
    text, count = re.subn(r"\bngen=\d+", f"ngen={ngen}", text)
    if count != 1:
        raise ValueError("Expected one ngen setting in the original analysis.")
    text = text.replace("  mcmc;", "  mcmc append=yes;")
    if "  mcmc append=yes;" not in text:
        raise ValueError("Unsupported MCMC command; no samples were modified.")
    (destination / ANALYSIS).write_text(text, encoding="utf-8")
    state.setdefault("history", []).append(dict(source=str(bundle), checkpoint_generation=generation,
                                                previous_target=state["ngen"], target=ngen))
    state.update(ngen=ngen, analysis_sha256=sha256(destination / ANALYSIS), command=command)
    execute_bayes(destination, command, state, run_checked)
    shutil.copy2(destination / (ANALYSIS + ".con.tre"), output_dir / "BI.tre")
    shutil.copy2(destination / ANALYSIS, output_dir / "BI.nex")
    for path in sample_paths(destination, state["nruns"]):
        if path.suffix == ".p":
            shutil.copy2(path, output_dir / path.name.removeprefix(ANALYSIS + "."))
    return PipelineResult(str(output_dir), records_written=1)
