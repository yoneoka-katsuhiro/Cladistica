"""Build and smoke-test a native CLI archive on the current operating system."""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from cladistica import __version__  # noqa: E402


def main() -> None:
    build_env = os.environ.copy()
    build_env["PYINSTALLER_CONFIG_DIR"] = str(ROOT / ".pyinstaller-cache")
    subprocess.run([
        sys.executable, "-m", "PyInstaller", "--noconfirm", "--clean", "--onedir",
        "--console", "--name", "cladistica", "--paths", str(ROOT),
        "--collect-submodules", "Bio.SeqIO", "--collect-submodules", "Bio.AlignIO",
        "--exclude-module", "matplotlib", "--exclude-module", "IPython",
        "--exclude-module", "tkinter", "--exclude-module", "pytest",
        str(ROOT / "packaging/entrypoint.py"),
    ], cwd=ROOT, check=True, env=build_env)
    bundle = ROOT / "dist/cladistica"
    executable = bundle / ("cladistica.exe" if sys.platform == "win32" else "cladistica")
    subprocess.run([str(executable), "--version"], check=True)
    subprocess.run([str(executable), "self-check"], check=True)
    with tempfile.TemporaryDirectory(prefix="cladistica_frozen_") as temporary:
        workspace = Path(temporary)
        inputs = workspace / "input"
        inputs.mkdir()
        shutil.copy2(ROOT / "examples/tiny.fasta", inputs / "tiny.fasta")
        subprocess.run([str(executable), "tree", "--input", str(inputs), "--output",
                        str(workspace / "output"), "--dry-run", "--no-progress"], cwd=workspace, check=True)
        assert (workspace / "output/00_inputs/mrbayes_analysis.nex").is_file()
    for name in ("README.md", "README-ja.md", "LICENSE", "CITATIONS.md", "CITATION.cff", "CHANGELOG.md"):
        shutil.copy2(ROOT / name, bundle / name)
    for name in ("docs", "examples"):
        shutil.copytree(ROOT / name, bundle / name, dirs_exist_ok=True)
    licenses = bundle / "THIRD_PARTY_LICENSES"
    licenses.mkdir(exist_ok=True)
    dependencies = []
    for distribution in sorted(importlib.metadata.distributions(), key=lambda dist: dist.metadata["Name"].lower()):
        name = distribution.metadata["Name"]
        dependencies.append({"name": name, "version": distribution.version})
        for item in distribution.files or []:
            if any(part.lower().startswith(("license", "copying", "notice")) for part in item.parts):
                source = Path(distribution.locate_file(item))
                if source.is_file():
                    target = licenses / name / str(item).replace("../", "parent/")
                    target.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(source, target)
    for path in (Path(sys.base_prefix) / "LICENSE.txt", Path(sys.base_prefix) / "LICENSE"):
        if path.is_file():
            shutil.copy2(path, licenses / "Python-LICENSE.txt")
            break
    (bundle / "BUILD_INFO.json").write_text(json.dumps({
        "cladistica_version": __version__, "python": sys.version,
        "platform": platform.platform(), "architecture": platform.machine(),
        "build_environment_packages": dependencies,
        "external_tools_included": False,
    }, indent=2) + "\n", encoding="utf-8")
    architecture = {"aarch64": "arm64", "AMD64": "x64", "x86_64": "x64"}.get(platform.machine(), platform.machine())
    system = "macOS" if sys.platform == "darwin" else "Windows" if sys.platform == "win32" else "Linux"
    basename = f"Cladistica-{__version__}-{system}-{architecture}"
    archive = Path(shutil.make_archive(str(ROOT / "dist" / basename), "zip", root_dir=ROOT / "dist", base_dir="cladistica"))
    checksum = hashlib.sha256(archive.read_bytes()).hexdigest()
    archive.with_suffix(".zip.sha256").write_text(f"{checksum}  {archive.name}\n", encoding="utf-8")
    print(archive)


if __name__ == "__main__":
    main()
