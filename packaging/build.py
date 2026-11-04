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
import sysconfig
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from cladistica import __version__  # noqa: E402


def copy_runtime_licenses(bundle: Path) -> None:
    """Preserve interpreter and native-library notices without confusing an EULA with PSF."""
    destination = bundle / "THIRD_PARTY_LICENSES"
    destination.mkdir(exist_ok=True)
    base = Path(sys.base_prefix)
    python_license = None
    for directory in (Path(sysconfig.get_path("stdlib")), base):
        for name in ("LICENSE.txt", "LICENSE"):
            source = directory / name
            if source.is_file() and "Python Software Foundation" in source.read_text(encoding="utf-8", errors="replace"):
                python_license = source
                break
        if python_license is not None:
            break
    if python_license is None:
        raise RuntimeError("Python PSF license was not found; provide it before distributing this build.")
    shutil.copy2(python_license, destination / "Python-LICENSE.txt")
    for name in ("LICENSE.txt", "LICENSE"):
        source = base / name
        if source.is_file() and source.resolve() != python_license.resolve():
            shutil.copy2(source, destination / "Python-distribution-LICENSE.txt")
            break

    metadata_dir = base / "conda-meta"
    if not metadata_dir.is_dir():
        return
    binary_names = {path.name for path in (bundle / "_internal").iterdir()
                    if path.is_file() and (path.suffix.lower() in (".dylib", ".dll") or ".so" in path.name)}
    for metadata in sorted(metadata_dir.glob("*.json")):
        record = json.loads(metadata.read_text(encoding="utf-8"))
        matched = sorted({Path(name).name for name in record.get("files", [])} & binary_names)
        if not matched:
            continue
        cache = Path(record.get("extracted_package_dir") or record.get("link", {}).get("source")
                     or base / "pkgs" / metadata.stem)
        notices = cache / "info" / "licenses"
        target = destination / "conda-runtime" / metadata.stem
        about_path = cache / "info" / "about.json"
        about = json.loads(about_path.read_text(encoding="utf-8")) if about_path.is_file() else {}
        if notices.is_dir():
            shutil.copytree(notices, target, dirs_exist_ok=True)
        elif record.get("name") in ("sqlite", "libsqlite") and record.get("license") == "blessing":
            # Some SQLite packages provide their public-domain declaration and
            # upstream license URL in about.json instead of info/licenses.
            if not about.get("license_url") or "public-domain" not in about.get("description", "").lower():
                raise RuntimeError(f"Missing public-domain declaration for {metadata.stem}")
            target.mkdir(parents=True, exist_ok=True)
        else:
            raise RuntimeError(f"Missing native-library license notices for {metadata.stem}")
        summary = {key: record.get(key) or about.get(key)
                   for key in ("name", "version", "build", "license", "home", "license_url")}
        summary["bundled_library_names"] = matched
        if record.get("name") in ("sqlite", "libsqlite"):
            summary["upstream_description"] = about.get("description", "")
        (target / "PACKAGE_LICENSE_INFO.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    build_env = os.environ.copy()
    build_env["PYINSTALLER_CONFIG_DIR"] = str(ROOT / ".pyinstaller-cache")
    subprocess.run([
        sys.executable, "-m", "PyInstaller", "--noconfirm", "--clean", "--onedir",
        "--console", "--name", "cladistica", "--paths", str(ROOT),
        "--collect-submodules", "Bio.SeqIO", "--collect-submodules", "Bio.AlignIO",
        "--exclude-module", "matplotlib", "--exclude-module", "IPython",
        "--exclude-module", "tkinter", "--exclude-module", "pytest",
        "--exclude-module", "readline",
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
    copy_runtime_licenses(bundle)
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, stderr=subprocess.DEVNULL).strip()
        uncommitted_changes = bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True, stderr=subprocess.DEVNULL).strip())
    except (OSError, subprocess.CalledProcessError):
        commit, uncommitted_changes = "unavailable", None
    (bundle / "BUILD_INFO.json").write_text(json.dumps({
        "cladistica_version": __version__, "python": sys.version,
        "commit": commit, "uncommitted_changes": uncommitted_changes,
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
