# Build and distribution / ビルドと配布

Native packages are command-line distributions. They include the Python runtime and Python dependencies, examples, both manuals, `BUILD_INFO.json`, and third-party license texts. MUSCLE, IQ-TREE and MrBayes remain separate installations.

Build on the target OS; PyInstaller does not cross-compile Windows packages on macOS. The GitHub Actions `Native packages` workflow uses `macos-15` (Apple Silicon) and `windows-2022` (x64), runs the offline suite, builds the CLI, checks `--version` and GenBank/Excel round-trips through `self-check`, and executes a frozen `tree --dry-run` on the example matrix. It uploads ZIP/checksum artifacts; it does not publish a release automatically.

```bash
python -m venv .venv-build
# Activate .venv-build for your operating system.
python -m pip install --upgrade pip
python -m pip install -c packaging/constraints-release.txt -e . "pyinstaller>=6,<7"
python -m unittest discover -s tests -v
python packaging/build.py
```

Release CI uses Python 3.12 and `packaging/constraints-release.txt`; `BUILD_INFO.json` records the actual resolved environment. Local builds with another Python version are development builds and record that version explicitly.

The ZIP and adjacent `.sha256` file appear in `dist/`. `BUILD_INFO.json` records Python, architecture, platform and the installed build-environment package names/versions. It does not expose editable-install filesystem paths. The inventory may contain build-only packages, which are not necessarily part of the executable. License copies are collected from installed distributions; review them before redistributing a customized bundle. The project's own license is MIT.

Before publishing, inspect the platform CI result, download the artifact, check the checksum, and verify its package contents. A local Apple Silicon smoke test and real MrBayes continuation tests do not establish Windows external-tool behavior. Report that distinction in the release notes. The CLI is unsigned/not notarized; it is not an Apple Developer ID or Microsoft Authenticode release. Do not describe it as a GUI installer, Intel/universal Mac build, or self-contained phylogenetics suite.

## 日本語

配布物はコマンドライン版です。Python環境、依存、例題、英日マニュアル、バージョン一覧、第三者ライセンスを同梱します。MUSCLE・IQ-TREE・MrBayesは別導入です。

各OS上でビルドしてください。GitHub ActionsはmacOS 15のarm64とWindows Server 2022のx64でテスト・ビルド・CLI起動・dry-runを実行し、ZIPとSHA-256をartifactへ保存します。リリースは自動公開しません。Macで行った実MrBayes再開試験と、Windowsでの自動ビルド確認は区別して報告してください。物理Windows PCで外部解析ソフトを含む全工程を検証したという意味ではありません。

macOS版は未公証、Windows版は未署名です。全フォルダを展開し、実行ファイルだけを移動しないでください。

## References

- [PyInstaller documentation](https://pyinstaller.org/en/stable/) — native packaging.
- [GitHub runner images](https://github.com/actions/runner-images) — runner labels and architectures. macOS 14 GitHub runners retired in November 2026; the workflow uses macOS 15.
