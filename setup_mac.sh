#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV_DIR="${PROJECT_DIR}/.venv"
PYTHON_BIN="${PYTHON_BIN:-}"

if [[ -z "${PYTHON_BIN}" ]]; then
  for candidate in python3.13 python3.12 python3.11 python3.10 python3 python; do
    if command -v "${candidate}" >/dev/null 2>&1; then
      if "${candidate}" -c 'import sys; raise SystemExit(0 if sys.version_info >= (3, 10) else 1)'; then
        PYTHON_BIN="${candidate}"
        break
      fi
    fi
  done
fi

if [[ -z "${PYTHON_BIN}" ]] || ! "${PYTHON_BIN}" -c 'import sys; raise SystemExit(0 if sys.version_info >= (3, 10) else 1)'; then
  echo "ERROR: Cladistica requires Python 3.10 or newer." >&2
  echo "Install a newer Python or run with PYTHON_BIN=/path/to/python bash setup_mac.sh" >&2
  exit 1
fi

"${PYTHON_BIN}" -m venv "${VENV_DIR}"
"${VENV_DIR}/bin/python" -m pip install --upgrade pip
"${VENV_DIR}/bin/python" -m pip install -e "${PROJECT_DIR}"

echo "Cladistica environment ready:"
echo "  bash run_cladistica.sh examples"
