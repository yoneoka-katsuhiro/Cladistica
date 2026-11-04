"""PyInstaller entry point; use the installed package's public CLI."""
from cladistica.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
