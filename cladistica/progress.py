from __future__ import annotations

import sys
import threading
from collections.abc import Iterable
from dataclasses import dataclass
from typing import TextIO


FLOWER_FRAMES = [
    ("       .       ", "       |       ", "      / \\      "),
    ("      \\|/      ", "       |       ", "      / \\      "),
    ("    -- * --    ", "      \\|/      ", "      / \\      "),
    ("     \\ * /     ", "    -- | --    ", "      / \\      "),
]


@dataclass
class Stage:
    key: str
    label: str
    state: str = "pending"
    detail: str = ""


class PipelineProgress:
    def __init__(
        self,
        stages: Iterable[tuple[str, str]],
        *,
        enabled: bool = True,
        stream: TextIO | None = None,
    ) -> None:
        self.stream = stream or sys.stderr
        self.enabled = enabled
        self.animated = enabled and self.stream.isatty()
        self.stages = {key: Stage(key, label) for key, label in stages}
        self._frame = 0
        self._rendered = False
        self._stop = threading.Event()
        self._lock = threading.Lock()
        self._thread: threading.Thread | None = None
        self._last_plain_message = ""

    def __enter__(self) -> PipelineProgress:
        if self.animated:
            self.stream.write("\033[?25l")
            self._render()
            self._thread = threading.Thread(target=self._animate, daemon=True)
            self._thread.start()
        return self

    def __exit__(self, exc_type: object, exc: object, traceback: object) -> None:
        if exc is not None:
            active = next((stage for stage in self.stages.values() if stage.state == "running"), None)
            if active:
                self.fail(active.key, str(exc))
        self.close()

    def start(self, key: str, detail: str = "") -> None:
        self._set(key, "running", detail)

    def update(self, key: str, detail: str) -> None:
        self._set(key, "running", detail)

    def succeed(self, key: str, detail: str = "") -> None:
        self._set(key, "success", detail)

    def skip(self, key: str, detail: str = "") -> None:
        self._set(key, "skipped", detail)

    def fail(self, key: str, detail: str = "") -> None:
        self._set(key, "failed", detail)

    def close(self) -> None:
        self._stop.set()
        if self._thread:
            self._thread.join(timeout=1)
        if self.animated:
            with self._lock:
                self._render()
                self.stream.write("\033[?25h")
                self.stream.flush()

    def _set(self, key: str, state: str, detail: str) -> None:
        stage = self.stages.get(key)
        if not stage or not self.enabled:
            return
        with self._lock:
            stage.state = state
            stage.detail = detail
            if self.animated:
                self._render()
            else:
                marker = {
                    "running": ">>",
                    "success": "OK",
                    "skipped": "--",
                    "failed": "!!",
                }.get(state, "  ")
                message = f"[{marker}] {stage.label}"
                if detail:
                    message += f": {detail}"
                if message != self._last_plain_message:
                    self.stream.write(message + "\n")
                    self.stream.flush()
                    self._last_plain_message = message

    def _animate(self) -> None:
        while not self._stop.wait(0.18):
            with self._lock:
                self._frame = (self._frame + 1) % len(FLOWER_FRAMES)
                self._render()

    def _render(self) -> None:
        lines = [
            *FLOWER_FRAMES[self._frame],
            "  Cladistica is growing",
        ]
        markers = {
            "pending": "[  ]",
            "running": "[>>]",
            "success": "[OK]",
            "skipped": "[--]",
            "failed": "[!!]",
        }
        for stage in self.stages.values():
            line = f"{markers[stage.state]} {stage.label}"
            if stage.detail:
                line += f" - {stage.detail}"
            lines.append(line)
        if self._rendered:
            self.stream.write(f"\033[{len(lines)}A")
        for line in lines:
            self.stream.write(f"\r\033[2K{line}\n")
        self.stream.flush()
        self._rendered = True
