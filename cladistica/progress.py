from __future__ import annotations

import sys
from collections.abc import Iterable
from dataclasses import dataclass
from typing import TextIO


PROGRESS_BAR_WIDTH = 32


@dataclass
class Stage:
    key: str
    label: str
    state: str = "pending"
    detail: str = ""
    percent: float | None = None
    approximate: bool = False


class PipelineProgress:
    """Display one quiet progress-bar line for the currently active stage."""

    def __init__(
        self,
        stages: Iterable[tuple[str, str]],
        *,
        enabled: bool = True,
        stream: TextIO | None = None,
    ) -> None:
        self.stream = stream or sys.stderr
        self.enabled = enabled
        self.interactive = enabled and self.stream.isatty()
        self.stages = {key: Stage(key, label) for key, label in stages}
        self._active_key: str | None = None
        self._line_open = False
        self._last_line_length = 0
        self._last_plain_message = ""
        self._last_progress_signature: tuple[str, int, bool] | None = None

    def __enter__(self) -> PipelineProgress:
        return self

    def __exit__(self, exc_type: object, exc: object, traceback: object) -> None:
        if exc is not None:
            active = next(
                (stage for stage in self.stages.values() if stage.state == "running"),
                None,
            )
            if active:
                detail = str(exc)
                if not detail:
                    detail = getattr(exc_type, "__name__", "interrupted")
                self.fail(active.key, detail)
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

    def set_progress(
        self,
        key: str,
        percent: float,
        detail: str = "",
        *,
        approximate: bool = False,
    ) -> None:
        stage = self.stages.get(key)
        if not stage or not self.enabled:
            return
        stage.state = "running"
        stage.percent = min(100.0, max(0.0, float(percent)))
        stage.approximate = approximate
        if detail:
            stage.detail = detail
        self._active_key = key
        if self.interactive:
            signature = (key, round(stage.percent), approximate)
            if signature == self._last_progress_signature:
                return
            self._last_progress_signature = signature
            self._write_interactive(stage, final=False)

    def close(self) -> None:
        if self.enabled and self.interactive and self._line_open:
            self.stream.write("\n")
            self.stream.flush()
            self._line_open = False

    def _set(self, key: str, state: str, detail: str) -> None:
        stage = self.stages.get(key)
        if not stage or not self.enabled:
            return
        stage.state = state
        stage.detail = detail
        if state == "running" and stage.percent is None:
            stage.percent = 0.0
        elif state == "success":
            stage.percent = 100.0
            stage.approximate = False
        self._active_key = key if state == "running" else None

        if self.interactive:
            if state == "running":
                self._last_progress_signature = (
                    key,
                    round(stage.percent or 0.0),
                    stage.approximate,
                )
            else:
                self._last_progress_signature = None
            self._write_interactive(stage, final=state != "running")
        else:
            self._write_plain(stage)

    def _write_interactive(self, stage: Stage, *, final: bool) -> None:
        line = self._format_stage(stage)
        width = max(self._last_line_length, len(line))
        self.stream.write("\r" + line.ljust(width))
        if final:
            self.stream.write("\n")
            self._line_open = False
            self._last_line_length = 0
        else:
            self._line_open = True
            self._last_line_length = width
        self.stream.flush()

    def _write_plain(self, stage: Stage) -> None:
        marker = {
            "running": ">>",
            "success": "OK",
            "skipped": "--",
            "failed": "!!",
        }.get(stage.state, "  ")
        message = f"[{marker}] {stage.label}"
        if stage.detail:
            message += f": {stage.detail}"
        if message != self._last_plain_message:
            self.stream.write(message + "\n")
            self.stream.flush()
            self._last_plain_message = message

    @staticmethod
    def _format_stage(stage: Stage) -> str:
        if stage.state == "skipped":
            bar = "-" * PROGRESS_BAR_WIDTH
            line = f"[{bar}] SKIP {stage.label}"
        else:
            percent = stage.percent if stage.percent is not None else 0.0
            filled = round(PROGRESS_BAR_WIDTH * percent / 100)
            bar = "=" * filled + "-" * (PROGRESS_BAR_WIDTH - filled)
            percent_text = f"{percent:3.0f}%"
            if stage.approximate and stage.state == "running":
                percent_text = f"~{percent:.0f}%"
            if stage.state == "failed":
                percent_text = f"!{percent:.0f}%"
            line = f"[{bar}] {percent_text} {stage.label}"
        if stage.detail:
            line += f" - {stage.detail}"
        return line
