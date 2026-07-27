from __future__ import annotations

import re
import shutil
import signal
import sys
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from types import FrameType
from typing import TextIO


PROGRESS_BAR_WIDTH = 32
COMPACT_PROGRESS_BAR_WIDTH = 20
MINIMUM_PROGRESS_BAR_WIDTH = 6
ERASE_CURRENT_LINE = "\r\x1b[2K"

SHORT_STAGE_LABELS = {
    "accessions": "NCBI survey",
    "download": "Download",
    "align": "Alignment",
    "concat": "Concatenate",
    "model": "ModelFinder",
    "bootstrap": "Bootstrap",
    "ml": "ML tree",
    "bi": "MrBayes",
    "bi_summary": "BI summary",
    "package": "Final output",
}


@dataclass
class Stage:
    key: str
    label: str
    state: str = "pending"
    detail: str = ""
    percent: float | None = None
    approximate: bool = False


def default_terminal_width() -> int:
    return shutil.get_terminal_size(fallback=(80, 24)).columns


def compact_detail(detail: str) -> str:
    """Keep the useful count while removing a repeated unit at narrow widths."""
    return re.sub(
        r"\s+(?:bootstrap\s+)?(?:replicates?|generations?|sequences?|markers?)$",
        "",
        detail.strip(),
        flags=re.IGNORECASE,
    )


def ellipsize(text: str, width: int) -> str:
    if width <= 0:
        return ""
    if len(text) <= width:
        return text
    if width <= 3:
        return "." * width
    return text[: width - 3].rstrip() + "..."


class PipelineProgress:
    """Display one resize-safe progress line for the currently active stage."""

    def __init__(
        self,
        stages: Iterable[tuple[str, str]],
        *,
        enabled: bool = True,
        stream: TextIO | None = None,
        terminal_width_provider: Callable[[], int] | None = None,
    ) -> None:
        self.stream = stream or sys.stderr
        self.enabled = enabled
        self.interactive = enabled and self.stream.isatty()
        self.stages = {key: Stage(key, label) for key, label in stages}
        self._terminal_width_provider = (
            terminal_width_provider or default_terminal_width
        )
        self._active_key: str | None = None
        self._line_open = False
        self._last_terminal_columns: int | None = None
        self._last_progress_signature: tuple[str, int, bool] | None = None
        self._last_plain_message = ""
        self._previous_sigwinch_handler: object | None = None
        self._sigwinch_installed = False
        self._in_render = False
        self._resize_pending = False

    def __enter__(self) -> PipelineProgress:
        self._install_sigwinch_handler()
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
        if not self.interactive:
            return

        columns = self._terminal_columns()
        signature = (key, round(stage.percent), approximate)
        if (
            signature == self._last_progress_signature
            and columns == self._last_terminal_columns
        ):
            return
        self._last_progress_signature = signature
        self._write_interactive(stage, final=False, columns=columns)

    def close(self) -> None:
        if self.enabled and self.interactive and self._line_open:
            self.stream.write(ERASE_CURRENT_LINE + "\n")
            self.stream.flush()
            self._line_open = False
        self._restore_sigwinch_handler()

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

    def _write_interactive(
        self,
        stage: Stage,
        *,
        final: bool,
        columns: int | None = None,
    ) -> None:
        if self._in_render:
            self._resize_pending = True
            return

        self._in_render = True
        self._resize_pending = False
        try:
            current_columns = columns or self._terminal_columns()
            if (
                self._line_open
                and self._last_terminal_columns is not None
                and current_columns != self._last_terminal_columns
            ):
                # End the old terminal-width logical line before drawing at the
                # new width. This prevents terminal reflow from joining both
                # renderings.
                self.stream.write(ERASE_CURRENT_LINE + "\n")
                self._line_open = False

            max_width = max(1, current_columns - 1)
            rendered_line = self._format_stage(stage, max_width)
            self.stream.write(ERASE_CURRENT_LINE + rendered_line)
            if final:
                self.stream.write("\n")
                self._line_open = False
                self._last_terminal_columns = None
            else:
                self._line_open = True
                self._last_terminal_columns = current_columns
            self.stream.flush()
        finally:
            self._in_render = False

        if self._resize_pending and self._line_open and self._active_key:
            self._resize_pending = False
            active = self.stages.get(self._active_key)
            if active:
                self._write_interactive(active, final=False)

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

    def _terminal_columns(self) -> int:
        try:
            return max(2, int(self._terminal_width_provider()))
        except (OSError, TypeError, ValueError):
            return 80

    def _install_sigwinch_handler(self) -> None:
        if not self.interactive or not hasattr(signal, "SIGWINCH"):
            return
        try:
            self._previous_sigwinch_handler = signal.getsignal(signal.SIGWINCH)
            signal.signal(signal.SIGWINCH, self._handle_sigwinch)
            self._sigwinch_installed = True
        except (OSError, RuntimeError, ValueError):
            self._previous_sigwinch_handler = None
            self._sigwinch_installed = False

    def _restore_sigwinch_handler(self) -> None:
        if not self._sigwinch_installed or not hasattr(signal, "SIGWINCH"):
            return
        try:
            signal.signal(signal.SIGWINCH, self._previous_sigwinch_handler)
        except (OSError, RuntimeError, TypeError, ValueError):
            pass
        self._sigwinch_installed = False

    def _handle_sigwinch(
        self,
        signum: int,
        frame: FrameType | None,
    ) -> None:
        if self._in_render:
            self._resize_pending = True
            return
        if not self._line_open or not self._active_key:
            return
        stage = self.stages.get(self._active_key)
        if stage:
            try:
                self._write_interactive(stage, final=False)
            except (OSError, RuntimeError, ValueError):
                # A signal can interrupt an in-progress stream write. The next
                # progress update still detects the changed terminal width.
                self._resize_pending = True

    @staticmethod
    def _format_stage(stage: Stage, max_width: int) -> str:
        percent = stage.percent if stage.percent is not None else 0.0
        if stage.state == "skipped":
            percent_text = "SKIP"
        elif stage.state == "failed":
            percent_text = f"!{percent:.0f}%"
        elif stage.approximate and stage.state == "running":
            percent_text = f"~{percent:.0f}%"
        else:
            percent_text = f"{percent:.0f}%"

        short_label = SHORT_STAGE_LABELS.get(stage.key, stage.label)
        if max_width >= 90:
            label = stage.label
            detail = stage.detail.strip()
            preferred_bar_width = PROGRESS_BAR_WIDTH
        elif max_width >= 55:
            label = short_label
            detail = compact_detail(stage.detail)
            preferred_bar_width = COMPACT_PROGRESS_BAR_WIDTH
        elif max_width >= 25:
            label = short_label
            detail = compact_detail(stage.detail)
            preferred_bar_width = 0
        else:
            label = short_label
            detail = ""
            preferred_bar_width = 0

        text = PipelineProgress._fit_text(
            percent_text,
            label,
            detail,
            max_width,
        )
        if preferred_bar_width:
            available = max_width - len(text) - 3
            if available >= MINIMUM_PROGRESS_BAR_WIDTH:
                bar_width = min(preferred_bar_width, available)
                filled = round(bar_width * percent / 100)
                bar = "=" * filled + "-" * (bar_width - filled)
                text = f"[{bar}] {text}"
        return text[:max_width]

    @staticmethod
    def _fit_text(
        percent_text: str,
        label: str,
        detail: str,
        max_width: int,
    ) -> str:
        base = f"{percent_text} {label}"
        if len(base) > max_width:
            return percent_text[:max_width]

        if not detail:
            return base
        detail_prefix = base + " - "
        available_detail = max_width - len(detail_prefix)
        if available_detail >= 4:
            return detail_prefix + ellipsize(detail, available_detail)
        return base
