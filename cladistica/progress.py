from __future__ import annotations

import sys
import threading
import time
from collections.abc import Iterable
from dataclasses import dataclass
from typing import TextIO


STREAM_WIDTH = 40
STREAM_HEIGHT = 5
ANIMATION_INTERVAL_SECONDS = 0.10
MOVE_EVERY_FRAMES = 2
FALL_EVERY_FRAMES = 3
MINIMUM_FALL_SECONDS = 2.4
PROGRESS_BAR_WIDTH = 32

DEFAULT_STAGE_ITEMS = {
    "survey": (
        "NCBI GenBank",
        "DDBJ",
        "ENA",
        "INSDC",
        "Entrez",
        "accession",
        "taxonomy",
    ),
    "accessions": (
        "NCBI GenBank",
        "DDBJ",
        "ENA",
        "INSDC",
        "Entrez",
        "accession",
        "1 taxon = 1 sample",
        "marker coverage",
    ),
    "download": (
        "GenBank FASTA",
        "accession",
        "CDS",
        "noncoding",
        "feature extraction",
        "sequence QC",
    ),
    "align": (
        "MUSCLE",
        "multiple sequence alignment",
        "terminal gap",
        "ambiguous base",
        "rbcL",
        "matK",
    ),
    "concat": (
        "concatenation",
        "partition coordinates",
        "missing data = ?",
        "charset",
        "aligned matrix",
    ),
    "model": (
        "ModelFinder",
        "IQ-TREE",
        "JC",
        "F81",
        "K2P",
        "K80",
        "HKY",
        "HKY85",
        "TN",
        "TN93",
        "K3P",
        "K81",
        "K81u",
        "TPM2",
        "TPM2u",
        "TPM3",
        "TPM3u",
        "TIM",
        "TIMe",
        "TIM2",
        "TIM2e",
        "TIM3",
        "TIM3e",
        "TVM",
        "TVMe",
        "SYM",
        "GTR",
        "+F",
        "+I",
        "+G4",
        "+R",
        "BIC",
        "AIC",
        "AICc",
    ),
    "bootstrap": (
        "bootstrap replicate",
        "resampled alignment",
        "bootstrap consensus",
        "branch support",
        "IQ-TREE",
    ),
    "ml": (
        "maximum likelihood",
        "tree search",
        "log-likelihood",
        "NNI",
        "candidate tree",
        "best tree",
        "IQ-TREE",
        "bootstrap",
    ),
    "bi": (
        "MrBayes",
        "MCMC",
        "run1",
        "run2",
        "chain 1",
        "chain 2",
        "chain 3",
        "chain 4",
        "generation",
        "split frequency",
        "ESS",
        "PSRF",
    ),
    "bi_summary": (
        "sump",
        "sumt",
        "parameter summary",
        "tree consensus",
        "split frequency",
        "ESS",
        "PSRF",
        "run1.p",
        "run2.p",
        "BI.tre",
    ),
    "package": (
        "accession_all.csv",
        "accession_selected.csv",
        "concatenated.fasta",
        "partitions.txt",
        "ML.tre",
        "BI.tre",
        "run1.p",
        "run2.p",
        "summly.txt",
        "run.log",
    ),
}


@dataclass
class Stage:
    key: str
    label: str
    state: str = "pending"
    detail: str = ""
    percent: float | None = None
    approximate: bool = False


@dataclass
class StreamLane:
    text: str
    x: int


def normalize_stream_item(value: object) -> str:
    text = " ".join(str(value).replace("\x1b", "").split())
    return text[:80]


def sequence_windows(sequence: str) -> list[str]:
    bases = "".join(
        character
        for character in sequence.upper()
        if character in "ACGTRYSWKMBDHVN?-"
    )
    if not bases:
        return []
    if len(bases) <= STREAM_WIDTH:
        return [bases]
    starts = (0, (len(bases) - STREAM_WIDTH) // 2, len(bases) - STREAM_WIDTH)
    return list(dict.fromkeys(bases[start : start + STREAM_WIDTH] for start in starts))


class PipelineProgress:
    def __init__(
        self,
        stages: Iterable[tuple[str, str]],
        *,
        enabled: bool = True,
        stream: TextIO | None = None,
    ) -> None:
        stage_list = list(stages)
        self.stream = stream or sys.stderr
        self.enabled = enabled
        self.animated = enabled and self.stream.isatty()
        self.stages = {key: Stage(key, label) for key, label in stage_list}
        self._global_items: list[str] = []
        self._global_known: set[str] = set()
        self._stage_items: dict[str, list[str]] = {
            key: list(DEFAULT_STAGE_ITEMS.get(key, (label,)))
            for key, label in stage_list
        }
        self._stage_known: dict[str, set[str]] = {
            key: set(items) for key, items in self._stage_items.items()
        }
        self._active_key: str | None = None
        self._item_cursor = 0
        self._lanes: list[StreamLane] = []
        self._particles: list[list[str]] | None = None
        self._falling = False
        self._fall_started = 0.0
        self._frame = 0
        self._rendered = False
        self._stop = threading.Event()
        self._lock = threading.Lock()
        self._thread: threading.Thread | None = None
        self._last_plain_message = ""

    def __enter__(self) -> PipelineProgress:
        if self.animated:
            with self._lock:
                self._initialize_lanes()
                self.stream.write("\033[?25l")
                self._render()
            self._thread = threading.Thread(target=self._animate, daemon=True)
            self._thread.start()
        return self

    def __exit__(self, exc_type: object, exc: object, traceback: object) -> None:
        if exc is not None:
            active = next(
                (stage for stage in self.stages.values() if stage.state == "running"),
                None,
            )
            if active:
                self.fail(active.key, str(exc))
            elif self.animated:
                with self._lock:
                    self._start_fall()
        self.close()

    def feed(self, *items: object) -> None:
        self.feed_for(self._active_key, *items)

    def feed_for(self, stage_key: str | None, *items: object) -> None:
        normalized = [normalize_stream_item(item) for item in items]
        with self._lock:
            if stage_key and stage_key in self.stages:
                target = self._stage_items.setdefault(stage_key, [])
                known = self._stage_known.setdefault(stage_key, set())
            else:
                target = self._global_items
                known = self._global_known
            for item in normalized:
                if item and item not in known:
                    known.add(item)
                    target.append(item)

    def feed_sequence(self, label: object, sequence: str) -> None:
        self.feed(label, *sequence_windows(sequence))

    def feed_sequence_for(self, stage_key: str, label: object, sequence: str) -> None:
        self.feed_for(stage_key, label, *sequence_windows(sequence))

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
        with self._lock:
            stage.state = "running"
            stage.percent = min(100.0, max(0.0, float(percent)))
            stage.approximate = approximate
            if detail:
                stage.detail = detail
            if self.animated:
                if key != self._active_key:
                    self._active_key = key
                    self._reset_lanes()
                self._render()

    def close(self) -> None:
        if self._falling and self._thread:
            remaining = MINIMUM_FALL_SECONDS - (time.monotonic() - self._fall_started)
            if remaining > 0:
                time.sleep(remaining)
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
            if state == "running" and stage.percent is None:
                stage.percent = 0.0
            elif state == "success":
                stage.percent = 100.0
                stage.approximate = False
            if self.animated:
                if state == "running" and key != self._active_key:
                    self._active_key = key
                    self._reset_lanes()
                if state == "failed":
                    self._start_fall()
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

    def _initialize_lanes(self) -> None:
        if self._lanes:
            return
        if not self._global_items:
            self._global_items.extend(stage.label for stage in self.stages.values())
            self._global_known.update(self._global_items)
        self._reset_lanes()

    def _reset_lanes(self) -> None:
        self._item_cursor = 0
        self._lanes = [
            StreamLane(self._next_item(), -(row * 9))
            for row in range(STREAM_HEIGHT)
        ]

    def _next_item(self) -> str:
        choices = [
            *self._stage_items.get(self._active_key or "", []),
            *self._global_items,
        ] or ["ACGTN"]
        item = choices[self._item_cursor % len(choices)]
        self._item_cursor += 1
        return item

    def _animate(self) -> None:
        while not self._stop.wait(ANIMATION_INTERVAL_SECONDS):
            with self._lock:
                self._frame += 1
                if self._falling:
                    if self._frame % FALL_EVERY_FRAMES == 0:
                        self._advance_fall()
                elif self._frame % MOVE_EVERY_FRAMES == 0:
                    self._advance_stream()
                self._render()

    def _advance_stream(self) -> None:
        for lane in self._lanes:
            lane.x += 1
            if lane.x >= STREAM_WIDTH:
                lane.text = self._next_item()
                lane.x = -len(lane.text) - 6

    def _stream_canvas(self) -> list[list[str]]:
        canvas = [[" "] * STREAM_WIDTH for _ in range(STREAM_HEIGHT)]
        for row, lane in enumerate(self._lanes):
            for offset, character in enumerate(lane.text):
                x = lane.x + offset
                if 0 <= x < STREAM_WIDTH:
                    canvas[row][x] = character
        return canvas

    def _start_fall(self) -> None:
        if self._falling:
            return
        self._falling = True
        self._fall_started = time.monotonic()
        self._particles = self._stream_canvas()

    def _advance_fall(self) -> None:
        if self._particles is None:
            return
        grid = self._particles
        for y in range(STREAM_HEIGHT - 2, -1, -1):
            columns = range(STREAM_WIDTH)
            if (self._frame // FALL_EVERY_FRAMES + y) % 2:
                columns = range(STREAM_WIDTH - 1, -1, -1)
            for x in columns:
                character = grid[y][x]
                if character == " ":
                    continue
                destinations = [(x, y + 1)]
                direction = -1 if (x + y + self._frame) % 2 else 1
                destinations.extend(
                    [(x + direction, y + 1), (x - direction, y + 1)]
                )
                for next_x, next_y in destinations:
                    if 0 <= next_x < STREAM_WIDTH and grid[next_y][next_x] == " ":
                        grid[next_y][next_x] = character
                        grid[y][x] = " "
                        break

    def _render(self) -> None:
        canvas = self._particles if self._falling and self._particles else self._stream_canvas()
        mode = "FALLING / ANALYSIS STOPPED" if self._falling else "RUNNING"
        lines = [
            f"+{'-' * STREAM_WIDTH}+",
            *(f"|{''.join(row)}|" for row in canvas),
            f"+{'-' * STREAM_WIDTH}+",
            f"Sequence stream [{mode}]",
        ]
        for stage in self.stages.values():
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
            lines.append(line)
        if self._rendered:
            self.stream.write(f"\033[{len(lines)}A")
        for line in lines:
            self.stream.write(f"\r\033[2K{line}\n")
        self.stream.flush()
        self._rendered = True
