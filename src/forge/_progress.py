"""Shared console progress/spinner utility for forge subcommands.

Provides a Rich-based status line that redraws at a constant wall-clock
rate, independent of whatever work is actually happening.

This matters whenever a step is made of few, slow units of work -- e.g. one
multi-minute BLAST search in `blast-species`, or one slow tree parse in
`filter`. If the display only redraws when a unit of work finishes, it can
appear frozen for minutes at a time even though the script is running
correctly. A background "ticker" thread separates render cadence from work
cadence, so the spinner always looks alive as long as the process is.

Usage:

    console = Console(stderr=True)
    done = 0
    total = len(items)

    def render(spin: str) -> Text:
        return Text(f"{spin} {done}/{total} done")

    with TickingStatus(console, render, interval=0.15) as status:
        for item in items:
            do_slow_thing(item)   # may take a long time; ticker keeps spinning
            done += 1
            status.refresh()      # optional: reflect the change immediately,
                                  # instead of waiting for the next tick

`render_fn` may be called from the ticker thread OR from a manual
`.refresh()` call, so it should only *read* shared state, never mutate it.
"""
from __future__ import annotations

import itertools
import threading
from typing import Callable, Optional

from rich.console import Console
from rich.live import Live

SPINNER_FRAMES = "|/-\\"


class TickingStatus:
    """A Live status line that redraws on a fixed wall-clock interval,
    regardless of how often (or rarely) the underlying work reports in.
    """

    def __init__(
        self,
        console: Console,
        render_fn: Callable[[str], object],
        interval: float = 0.15,
        refresh_per_second: int = 10,
        transient: bool = True,
    ) -> None:
        self._console = console
        self._render_fn = render_fn
        self._interval = interval
        self._spinner_cycle = itertools.cycle(SPINNER_FRAMES)
        self._spinner_lock = threading.Lock()
        self._render_lock = threading.Lock()
        self._stop_event = threading.Event()
        self._thread: Optional[threading.Thread] = None
        self._live = Live(
            self._render(),
            console=console,
            refresh_per_second=refresh_per_second,
            transient=transient,
            auto_refresh=False,
        )

    def _next_spin(self) -> str:
        with self._spinner_lock:
            return next(self._spinner_cycle)

    def _render(self):
        return self._render_fn(self._next_spin())

    def refresh(self) -> None:
        """Force an immediate redraw. Safe to call from any thread."""
        renderable = self._render()
        with self._render_lock:
            self._live.update(renderable, refresh=True)

    def _tick_loop(self) -> None:
        while not self._stop_event.is_set():
            self.refresh()
            self._stop_event.wait(self._interval)

    def __enter__(self) -> "TickingStatus":
        self._live.__enter__()
        self._thread = threading.Thread(
            target=self._tick_loop, name="forge-ticking-status", daemon=True
        )
        self._thread.start()
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self._stop_event.set()
        if self._thread is not None:
            self._thread.join(timeout=self._interval * 4)
        self._live.__exit__(exc_type, exc, tb)
