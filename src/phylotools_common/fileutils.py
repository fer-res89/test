"""Generic, tree-library-agnostic file/path utilities shared across tools."""
from __future__ import annotations

import re
from pathlib import Path
from typing import List, Tuple


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def newly_created_dir_was(path: Path) -> bool:
    return not path.exists()


def describe_and_print_new_dirs(new_dirs: List[Tuple[Path, str]]) -> None:
    if not new_dirs:
        return
    print("\nNew directories created:")
    for p, desc in new_dirs:
        print(f" - {p.as_posix()}: {desc}")


def resolve_path(p: str | None) -> Path | None:
    if not p:
        return None
    return Path(p).expanduser().resolve()


def sanitize_id(s: str) -> str:
    """
    Make an ID safe for filenames while preserving uniqueness reasonably well.
    """
    s = s.strip()
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"[\\/]", "_", s)
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", s)
    return s.strip("_") or "tree"
