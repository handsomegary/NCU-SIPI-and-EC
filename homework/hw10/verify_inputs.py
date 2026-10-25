"""Verify the exact supplied HW10 inputs without importing legacy templates."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


def verify(directory: Path | None = None) -> dict:
    directory = directory or Path(__file__).resolve().parent / "instructions"
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
    checked = []
    for item in manifest["files"]:
        path = directory / item["filename"]
        data = path.read_bytes()
        if len(data) != item["bytes"]:
            raise ValueError(f"Size mismatch: {path.name}")
        if hashlib.sha256(data).hexdigest() != item["sha256"]:
            raise ValueError(f"SHA-256 mismatch: {path.name}")
        checked.append(path.name)
    return {"status": "passed", "checked_files": checked}


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
