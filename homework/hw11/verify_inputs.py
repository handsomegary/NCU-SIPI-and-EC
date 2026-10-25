"""Verify all supplied HW11 inputs against their recorded SHA-256 digests."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


def main() -> None:
    root = Path(__file__).resolve().parent
    manifest = json.loads((root / "source_manifest.json").read_text(encoding="utf-8"))
    records = manifest["files"]
    if sum(item["source_group"] == "Instruction" for item in records) != 4:
        raise ValueError("Expected all four supplied Instruction files")
    for item in records:
        path = (root / item["path"]).resolve()
        if not path.is_relative_to(root):
            raise ValueError("Manifest path must remain inside the homework package")
        with path.open("rb") as stream:
            digest = hashlib.file_digest(stream, "sha256").hexdigest()
        if digest != item["staged_sha256"]:
            raise ValueError(f"Checksum mismatch: {item['path']}")
        if path.stat().st_size != item["bytes"]:
            raise ValueError(f"Size mismatch: {item['path']}")
        if not item["modified"] and digest != item["original_sha256"]:
            raise ValueError(f"Original-byte mismatch: {item['path']}")
    print(f"Verified {len(records)} supplied files, including all four Instruction files.")


if __name__ == "__main__":
    main()
