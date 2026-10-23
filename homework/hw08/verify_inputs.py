"""Verify the supplied HW8 assignment, template, helper, channel data, and note."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def digest_stream(stream) -> tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    for chunk in iter(lambda: stream.read(1024 * 1024), b""):
        digest.update(chunk)
        size += len(chunk)
    return digest.hexdigest(), size


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir", type=Path,
        default=Path(__file__).resolve().parent / "instructions",
        help="Directory containing manifest.json and the supplied references.",
    )
    args = parser.parse_args()
    manifest = json.loads((args.data_dir / "manifest.json").read_text(encoding="utf-8"))
    for entry in manifest["files"]:
        path = args.data_dir / entry["path"]
        with path.open("rb") as stream:
            actual_hash, actual_size = digest_stream(stream)
        if (actual_hash, actual_size) != (entry["sha256"], entry["bytes"]):
            raise SystemExit(f"Integrity check failed: {entry['path']}")
        print(f"OK {entry['path']}")
    print("All supplied input files match the recorded checksums.")


if __name__ == "__main__":
    main()
