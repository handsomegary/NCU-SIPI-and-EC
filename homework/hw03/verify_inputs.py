"""Verify the supplied HW3 files, including the lossless ADS gzip round trip."""
from __future__ import annotations

import argparse
import gzip
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
    parser.add_argument("--data-dir", type=Path, default=Path(__file__).resolve().parent / "instructions")
    parser.add_argument("--quick", action="store_true", help="Check stored hashes without decompressing ADS.")
    args = parser.parse_args()
    manifest = json.loads((args.data_dir / "manifest.json").read_text(encoding="utf-8"))
    for entry in manifest["files"]:
        path = args.data_dir / entry["path"]
        compressed = entry["storage"] == "lossless-gzip"
        expected_hash = entry["stored_sha256"] if compressed else entry["sha256"]
        expected_size = entry["stored_bytes"] if compressed else entry["bytes"]
        with path.open("rb") as stream:
            actual_hash, actual_size = digest_stream(stream)
        if (actual_hash, actual_size) != (expected_hash, expected_size):
            raise SystemExit(f"Integrity check failed: {entry['path']}")
        if compressed and not args.quick:
            with gzip.open(path, "rb") as stream:
                restored_hash, restored_size = digest_stream(stream)
            if (restored_hash, restored_size) != (entry["sha256"], entry["bytes"]):
                raise SystemExit(f"Decompressed integrity check failed: {entry['original_name']}")
        print(f"OK {entry['path']}")
    print("All supplied input files match the recorded checksums.")


if __name__ == "__main__":
    main()
