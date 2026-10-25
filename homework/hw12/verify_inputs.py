"""Check the complete supplied HW12 instruction files against their SHA-256 manifest."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path


def main() -> None:
    folder = Path(__file__).resolve().parent / "instructions"
    manifest = json.loads((folder / "manifest.json").read_text(encoding="utf-8"))
    for item in manifest["files"]:
        path = folder / item["path"]
        data = path.read_bytes()
        if len(data) != item["bytes"] or hashlib.sha256(data).hexdigest() != item["sha256"]:
            raise SystemExit(f"Input integrity check failed: {item['path']}")
        print(f"OK {item['path']}")
    print("All three supplied instruction files match the original bytes.")


if __name__ == "__main__":
    main()
