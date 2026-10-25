"""Check every supplied HW9 reference against the recorded SHA-256 manifest."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXPECTED = {
    "instructions/HW9_DFE.pdf", "instructions/HW9_DFE_sample.py",
    "data/HW7_stripline_50cm.s2p", "instructions/snpreader.py",
    "instructions/cautions.txt",
}


def verify(manifest=None):
    manifest = Path(manifest) if manifest else HERE / "instructions" / "manifest.json"
    entries = json.loads(manifest.read_text(encoding="utf-8"))["files"]
    names = [entry["path"] for entry in entries]
    if len(names) != len(set(names)) or set(names) != EXPECTED:
        raise AssertionError("Manifest must list the five supplied HW9 references exactly once")
    checked = []
    for entry in entries:
        relative = Path(entry["path"])
        target = (HERE / relative).resolve()
        if relative.is_absolute() or not target.is_relative_to(HERE.resolve()):
            raise AssertionError("Manifest contains a path outside the HW9 directory")
        content = target.read_bytes()
        actual = hashlib.sha256(content).hexdigest()
        if actual != entry["sha256"]:
            raise AssertionError(f"SHA-256 mismatch: {entry['path']}")
        if "bytes" in entry and len(content) != entry["bytes"]:
            raise AssertionError(f"Byte-count mismatch: {entry['path']}")
        checked.append({"file": entry["path"], "sha256": actual, "bytes": len(content)})
    return {"passed": True, "verified_files": checked}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path)
    args = parser.parse_args(argv)
    print(json.dumps(verify(args.manifest), indent=2))


if __name__ == "__main__":
    main()
