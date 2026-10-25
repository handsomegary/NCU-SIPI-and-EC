"""Verify the supplied HW14 data and privacy-preserving Touchstone copy."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def verify() -> dict:
    folder = ROOT / "instructions"
    manifest = json.loads((folder / "input_manifest.json").read_text(encoding="utf-8"))
    channel = (folder / manifest["channel"]["file"]).read_bytes()
    assignment = (folder / manifest["assignment"]["file"]).read_bytes()
    if sha256(channel) != manifest["channel"]["published_sha256"]:
        raise ValueError("Published channel checksum differs from the recorded input")
    if sha256(assignment) != manifest["assignment"]["sha256"]:
        raise ValueError("Assignment checksum differs from the supplied original")
    lines = channel.decode("ascii").splitlines()
    options = [line for line in lines if line.lstrip().startswith("#")]
    if len(options) != 1 or options[0].split() != ["#", "GHz", "S", "MA", "R", "50.000000"]:
        raise ValueError("Unexpected Touchstone frequency/representation/reference units")
    tokens = " ".join(
        line.split("!", 1)[0].strip() for line in lines
        if not line.lstrip().startswith("#") and line.split("!", 1)[0].strip()
    ).split()
    if len(tokens) != manifest["channel"]["numeric_token_count"]:
        raise ValueError("Channel numerical token count changed")
    if sha256(" ".join(tokens).encode("ascii")) != manifest["channel"]["numeric_token_sha256"]:
        raise ValueError("Channel numerical tokens differ from the original")
    if len(tokens) % 33:
        raise ValueError("Incomplete four-port Touchstone record")
    frequencies = [float(tokens[i]) for i in range(0, len(tokens), 33)]
    if len(frequencies) != 8001 or any(abs(f - i * 0.01) > 1e-10 for i, f in enumerate(frequencies)):
        raise ValueError("Expected 8001 uniformly spaced frequencies from 0 to 80 GHz")
    return {
        "passed": True,
        "frequency_samples": len(frequencies),
        "frequency_range_ghz": [frequencies[0], frequencies[-1]],
        "frequency_step_mhz": 10,
        "numerical_tokens_unchanged": True,
        "assignment_byte_identical": True,
        "published_channel_sha256": sha256(channel),
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
