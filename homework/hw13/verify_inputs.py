"""Verify staged HW13 input bytes and schemas using only the standard library."""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
CASES = ('case1_no_tx_jitter', 'case2_sj', 'case3_rj', 'case4_sj_rj')
SIDES = ('tx', 'rx')
LEGACY_COLUMNS = ('jtx1', 'jtx2', 'jtx3', 'jtx4',
                  'jty1', 'jty2', 'jty3', 'jty4', 'tjp')
VERIFIED_COLUMNS = ('case', 'side', 'time_seconds', 'pmf', 'pdf_per_second')


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path: Path, expected_columns: tuple[str, ...]) -> list[dict[str, str]]:
    with path.open(newline='', encoding='utf-8-sig') as stream:
        reader = csv.DictReader(stream)
        require(tuple(reader.fieldnames or ()) == expected_columns,
                f'Unexpected column names: {path.name}')
        return list(reader)


def verify() -> dict[str, object]:
    manifest = json.loads((ROOT / 'instructions/manifest.json').read_text(encoding='utf-8'))
    require(manifest['schema_version'] == 1, 'Unsupported manifest schema.')
    expected_files = {
        'instructions/HW13_BER_extrapolation.pdf',
        'instructions/HW13_BER_extrapolation_sample.py',
        'inputs/jitter_variables_for_hw13.csv',
        'inputs/hw12_verified_jitter.csv',
        'inputs/hw12_source_summary.json',
        'reference/legacy_archived_analysis.py.txt',
        'reference/legacy_adjacent_analysis.py.txt',
        'reference/README.md',
    }
    require(set(manifest['files']) == expected_files, 'Unexpected manifest file set.')
    for relative, metadata in manifest['files'].items():
        path = (ROOT / relative).resolve()
        require(path.is_relative_to(ROOT), 'Manifest path escaped the HW13 directory.')
        require(path.is_file(), f'Missing input: {relative}')
        require(path.stat().st_size == metadata['bytes'], f'Size mismatch: {relative}')
        require(sha256(path) == metadata['sha256'], f'Hash mismatch: {relative}')

    legacy = read_csv(ROOT / 'inputs/jitter_variables_for_hw13.csv', LEGACY_COLUMNS)
    require(len(legacy) == 200, 'Legacy CSV must contain 200 rows.')
    for index, row in enumerate(legacy):
        values = [float(row[column]) for column in LEGACY_COLUMNS]
        require(all(math.isfinite(value) for value in values), 'Non-finite legacy value.')
        require(all(value >= 0 for value in values[:-1]), 'Negative legacy histogram.')
        require(math.isclose(values[-1], -12.5 + index * 25 / 199, abs_tol=1e-12),
                'Unexpected legacy coordinate grid.')
    for column in LEGACY_COLUMNS[:-1]:
        require(math.fsum(float(row[column]) for row in legacy) > 0,
                f'Empty legacy distribution: {column}')

    summary = json.loads((ROOT / 'inputs/hw12_source_summary.json').read_text(encoding='utf-8'))
    require(summary['source']['sha256'] == sha256(ROOT / 'inputs/hw12_verified_jitter.csv'),
            'Source summary does not identify the copied verified input.')
    require(summary['parameters']['analyzed_bits'] == 1048576, 'Wrong HW12 bit count.')
    require(summary['parameters']['points_per_ui'] == 200, 'Wrong samples per UI.')
    require(all(summary['checks'].values()), 'HW12 source checks did not all pass.')
    dt = float(summary['parameters']['dt_seconds'])
    require(math.isclose(dt, 1.25e-13, rel_tol=0, abs_tol=1e-28), 'Wrong HW12 time step.')

    verified = read_csv(ROOT / 'inputs/hw12_verified_jitter.csv', VERIFIED_COLUMNS)
    require(len(verified) == 1600, 'Verified CSV must contain 1600 rows.')
    groups = {(case, side): [] for case in CASES for side in SIDES}
    for row in verified:
        group = (row['case'], row['side'])
        require(group in groups, 'Unexpected case or side in verified input.')
        values = tuple(float(row[key]) for key in VERIFIED_COLUMNS[2:])
        require(all(math.isfinite(value) for value in values), 'Non-finite verified value.')
        require(values[1] >= 0 and values[2] >= 0, 'Negative verified probability.')
        groups[group].append(values)
    for (case, side), values in groups.items():
        require(len(values) == 200, f'Wrong bin count: {case}/{side}')
        require(math.isclose(math.fsum(value[1] for value in values), 1, abs_tol=1e-12),
                f'PMF not normalized: {case}/{side}')
        count = float(summary['cases'][case][side]['observation_count'])
        require(count > 0 and count.is_integer(), f'Invalid observation count: {case}/{side}')
        for index, (time, pmf, pdf) in enumerate(values):
            require(math.isclose(time, (index - 100) * dt, rel_tol=0, abs_tol=1e-25),
                    f'Inconsistent verified grid: {case}/{side}')
            require(math.isclose(pdf * dt, pmf, rel_tol=1e-12, abs_tol=1e-15),
                    f'PDF unit mismatch: {case}/{side}')
            require(abs(pmf * count - round(pmf * count)) <= 1e-7,
                    f'PMF/count mismatch: {case}/{side}')
        require(any(pmf == 0 for _, pmf, _ in values),
                f'Expected preserved zero cells: {case}/{side}')
    return {
        'status': 'passed',
        'files_verified': len(expected_files),
        'legacy_rows': len(legacy),
        'verified_rows': len(verified),
        'distributions': len(groups),
        'verified_source_sha256': summary['source']['sha256'],
    }


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
