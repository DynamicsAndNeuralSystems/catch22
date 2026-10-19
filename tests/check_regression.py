#!/usr/bin/env python3
"""Regression test: compare features computed by the compiled catch22 binary
against the expected outputs committed in testData/.

Each testData/<group>/<name>.txt input is run through the binary, writing to a
temporary directory (the committed <name>_output.txt files are never
modified), and every feature in the committed output must be reproduced:
exactly for integer-valued features, and to a relative tolerance otherwise
(to allow for floating-point differences across platforms and compilers).

Usage: python3 tests/check_regression.py [path/to/run_features] [--tol 1e-9]
Uses only the Python standard library.
"""
import argparse
import math
import os
import subprocess
import sys
import tempfile

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')

# Features whose output is an integer (a lag or a run length): compared exactly
INTEGER_FEATURES = {
    'CO_FirstMin_ac',
    'PD_PeriodicityWang_th0_01',
    'IN_AutoMutualInfoStats_40_gaussian_fmmi',
    'SB_BinaryStats_mean_longstretch1',
    'SB_BinaryStats_diff_longstretch0',
}


def read_output(path):
    """{feature name: value} from a run_features output file."""
    values = {}
    with open(path) as f:
        for line in f:
            parts = [p.strip() for p in line.split(',')]
            if len(parts) >= 2 and parts[1]:
                values[parts[1]] = float(parts[0])
    return values


def run_binary(binary, infile, outfile):
    # The binary asks on stdin whether to compute catch24 (1) or catch22 (0)
    result = subprocess.run([binary, infile, outfile], input='1\n',
                            capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f'{binary} exited with code {result.returncode} on {infile}')


def matches(name, got, expected, tol):
    if math.isnan(expected) or math.isnan(got):
        return math.isnan(expected) and math.isnan(got)
    if name in INTEGER_FEATURES or math.isinf(expected):
        return got == expected
    return abs(got - expected) <= tol * max(1.0, abs(expected))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('binary', nargs='?', default=os.path.join(ROOT, 'C', 'run_features'))
    parser.add_argument('--tol', type=float, default=1e-9,
                        help='relative tolerance for real-valued features (default 1e-9)')
    args = parser.parse_args()

    test_dir = os.path.join(ROOT, 'testData')
    inputs = []
    for group in sorted(os.listdir(test_dir)):
        group_dir = os.path.join(test_dir, group)
        if not os.path.isdir(group_dir):
            continue
        for name in sorted(os.listdir(group_dir)):
            expected = os.path.join(group_dir, name[:-4] + '_output.txt')
            if name.endswith('.txt') and not name.endswith('_output.txt') and os.path.exists(expected):
                inputs.append((os.path.join(group_dir, name), expected))

    failures = []
    n_values = 0
    with tempfile.TemporaryDirectory() as tmp:
        for infile, expected_file in inputs:
            actual_file = os.path.join(tmp, 'output.txt')
            run_binary(args.binary, infile, actual_file)
            expected = read_output(expected_file)
            actual = read_output(actual_file)
            for name, exp_value in expected.items():
                n_values += 1
                if name not in actual:
                    failures.append(f'{os.path.relpath(infile, ROOT)}: {name} missing from output')
                elif not matches(name, actual[name], exp_value, args.tol):
                    failures.append(f'{os.path.relpath(infile, ROOT)}: {name} = {actual[name]!r}, expected {exp_value!r}')

    print(f'Checked {n_values} feature values across {len(inputs)} time series.')
    if failures:
        print(f'{len(failures)} mismatches:')
        for failure in failures[:50]:
            print('  ' + failure)
        if len(failures) > 50:
            print(f'  ... and {len(failures) - 50} more')
        sys.exit(1)
    print('All match the expected outputs.')


if __name__ == '__main__':
    main()
