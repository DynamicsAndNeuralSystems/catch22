#!/usr/bin/env python3
"""Reference test: check that the C implementation matches the MATLAB
implementation (Matlab/) for all 22 features.

tests/reference/matlab_reference.csv holds feature values computed by
tests/reference/make_matlab_reference.m on the series in
tests/reference/series/ (sines, AR processes, noise, random walks,
integer-valued series with ties, a logistic map, and a short series).
Integer-valued features must match exactly; real-valued features to a
relative tolerance (C and MATLAB agree to ~1e-12 in practice).

Usage: python3 tests/check_matlab_reference.py [path/to/run_features] [--tol 1e-9]
Uses only the Python standard library.
"""
import argparse
import csv
import os
import sys
import tempfile
from collections import defaultdict

from check_regression import ROOT, matches, read_output, run_binary


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('binary', nargs='?', default=os.path.join(ROOT, 'C', 'run_features'))
    parser.add_argument('--tol', type=float, default=1e-9,
                        help='relative tolerance for real-valued features (default 1e-9)')
    args = parser.parse_args()

    ref_dir = os.path.join(ROOT, 'tests', 'reference')
    reference = defaultdict(dict)
    with open(os.path.join(ref_dir, 'matlab_reference.csv')) as f:
        for row in csv.DictReader(f):
            reference[row['series']][row['feature']] = float(row['value'])

    failures = []
    n_values = 0
    with tempfile.TemporaryDirectory() as tmp:
        for series, expected in sorted(reference.items()):
            actual_file = os.path.join(tmp, 'output.txt')
            run_binary(args.binary, os.path.join(ref_dir, 'series', series + '.txt'), actual_file)
            actual = read_output(actual_file)
            for name, exp_value in expected.items():
                n_values += 1
                if name not in actual:
                    failures.append(f'{series}: {name} missing from output')
                elif not matches(name, actual[name], exp_value, args.tol):
                    failures.append(f'{series}: {name}: C = {actual[name]!r}, MATLAB = {exp_value!r}')

    print(f'Checked {n_values} feature values across {len(reference)} reference series against MATLAB.')
    if failures:
        print(f'{len(failures)} mismatches:')
        for failure in failures:
            print('  ' + failure)
        sys.exit(1)
    print('C matches MATLAB for all features.')


if __name__ == '__main__':
    main()
