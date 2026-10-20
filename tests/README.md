# Tests

CI (`.github/workflows/ci.yml`) runs three tests on every platform and compiler the wrappers build on (Linux x86-64 and arm64 with GCC and clang, macOS arm64 and x86-64, Windows x86-64 and arm64 with MSVC and MinGW), with the flags they use (`-std=c99` as pycatch22, `-std=c11 -fPIC` with GCC 9 as catch22_jll) and at `-O0`, `-O3` and 32-bit; then again under AddressSanitizer and UndefinedBehaviorSanitizer (which also catch memory leaks) and under valgrind (uninitialised reads). `tests/build.sh` is the one compile line for all of them (`CC` and `CFLAGS` from the environment, MSVC's `cl` included). Run the tests locally from the repository root:

```bash
gcc -Wall -Wextra -Werror -O2 C/*.c -o C/run_features -lm

# 1. Analytic unit tests: constructed series with known answers
gcc -O2 -IC tests/test_analytic.c $(ls C/*.c | grep -v main.c) -o test_analytic -lm
./test_analytic

# 2. Regression test: every series in testData/ against its committed *_output.txt
python3 tests/check_regression.py C/run_features

# 3. Reference test: the C matches hctsa (each feature's hctsa call)
python3 tests/check_hctsa_reference.py C/run_features
```

The Python scripts use only the standard library.

## What each test checks

- **`test_analytic.c`**: feature values on time series whose correct output follows from first principles. For example, a sine wave of period _P_ has `PD_PeriodicityWang_th0_01` = _P_ and `CO_FirstMin_ac` = ⌊_P_/2⌋; a ramp has `SB_BinaryStats_diff_longstretch0` = 0; white noise has `SB_MotifThree_quantile_hh` ≈ log 9.
- **`check_regression.py`**: catches any unintended change in feature values. Integer-valued features must match exactly; real-valued features to a relative tolerance of 1e-9 (to allow for floating-point differences across platforms and compilers). The committed expected outputs are never modified.
- **`check_hctsa_reference.py`**: checks the C (run with `--catch24 --no-mask`) against `reference/hctsa_reference.csv`. That file holds, for each series in `reference/series/`, the value of each feature's hctsa call. The calls are listed by `run_features --list` (in `reference/features.tsv`; `C/catch22.c` is the source), and `reference/hctsa_reference_version.txt` records the hctsa commit. All 24 features agree to 1e-9.

## Feature values unchanged (pull requests)

On a pull request, CI also builds `tests/print_features.c` against the base branch and the head and diffs the two outputs (every feature of every series in `testData/`, with and without masking, as hex floats: bit for bit), on each platform and compiler. The committed outputs are generated on Linux and agree with other platforms only to the regression tolerance, so this is the check that a change is exactly value-preserving where it claims to be. A pull request that changes values on purpose fails it, and says so.

## Updating expected values after an intentional change

If a change is meant to alter feature values:

1. Regenerate the expected outputs: `bash testData/runtests.sh` (this overwrites the `*_output.txt` files in `testData/`). Check with `git diff` that only the intended features changed. Regenerating on a different platform (e.g., macOS rather than the Linux CI) can change the last digits of real-valued features; keep only the intended changes.
   The committed outputs (catch24, masked) were generated on Linux (GCC 11).
2. If a feature's hctsa call changed (or for a new hctsa release), regenerate the reference values: `C/run_features --list > tests/reference/features.tsv`, then in MATLAB, from `tests/reference`, `make_hctsa_reference('/path/to/hctsa')`.
3. Update the expected values in `test_analytic.c` if they are affected.
4. Note the change in the README, since results are then not comparable across versions.
