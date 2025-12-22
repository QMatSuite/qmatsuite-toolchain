# Quantum ESPRESSO Test Harness

A version-aware Python test driver for running QE test-suite tests across different platforms and QE versions.

## Overview

This test harness provides a clean, Python-based interface to the Quantum ESPRESSO test-suite. It:

- **Knows where QE is built** (per platform / toolchain)
- **Knows how to locate** the official test-suite for a given QE version
- **Runs tests** (optionally filtered) and summarizes results
- **Works from both CI and local development**

## Relationship to QE's Official Test-Suite

This harness is a wrapper around QE's official test-suite. It:

- Uses the same test inputs and reference outputs from QE's `test-suite` directory
- Runs the same QE executables (`pw.x`, `ph.x`, etc.)
- Provides a simpler, version-aware interface for CI and automation

The test-suite itself is maintained by the QE project and should be obtained from the QE source tree or downloaded separately.

## Usage

### Local Development

```bash
# List available tests
python scripts/qe_run_tests.py \
  --qe-version qe-7.5 \
  --build-dir "$RUNNER_TEMP/qe-install" \
  --list

# Run a specific test group
python scripts/qe_run_tests.py \
  --qe-version qe-7.5 \
  --build-dir "$RUNNER_TEMP/qe-install" \
  --group pw

# Run tests matching a pattern
python scripts/qe_run_tests.py \
  --qe-version qe-7.5 \
  --build-dir "$RUNNER_TEMP/qe-install" \
  --pattern "graphene*"

# Run with verbose output
python scripts/qe_run_tests.py \
  --qe-version qe-7.5 \
  --build-dir "$RUNNER_TEMP/qe-install" \
  --group pw \
  --verbose
```

### CI Integration

```yaml
# Example GitHub Actions step
- name: Run QE Tests
  run: |
    python scripts/qe_run_tests.py \
      --qe-version qe-7.5 \
      --build-dir "${{ env.QE_BUILD_DIR }}" \
      --group pw \
      --group ph
```

### Auto-Detection

The harness can auto-detect build directories and test-suite locations based on the toolchain structure:

```
toolchains/quantum_espresso/
  linux/
    upstream/q-e-qe-7.5/bin/     # Auto-detected build dir
  windows/
    mingw/upstream/q-e-qe-7.5/bin/
    oneapi/upstream/qe/bin/
  tests/
    test-suite/                  # Auto-detected test-suite
```

If auto-detection fails, specify paths explicitly:

```bash
python scripts/qe_run_tests.py \
  --qe-version qe-7.5 \
  --build-dir /path/to/qe/bin \
  --test-root /path/to/test-suite
```

## Configuration

Test groups and executable mappings are defined in per-version YAML config files:

- `configs/qe-7.5.yaml` - Configuration for QE 7.5
- `configs/qe-7.6.yaml` - Configuration for QE 7.6 (future)

Each config file defines:
- **groups**: Test groups (e.g., `pw`, `ph`) and their locations
- **patterns**: File patterns to match test inputs
- **executable_map**: Mapping from test groups to QE executables

## Test Discovery

Tests are discovered by:
1. Loading the per-version config file
2. Scanning test-suite directories for input files (`.in` files)
3. Matching against configured patterns

Test cases are represented as `QETestCase` objects with:
- `name`: Unique test identifier
- `group`: Test group (e.g., "pw", "ph")
- `workdir`: Directory where test runs
- `input_file`: Primary input file
- `reference_dir`: Where reference outputs live

## Test Execution

Each test is executed by:
1. Building command line with appropriate QE executable
2. Running in test's working directory
3. Checking for CRASH files
4. Comparing return codes (detailed diff logic can be added later)

Results are returned as `QETestResult` objects with:
- `status`: "PASS" | "FAIL" | "SKIPPED" | "ERROR"
- `details`: Short reason or error message
- `return_code`: Process return code
- `output_file`: Path to stdout
- `error_file`: Path to stderr

## Platform Support

The harness works on:
- **Linux**: Standard executable handling
- **macOS**: Standard executable handling
- **Windows**: Handles `.exe` extensions and DLL dependencies

Platform is auto-detected but can be overridden with `--platform`.

## Dependencies

- Python 3.7+
- PyYAML (`pip install pyyaml`)

## Architecture

```
toolchains/quantum_espresso/testing/
  __init__.py              # Package initialization
  qe_test_driver.py        # Core test driver (QETestSuite, QETestCase, etc.)
  qe_tests_cli.py          # CLI interface
  configs/
    qe-7.5.yaml            # Per-version configs
  README.md                # This file

scripts/
  qe_run_tests.py          # Thin wrapper for CI/human use
```

## Future Enhancements

- [ ] Parallel test execution (`--jobs` support)
- [ ] Detailed output comparison (beyond return codes)
- [ ] Automatic test-suite download
- [ ] Test result caching
- [ ] JUnit XML output for CI
- [ ] More sophisticated test filtering

## Troubleshooting

### "Config file not found"
Create a config file in `configs/` for your QE version, or use an existing version.

### "Test-suite not found"
Ensure the test-suite exists at one of:
- `toolchains/quantum_espresso/tests/test-suite/`
- `toolchains/quantum_espresso/<platform>/upstream/q-e-<version>/test-suite/`
- Or specify `--test-root` explicitly

### "Executable not found"
Verify `--build-dir` points to a directory containing QE binaries (`pw.x`, `ph.x`, etc.), or that `bin/` subdirectory exists.

### DLL errors on Windows
The harness automatically adds the build directory to PATH. If issues persist, ensure all required DLLs are in the build directory.



