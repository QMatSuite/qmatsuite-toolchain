# Migration Notes

## What Was Created

This refactor creates a new, clean test harness at `toolchains/quantum_espresso/testing/` that provides a version-aware, Python-based interface to the QE test-suite.

## Relationship to Existing Code

### Existing Test-Suite Code (tests/test-suite/)

The existing Python code in `toolchains/quantum_espresso/tests/test-suite/` serves a different purpose:
- `qe_common.py`, `qe_runner.py`, `qe_extractor.py` - Used by QE's testcode framework
- These are low-level utilities that replace shell scripts within the test-suite
- They are called by testcode.py to run individual tests

### New Test Harness (testing/)

The new test harness is a higher-level abstraction:
- Provides version-aware test discovery and execution
- Works with per-version YAML configs
- Integrates with toolchain build locations
- Designed for CI and automation

These two layers can coexist:
- The test-suite's internal Python scripts handle individual test execution
- The new harness orchestrates running groups of tests and provides a clean API

## What to Keep/Remove

### Keep
- `toolchains/quantum_espresso/tests/test-suite/` - The actual test-suite with all its Python utilities
- `toolchains/quantum_espresso/testing/` - The new test harness

### Can Remove (if desired)
- `toolchains/quantum_espresso/tests/scripts/run_qe_tests.py` - Old test runner, replaced by new harness
- `toolchains/quantum_espresso/windows/mingw/scripts/run_qe_tests.py` - Platform-specific old runner

However, these old scripts might still be useful for backward compatibility, so consider keeping them temporarily.

## Usage Migration

### Old Way
```bash
python toolchains/quantum_espresso/tests/scripts/run_qe_tests.py --category pw_all
```

### New Way
```bash
python scripts/qe_run_tests.py --qe-version qe-7.5 --build-dir /path/to/build --group pw
```

## Benefits of New Structure

1. **Version-aware**: Easy to add support for new QE versions (just add a YAML config)
2. **Config-driven**: Test groups and mappings in YAML, not hard-coded
3. **Toolchain-aware**: Auto-detects build directories based on platform
4. **Clean API**: Class-based design with clear separation of concerns
5. **CI-friendly**: Simple CLI interface perfect for automation

## Next Steps

1. Update CI workflows to use `scripts/qe_run_tests.py`
2. Add more test groups to `configs/qe-7.5.yaml` as needed
3. Create configs for other QE versions (qe-7.6, etc.) when needed
4. Consider adding parallel execution support
5. Consider adding detailed output comparison

