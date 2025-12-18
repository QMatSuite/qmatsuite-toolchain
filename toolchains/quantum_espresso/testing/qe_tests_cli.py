"""
CLI interface for Quantum ESPRESSO test harness.
"""

import argparse
import sys
from pathlib import Path
from typing import Optional

from .qe_test_driver import QETestSuite, PlatformInfo


def main(argv=None):
    """
    Main CLI entry point.
    
    Usage examples:
    
      python -m toolchains.quantum_espresso.testing.qe_tests_cli \\
        --qe-version qe-7.5 \\
        --build-dir /path/to/build \\
        --group pw
    
      python -m toolchains.quantum_espresso.testing.qe_tests_cli \\
        --qe-version qe-7.5 \\
        --build-dir /path/to/build \\
        --pattern graphene*
    """
    parser = argparse.ArgumentParser(
        description='Run Quantum ESPRESSO test-suite tests',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument(
        '--qe-version',
        required=True,
        help='QE version (e.g., qe-7.5)'
    )
    
    parser.add_argument(
        '--build-dir',
        type=Path,
        help='Directory containing QE binaries (default: auto-detect from toolchain)'
    )
    
    parser.add_argument(
        '--test-root',
        type=Path,
        help='Path to test-suite root (default: auto-detect)'
    )
    
    parser.add_argument(
        '--platform',
        choices=['linux', 'macos', 'windows'],
        help='Platform override (default: auto-detect)'
    )
    
    parser.add_argument(
        '--group',
        action='append',
        dest='groups',
        help='Test group to run (can be specified multiple times)'
    )
    
    parser.add_argument(
        '--pattern',
        help='Glob pattern to match test names'
    )
    
    parser.add_argument(
        '--jobs',
        type=int,
        default=1,
        help='Number of parallel jobs (currently ignored, runs serially)'
    )
    
    parser.add_argument(
        '--limit',
        type=int,
        help='Limit number of tests to run (useful for quick testing)'
    )
    
    parser.add_argument(
        '--list',
        action='store_true',
        help='List available tests without running them'
    )
    
    parser.add_argument(
        '--verbose',
        '-v',
        action='store_true',
        help='Verbose output'
    )
    
    args = parser.parse_args(argv)
    
    # Auto-detect build directory if not provided
    build_dir = args.build_dir
    if not build_dir:
        # Try to infer from toolchain structure
        # Create a temporary suite just for detection
        temp_suite = QETestSuite(
            qe_version=args.qe_version,
            build_dir=Path.cwd(),  # Temporary, will be updated
            test_root=args.test_root,
            platform=args.platform
        )
        inferred_build = temp_suite._find_toolchain_build_dir()
        if inferred_build:
            build_dir = inferred_build
            if args.verbose:
                print(f"Auto-detected build directory: {build_dir}")
        else:
            print("ERROR: --build-dir not specified and could not auto-detect")
            print("Please specify --build-dir explicitly")
            print("\nExpected locations:")
            print(f"  - toolchains/quantum_espresso/{args.platform or PlatformInfo.detect()}/upstream/q-e-{args.qe_version}/bin")
            print(f"  - toolchains/quantum_espresso/{args.platform or PlatformInfo.detect()}/upstream/qe/bin")
            if (args.platform or PlatformInfo.detect()) == 'windows':
                print(f"  - toolchains/quantum_espresso/windows/mingw/dist/win-mingw")
                print(f"  - toolchains/quantum_espresso/windows/oneapi/dist/win-oneapi")
            return 1
    
    # Create test suite
    try:
        suite = QETestSuite(
            qe_version=args.qe_version,
            build_dir=build_dir,
            test_root=args.test_root,
            platform=args.platform
        )
    except Exception as e:
        print(f"ERROR: Failed to initialize test suite: {e}")
        return 1
    
    # List tests if requested
    if args.list:
        try:
            tests = suite.discover_tests(groups=args.groups)
            print(f"\nFound {len(tests)} test(s):\n")
            for test in tests:
                print(f"  {test.name:40s} [{test.group:8s}] {test.input_file}")
            return 0
        except Exception as e:
            print(f"ERROR: Failed to discover tests: {e}")
            return 1
    
    # Run tests
    try:
        if args.verbose:
            print(f"QE Version: {args.qe_version}")
            print(f"Build Dir:  {build_dir}")
            print(f"Platform:   {suite.platform}")
            if args.groups:
                print(f"Groups:     {', '.join(args.groups)}")
            if args.pattern:
                print(f"Pattern:    {args.pattern}")
            print()
        
        results = suite.run_tests(
            pattern=args.pattern,
            groups=args.groups,
            jobs=args.jobs,
            limit=args.limit
        )
        
        # Print summary
        print(f"\n{'='*60}")
        print(f"Test Results: {len(results)} test(s)")
        print(f"{'='*60}\n")
        
        passed = sum(1 for r in results if r.status == "PASS")
        failed = sum(1 for r in results if r.status == "FAIL")
        errors = sum(1 for r in results if r.status == "ERROR")
        skipped = sum(1 for r in results if r.status == "SKIPPED")
        
        for result in results:
            status_symbol = {
                "PASS": "[PASS]",
                "FAIL": "[FAIL]",
                "ERROR": "[ERROR]",
                "SKIPPED": "[SKIP]"
            }.get(result.status, "[?]")
            
            print(f"{status_symbol} {result.case.name:40s} [{result.case.group:8s}] {result.status}")
            if result.details and args.verbose:
                print(f"    {result.details}")
        
        print(f"\n{'='*60}")
        print(f"Summary: {passed} passed, {failed} failed, {errors} errors, {skipped} skipped")
        print(f"{'='*60}\n")
        
        return 0 if (failed == 0 and errors == 0) else 1
        
    except Exception as e:
        print(f"ERROR: Failed to run tests: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())

