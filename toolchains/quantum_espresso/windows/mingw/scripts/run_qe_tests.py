#!/usr/bin/env python3
"""
Cross-platform script to run Quantum ESPRESSO test-suite.

This script automates the setup and execution of QE test-suite tests,
working on Windows, Linux, and macOS.

Usage:
    python run_qe_tests.py [options] [category]

Examples:
    python run_qe_tests.py --category pw_all
    python run_qe_tests.py --category pw_scf --nprocs 2
    python run_qe_tests.py --category pw_all cp_all --verbose
"""

import argparse
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path


def find_qe_root(start_path=None):
    """Find QE root directory by looking for test-suite directory."""
    if start_path is None:
        start_path = Path(__file__).resolve().parent
    
    # Walk up from script location
    current = Path(start_path).resolve()
    for _ in range(10):  # Limit search depth
        test_suite = current / "test-suite"
        if test_suite.exists() and (test_suite / "testcode").exists():
            return current
        if current.parent == current:
            break
        current = current.parent
    
    # Try common locations
    common_paths = [
        Path.cwd(),
        Path.cwd() / "upstream" / "qe",
        Path.cwd() / "upstream" / "q-e-qe-7.5",
    ]
    
    for path in common_paths:
        test_suite = path / "test-suite"
        if test_suite.exists() and (test_suite / "testcode").exists():
            return path.resolve()
    
    return None


def setup_environment(qe_root, qe_build=None, pseudo_dir=None, tmp_dir=None):
    """Set up environment variables for test-suite."""
    qe_root = Path(qe_root).resolve()
    
    # Default paths
    if qe_build is None:
        qe_build = qe_root
    else:
        qe_build = Path(qe_build).resolve()
    
    if pseudo_dir is None:
        pseudo_dir = qe_root / "pseudo"
    else:
        pseudo_dir = Path(pseudo_dir).resolve()
    
    if tmp_dir is None:
        # Use platform-appropriate temp directory
        if platform.system() == "Windows":
            tmp_dir = Path(os.environ.get("TEMP", os.environ.get("TMP", "C:\\temp")))
        else:
            tmp_dir = Path("/tmp")
        tmp_dir = tmp_dir / "qe-test-suite"
    else:
        tmp_dir = Path(tmp_dir).resolve()
    
    # Create temp directory
    tmp_dir.mkdir(parents=True, exist_ok=True)
    
    env = os.environ.copy()
    
    # Set QE environment variables
    env["ESPRESSO_ROOT"] = str(qe_root)
    env["ESPRESSO_BUILD"] = str(qe_build)
    env["ESPRESSO_PSEUDO"] = str(pseudo_dir)
    env["ESPRESSO_TMPDIR"] = str(tmp_dir)
    env["TESTCODE_DIR"] = str(qe_root / "test-suite" / "testcode")
    env["NETWORK_PSEUDO"] = "https://pseudopotentials.quantum-espresso.org/upf_files/"
    
    # Set OpenMP threads (default to 1)
    env.setdefault("OMP_NUM_THREADS", "1")
    
    return env, tmp_dir


def generate_userconfig(qe_root, test_suite_dir):
    """Generate userconfig from userconfig.tmp."""
    userconfig_tmp = test_suite_dir / "userconfig.tmp"
    userconfig = test_suite_dir / "userconfig"
    
    if not userconfig_tmp.exists():
        print(f"ERROR: userconfig.tmp not found at {userconfig_tmp}")
        return False
    
    # Read template
    with open(userconfig_tmp, "r", encoding="utf-8") as f:
        content = f.read()
    
    # Replace XXXXXX with actual QE root path
    # Handle both forward and backward slashes for cross-platform compatibility
    qe_root_str = str(qe_root).replace("\\", "/")
    content = content.replace("XXXXXX", qe_root_str)
    
    # Write userconfig
    with open(userconfig, "w", encoding="utf-8") as f:
        f.write(content)
    
    return True


def check_pseudos(test_suite_dir, env):
    """Check if pseudopotentials are available."""
    check_script = test_suite_dir / "check_pseudo.sh"
    if not check_script.exists():
        print("WARNING: check_pseudo.sh not found, skipping pseudo check")
        return True
    
    # On Windows, we need to run via bash if available
    if platform.system() == "Windows":
        # Try to find bash (MSYS2, Git Bash, WSL)
        bash_paths = [
            "C:\\tools\\msys64\\usr\\bin\\bash.exe",
            "C:\\msys64\\usr\\bin\\bash.exe",
            "bash.exe",
        ]
        bash = None
        for path in bash_paths:
            if shutil.which(path) or Path(path).exists():
                bash = path
                break
        
        if bash:
            cmd = [bash, str(check_script), "pw_"]
        else:
            print("WARNING: bash not found, skipping pseudo check")
            return True
    else:
        cmd = ["bash", str(check_script), "pw_"]
    
    try:
        result = subprocess.run(cmd, env=env, cwd=test_suite_dir, 
                              capture_output=True, text=True)
        return result.returncode == 0
    except Exception as e:
        print(f"WARNING: Failed to check pseudos: {e}")
        return True  # Continue anyway


def run_tests(test_suite_dir, category, env, nprocs=None, nbgrp=None, verbose=1):
    """Run test-suite tests using testcode.py."""
    testcode_py = test_suite_dir / "testcode" / "bin" / "testcode.py"
    
    if not testcode_py.exists():
        print(f"ERROR: testcode.py not found at {testcode_py}")
        return False
    
    # Set MPI environment variables if specified
    if nprocs:
        env["QE_USE_MPI"] = str(nprocs)
        if nbgrp:
            env["QE_USE_BGRP"] = str(nbgrp)
    else:
        env.pop("QE_USE_MPI", None)
        env.pop("QE_USE_BGRP", None)
    
    # Build command
    cmd = [sys.executable, str(testcode_py), "--verbose"]
    
    if verbose > 1:
        cmd.extend(["-" + "v" * (verbose - 1)])
    
    cmd.extend(["--category", category])
    
    print(f"Running tests: {' '.join(cmd)}")
    if nprocs:
        print(f"  Using {nprocs} MPI processes")
    if nbgrp:
        print(f"  Using {nbgrp} band groups")
    print()
    
    # Run testcode
    try:
        result = subprocess.run(cmd, env=env, cwd=test_suite_dir)
        return result.returncode == 0
    except KeyboardInterrupt:
        print("\nInterrupted by user")
        return False
    except Exception as e:
        print(f"ERROR: Failed to run tests: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Run Quantum ESPRESSO test-suite tests",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run all PW tests
  python run_qe_tests.py --category pw_all

  # Run specific test category with MPI
  python run_qe_tests.py --category pw_scf --nprocs 4

  # Run multiple categories
  python run_qe_tests.py --category pw_all cp_all

  # Specify custom QE root and build directories
  python run_qe_tests.py --qe-root /path/to/qe --qe-build /path/to/qe/build --category pw_all
        """
    )
    
    parser.add_argument(
        "--category", "-c",
        nargs="+",
        default=["pw_all"],
        help="Test category to run (e.g., pw_all, cp_all, ph_all). Can specify multiple. Default: pw_all"
    )
    
    parser.add_argument(
        "--qe-root",
        type=Path,
        help="QE root directory (auto-detected if not specified)"
    )
    
    parser.add_argument(
        "--qe-build",
        type=Path,
        help="QE build directory (defaults to qe-root if not specified)"
    )
    
    parser.add_argument(
        "--pseudo-dir",
        type=Path,
        help="Pseudopotential directory (defaults to qe-root/pseudo)"
    )
    
    parser.add_argument(
        "--tmp-dir",
        type=Path,
        help="Temporary directory for test outputs"
    )
    
    parser.add_argument(
        "--nprocs", "-n",
        type=int,
        help="Number of MPI processes to use"
    )
    
    parser.add_argument(
        "--nbgrp",
        type=int,
        help="Number of band groups (requires --nprocs)"
    )
    
    parser.add_argument(
        "--verbose", "-v",
        action="count",
        default=1,
        help="Increase verbosity (can be specified multiple times)"
    )
    
    parser.add_argument(
        "--skip-pseudo-check",
        action="store_true",
        help="Skip pseudopotential availability check"
    )
    
    args = parser.parse_args()
    
    # Find QE root
    if args.qe_root:
        qe_root = Path(args.qe_root).resolve()
    else:
        qe_root = find_qe_root()
        if qe_root is None:
            print("ERROR: Could not find QE root directory.")
            print("Please specify --qe-root or run from QE directory.")
            return 1
        print(f"Found QE root: {qe_root}")
    
    test_suite_dir = qe_root / "test-suite"
    if not test_suite_dir.exists():
        print(f"ERROR: test-suite directory not found at {test_suite_dir}")
        return 1
    
    # Setup environment
    env, tmp_dir = setup_environment(
        qe_root, args.qe_build, args.pseudo_dir, args.tmp_dir
    )
    
    print(f"QE Root: {qe_root}")
    print(f"QE Build: {env['ESPRESSO_BUILD']}")
    print(f"Pseudo Dir: {env['ESPRESSO_PSEUDO']}")
    print(f"Temp Dir: {tmp_dir}")
    print()
    
    # Generate userconfig
    if not generate_userconfig(qe_root, test_suite_dir):
        return 1
    
    # Check pseudos
    if not args.skip_pseudo_check:
        print("Checking pseudopotentials...")
        check_pseudos(test_suite_dir, env)
        print()
    
    # Run tests for each category
    success = True
    for category in args.category:
        print(f"{'='*60}")
        print(f"Running category: {category}")
        print(f"{'='*60}")
        
        if not run_tests(test_suite_dir, category, env, args.nprocs, args.nbgrp, args.verbose):
            success = False
            print(f"\nERROR: Tests failed for category {category}")
        else:
            print(f"\nSUCCESS: Tests completed for category {category}")
        print()
    
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())



