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
import tarfile
import urllib.request
from pathlib import Path


def download_test_suite(test_suite_dir, version="qe-7.5"):
    """Download test-suite from GitLab if it doesn't exist."""
    if test_suite_dir.exists() and (test_suite_dir / "testcode").exists():
        print(f"Test-suite already exists at {test_suite_dir}")
        print(f"DEBUG: Skipping download")
        return True
    
    # If directory exists but is incomplete, try to remove it
    if test_suite_dir.exists():
        print(f"DEBUG: Found incomplete test-suite directory, attempting to remove...")
        import time
        for attempt in range(5):
            try:
                shutil.rmtree(test_suite_dir)
                print(f"DEBUG: Successfully removed incomplete directory")
                break
            except (PermissionError, OSError) as e:
                if attempt < 4:
                    print(f"DEBUG: Retry {attempt + 1}/5: waiting 2 seconds...")
                    time.sleep(2)
                else:
                    print(f"ERROR: Could not remove test-suite directory after 5 attempts: {e}")
                    print(f"  Please manually close any programs using files in {test_suite_dir}")
                    print(f"  Or delete the directory manually and run again")
                    return False
    
    print(f"Downloading test-suite version {version}...")
    
    # GitLab URL for test-suite archive
    # For version 7.5, we need to get the test-suite from the qe-7.5 tag
    url = f"https://gitlab.com/QEF/q-e/-/archive/{version}/q-e-{version}.tar.gz"
    
    # Create parent directory
    test_suite_dir.parent.mkdir(parents=True, exist_ok=True)
    
    # Download to temporary file
    temp_tarball = test_suite_dir.parent / f"q-e-{version}.tar.gz"
    
    try:
        print(f"Downloading from {url}...")
        urllib.request.urlretrieve(url, temp_tarball)
        print("Download complete. Extracting...")
        
        # Extract tarball
        # The tarball structure is: q-e-{version}/test-suite/...
        extract_dir = test_suite_dir.parent / f"q-e-{version}-temp"
        
        print(f"Extracting to temporary directory: {extract_dir}")
        with tarfile.open(temp_tarball, "r:gz") as tar:
            # Extract everything to a temporary location
            tar.extractall(extract_dir.parent)
        
        # Find test-suite in extracted directory
        # The archive contains q-e-{version}/test-suite/...
        qe_version_dir = extract_dir.parent / f"q-e-{version}"
        extracted_test_suite = qe_version_dir / "test-suite"
        
        print(f"DEBUG: Looking for test-suite at: {extracted_test_suite}")
        if extracted_test_suite.exists() and (extracted_test_suite / "testcode").exists():
            print(f"DEBUG: Found test-suite, moving to {test_suite_dir}")
            # Move test-suite to target location
            if test_suite_dir.exists():
                print(f"DEBUG: Removing existing test-suite directory")
                # Try multiple times with delay in case files are locked
                import time
                for attempt in range(3):
                    try:
                        shutil.rmtree(test_suite_dir)
                        break
                    except PermissionError as e:
                        if attempt < 2:
                            print(f"DEBUG: Retry {attempt + 1}/3: waiting 1 second...")
                            time.sleep(1)
                        else:
                            print(f"ERROR: Could not remove test-suite directory: {e}")
                            print(f"  Please close any programs using files in {test_suite_dir}")
                            return False
            
            shutil.move(str(extracted_test_suite), str(test_suite_dir))
            print(f"DEBUG: Test-suite moved successfully")
            # Clean up extracted directory
            try:
                print(f"DEBUG: Cleaning up temporary extraction directory: {qe_version_dir}")
                shutil.rmtree(qe_version_dir, ignore_errors=True)
                if extract_dir.exists():
                    shutil.rmtree(extract_dir, ignore_errors=True)
            except Exception as e:
                print(f"DEBUG: Warning: Could not clean up temp directory: {e}")
        else:
            print(f"ERROR: test-suite not found in extracted archive")
            print(f"  Expected at: {extracted_test_suite}")
            # Try to list what was extracted
            if qe_version_dir.exists():
                print(f"  Extracted directory contents: {list(qe_version_dir.iterdir())[:10]}")
            if extract_dir.exists():
                print(f"  Temp directory contents: {list(extract_dir.iterdir())[:10]}")
            return False
        
        # Clean up tarball
        temp_tarball.unlink()
        
        print(f"Test-suite extracted to {test_suite_dir}")
        return True
        
    except Exception as e:
        print(f"ERROR: Failed to download test-suite: {e}")
        if temp_tarball.exists():
            temp_tarball.unlink()
        return False


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
    
    # Try common locations relative to script
    script_dir = Path(__file__).resolve().parent
    common_paths = [
        script_dir.parent.parent / "windows" / "mingw" / "upstream" / "qe",
        script_dir.parent.parent / "windows" / "mingw" / "upstream" / "q-e-qe-7.5",
        Path.cwd(),
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


def find_bash():
    """Find bash executable on Windows."""
    if platform.system() != "Windows":
        return "bash"
    
    bash_paths = [
        "C:\\tools\\msys64\\usr\\bin\\bash.exe",
        "C:\\msys64\\usr\\bin\\bash.exe",
        "bash.exe",
    ]
    
    for path in bash_paths:
        if shutil.which(path) or Path(path).exists():
            return path
    
    return None


def generate_userconfig(qe_root, test_suite_dir):
    """Generate userconfig from userconfig.tmp."""
    userconfig_tmp = test_suite_dir / "userconfig.tmp"
    userconfig = test_suite_dir / "userconfig"
    
    print(f"DEBUG: Generating userconfig from {userconfig_tmp}")
    print(f"DEBUG: Target userconfig: {userconfig}")
    
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
    
    # On Windows, create batch file wrappers for shell scripts
    if platform.system() == "Windows":
        bash = find_bash()
        if bash:
            bash_path = str(Path(bash).resolve())
            print(f"DEBUG: Found bash at: {bash_path}")
            # Create batch file wrappers for .sh scripts
            import re
            # Match lines like: exe = XXXXXX/test-suite/run-*.sh
            pattern = r'(exe\s*=\s*)([^\n]+/run-([^\n]+)\.sh)'
            def replace_exe(match):
                prefix = match.group(1)
                script_path = match.group(2)
                script_name = match.group(3)  # e.g., "pw", "cp", etc.
                # Create a batch file wrapper in the test-suite directory
                wrapper_name = f"run-{script_name}.bat"
                wrapper_path = test_suite_dir / wrapper_name
                
                # Write batch file that calls bash with the script
                wrapper_content = f"""@echo off
"{bash_path}" "{script_path}" %*
"""
                with open(wrapper_path, "w", encoding="utf-8") as f:
                    f.write(wrapper_content)
                
                # Return the wrapper path (use forward slashes for consistency)
                wrapper_path_str = str(wrapper_path).replace("\\", "/")
                return f"{prefix}{wrapper_path_str}"
            content = re.sub(pattern, replace_exe, content)
            print(f"DEBUG: Created batch file wrappers for .sh scripts")
        else:
            print("WARNING: bash not found, .sh scripts may not work on Windows")
    
    print(f"DEBUG: Replaced XXXXXX with: {qe_root_str}")
    
    # Write userconfig
    with open(userconfig, "w", encoding="utf-8") as f:
        f.write(content)
    
    print(f"DEBUG: userconfig written successfully")
    return True


def check_pseudos(test_suite_dir, env):
    """Check if pseudopotentials are available."""
    check_script = test_suite_dir / "check_pseudo.sh"
    if not check_script.exists():
        print("WARNING: check_pseudo.sh not found, skipping pseudo check")
        return True
    
    print(f"DEBUG: Checking pseudos with script: {check_script}")
    
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
                print(f"DEBUG: Found bash at: {bash}")
                break
        
        if bash:
            cmd = [bash, str(check_script), "pw_"]
        else:
            print("WARNING: bash not found, skipping pseudo check")
            return True
    else:
        cmd = ["bash", str(check_script), "pw_"]
    
    try:
        print(f"DEBUG: Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, env=env, cwd=test_suite_dir, 
                              capture_output=True, text=True)
        print(f"DEBUG: Pseudo check exit code: {result.returncode}")
        if result.stdout:
            print(f"DEBUG: Pseudo check stdout: {result.stdout[:500]}")
        if result.stderr:
            print(f"DEBUG: Pseudo check stderr: {result.stderr[:500]}")
        return result.returncode == 0
    except Exception as e:
        print(f"WARNING: Failed to check pseudos: {e}")
        import traceback
        traceback.print_exc()
        return True  # Continue anyway


def run_tests(test_suite_dir, category, env, nprocs=None, nbgrp=None, verbose=1, pw_exe=None):
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
    
    # Override PW executable if specified
    if pw_exe:
        # Resolve path - if relative, resolve relative to test_suite_dir or current working directory
        pw_exe_path = Path(pw_exe)
        if not pw_exe_path.is_absolute():
            # Try relative to test_suite_dir first
            test_suite_pw = test_suite_dir / pw_exe_path
            if test_suite_pw.exists():
                pw_exe_path = test_suite_pw.resolve()
            else:
                # Try relative to current working directory
                pw_exe_path = pw_exe_path.resolve()
        else:
            pw_exe_path = pw_exe_path.resolve()
        
        if not pw_exe_path.exists():
            print(f"ERROR: pw.x executable not found at {pw_exe_path}")
            print(f"  Tried resolving from: {Path(pw_exe)}")
            if not Path(pw_exe).is_absolute():
                print(f"  Also tried: {test_suite_dir / pw_exe_path}")
            return False
        
        # Set ESPRESSO_BUILD to point to the directory containing the executable
        # The run-pw.sh script expects ${ESPRESSO_BUILD}/bin/pw.x
        pw_bin_dir = pw_exe_path.parent
        
        # On Windows, create expected directory structure
        if platform.system() == "Windows":
            # Create a temporary bin directory structure that run-pw.sh expects
            temp_build_dir = test_suite_dir / "temp_build"
            temp_bin_dir = temp_build_dir / "bin"
            temp_bin_dir.mkdir(parents=True, exist_ok=True)
            
            print(f"DEBUG: Created temp build directory: {temp_build_dir}")
            print(f"DEBUG: Created temp bin directory: {temp_bin_dir}")
            
            # Create symlink or copy pw.x to expected location
            # On Windows, the executable might be pw.x.exe, but run-pw.sh expects pw.x
            expected_pw = temp_bin_dir / "pw.x"
            
            # Check if source and destination are the same (avoid self-copy)
            if pw_exe_path.resolve() == expected_pw.resolve():
                print(f"DEBUG: Source and destination are the same: {pw_exe_path}")
                print(f"DEBUG: Using existing file at {expected_pw}")
            else:
                # Remove existing file if it exists
                if expected_pw.exists():
                    expected_pw.unlink()
                    print(f"DEBUG: Removed existing {expected_pw}")
                
                # Try to create a symlink (works on Windows 10+ with developer mode)
                try:
                    if hasattr(os, "symlink"):
                        os.symlink(str(pw_exe_path), str(expected_pw))
                        print(f"DEBUG: Created symlink: {expected_pw} -> {pw_exe_path}")
                    else:
                        # Fallback: copy the executable
                        shutil.copy2(pw_exe_path, expected_pw)
                        print(f"DEBUG: Copied executable: {pw_exe_path} -> {expected_pw}")
                except Exception as e:
                    # Fallback: copy the executable
                    try:
                        shutil.copy2(pw_exe_path, expected_pw)
                        print(f"DEBUG: Copied executable (symlink failed: {e}): {pw_exe_path} -> {expected_pw}")
                    except Exception as e2:
                        print(f"ERROR: Failed to copy/symlink executable: {e2}")
                        print(f"  Source: {pw_exe_path}")
                        print(f"  Destination: {expected_pw}")
                        return False
            
            # Verify the file exists
            if not expected_pw.exists():
                print(f"ERROR: Failed to create {expected_pw}")
                return False
            
            env["ESPRESSO_BUILD"] = str(temp_build_dir)
            print(f"DEBUG: ESPRESSO_BUILD set to: {env['ESPRESSO_BUILD']}")
            print(f"DEBUG: Expected pw.x location: {expected_pw}")
            print(f"DEBUG: pw.x exists: {expected_pw.exists()}")
        else:
            env["ESPRESSO_BUILD"] = str(pw_bin_dir.parent)
        
        print(f"Using pw.x from: {pw_exe_path}")
        print(f"ESPRESSO_BUILD set to: {env['ESPRESSO_BUILD']}")
    
    # On Windows, we need to ensure Unix tools are available or use MSYS2
    if platform.system() == "Windows":
        # Add MSYS2/usr/bin to PATH for Unix tools (which, grep, etc.)
        msys2_paths = [
            "C:\\tools\\msys64\\usr\\bin",
            "C:\\msys64\\usr\\bin",
        ]
        bash_found = False
        for msys_path in msys2_paths:
            if Path(msys_path).exists():
                current_path = env.get("PATH", "")
                env["PATH"] = f"{msys_path};{current_path}"
                print(f"DEBUG: Added {msys_path} to PATH for Unix tools")
                # Check if bash is available
                bash_exe = Path(msys_path) / "bash.exe"
                if bash_exe.exists():
                    bash_found = True
                    print(f"DEBUG: Found bash at {bash_exe}")
                break
        
        if not bash_found:
            print("WARNING: bash not found in MSYS2 paths. Tests may fail.")
            print("  Please ensure MSYS2 is installed or bash is in PATH")
    
    # Build command
    cmd = [sys.executable, str(testcode_py), "--verbose"]
    
    if verbose > 1:
        cmd.extend(["-" + "v" * (verbose - 1)])
    
    cmd.extend(["--category", category])
    
    print(f"DEBUG: Running tests with command: {' '.join(cmd)}")
    print(f"DEBUG: Working directory: {test_suite_dir}")
    print(f"DEBUG: ESPRESSO_BUILD: {env.get('ESPRESSO_BUILD', 'NOT SET')}")
    print(f"DEBUG: ESPRESSO_ROOT: {env.get('ESPRESSO_ROOT', 'NOT SET')}")
    if nprocs:
        print(f"DEBUG: Using {nprocs} MPI processes")
    if nbgrp:
        print(f"DEBUG: Using {nbgrp} band groups")
    print()
    
    # Run testcode
    try:
        print("DEBUG: Starting subprocess.run...")
        result = subprocess.run(cmd, env=env, cwd=test_suite_dir, 
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              text=True, bufsize=1)
        
        # Print output in real-time
        if result.stdout:
            print(result.stdout)
        
        print(f"DEBUG: Process exited with code: {result.returncode}")
        return result.returncode == 0
    except KeyboardInterrupt:
        print("\nInterrupted by user")
        return False
    except Exception as e:
        print(f"ERROR: Failed to run tests: {e}")
        import traceback
        traceback.print_exc()
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

  # Use custom pw.x executable
  python run_qe_tests.py --pw-exe /path/to/pw.x --category pw_all
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
    
    parser.add_argument(
        "--pw-exe",
        type=Path,
        help="Path to pw.x executable (overrides ESPRESSO_BUILD/bin/pw.x)"
    )
    
    parser.add_argument(
        "--test-suite-version",
        default="qe-7.5",
        help="QE version for test-suite download (default: qe-7.5)"
    )
    
    args = parser.parse_args()
    
    # Determine test-suite location
    script_dir = Path(__file__).resolve().parent
    tests_dir = script_dir.parent
    test_suite_dir = tests_dir / "test-suite"
    
    # Download test-suite if needed
    if not test_suite_dir.exists() or not (test_suite_dir / "testcode").exists():
        print(f"Test-suite not found at {test_suite_dir}")
        if not download_test_suite(test_suite_dir, args.test_suite_version):
            print("ERROR: Failed to download test-suite")
            return 1
    
    # Find QE root (use test-suite parent if available, otherwise auto-detect)
    if args.qe_root:
        qe_root = Path(args.qe_root).resolve()
    else:
        # If we downloaded test-suite, we might not have a full QE source tree
        # In that case, use the test-suite directory as a reference
        qe_root = find_qe_root()
        if qe_root is None:
            # Use test-suite directory as QE root (tests can still run)
            qe_root = test_suite_dir.parent
            print(f"WARNING: Could not find full QE source tree, using {qe_root} as QE root")
        else:
            print(f"Found QE root: {qe_root}")
    
    # Setup environment
    env, tmp_dir = setup_environment(
        qe_root, args.qe_build, args.pseudo_dir, args.tmp_dir
    )
    
    print(f"QE Root: {qe_root}")
    print(f"QE Build: {env['ESPRESSO_BUILD']}")
    print(f"Pseudo Dir: {env['ESPRESSO_PSEUDO']}")
    print(f"Temp Dir: {tmp_dir}")
    print(f"Test Suite: {test_suite_dir}")
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
        
        if not run_tests(test_suite_dir, category, env, args.nprocs, args.nbgrp, args.verbose, args.pw_exe):
            success = False
            print(f"\nERROR: Tests failed for category {category}")
        else:
            print(f"\nSUCCESS: Tests completed for category {category}")
        print()
    
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())

