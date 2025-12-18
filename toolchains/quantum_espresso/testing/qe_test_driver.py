"""
Quantum ESPRESSO test driver.

Version-aware test harness for running QE test-suite tests.
"""

import os
import platform
import subprocess
import sys
import yaml
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Dict, List
import fnmatch


@dataclass
class QETestCase:
    """Represents a single QE test case."""
    name: str
    group: str           # e.g. "pw", "ph", etc.
    workdir: Path        # directory where the test will run
    input_file: Path     # primary input (e.g. *.in)
    reference_dir: Path  # where reference outputs live


@dataclass
class QETestResult:
    """Result of running a test case."""
    case: QETestCase
    status: str          # "PASS" | "FAIL" | "SKIPPED" | "ERROR"
    details: str = ""    # short reason or diff summary
    return_code: int = 0
    output_file: Optional[Path] = None
    error_file: Optional[Path] = None


class PlatformInfo:
    """Platform detection and utilities."""
    
    @staticmethod
    def detect() -> str:
        """Detect current platform."""
        system = platform.system()
        if system == 'Windows':
            return 'windows'
        elif system == 'Darwin':
            return 'macos'
        elif system == 'Linux':
            return 'linux'
        else:
            return 'linux'  # Default to linux for other Unix-like systems
    
    @staticmethod
    def is_windows() -> bool:
        """Check if running on Windows."""
        return platform.system() == 'Windows'
    
    @staticmethod
    def get_executable_extension() -> str:
        """Get executable extension for current platform."""
        return '.exe' if PlatformInfo.is_windows() else ''


def parse_environment_file(env_file: Path) -> Dict[str, str]:
    """
    Parse QE ENVIRONMENT file and return environment variables as a dictionary.
    
    This function parses shell script syntax like:
        export VAR=value
        export VAR=${OTHER_VAR}/path
    
    On Windows, this prevents the need to source .sh files which causes dialogs.
    
    Args:
        env_file: Path to ENVIRONMENT file
        
    Returns:
        Dictionary of environment variable name -> value
    """
    env_vars = {}
    
    if not env_file.exists():
        return env_vars
    
    try:
        with open(env_file, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # Parse export statements
        # Match: export VAR=value or export VAR=${OTHER}/path
        export_pattern = r'export\s+(\w+)=(.*?)(?=\n|$)'
        
        for match in re.finditer(export_pattern, content, re.MULTILINE):
            var_name = match.group(1)
            var_value = match.group(2).strip()
            
            # Remove quotes if present
            if var_value.startswith('"') and var_value.endswith('"'):
                var_value = var_value[1:-1]
            elif var_value.startswith("'") and var_value.endswith("'"):
                var_value = var_value[1:-1]
            
            # Handle special cases first (like ${shell pwd} which is bash command substitution)
            # Replace ${shell pwd} with the directory containing the ENVIRONMENT file
            var_value = re.sub(r'\$\{shell\s+pwd\}', str(env_file.parent), var_value)
            
            # Handle variable substitution ${VAR} or $VAR
            def substitute_var(m):
                var = m.group(1) or m.group(2)
                # Try to resolve from already parsed vars or environment
                if var in env_vars:
                    return env_vars[var]
                elif var in os.environ:
                    return os.environ[var]
                else:
                    # For ${PWD}, try to get current directory
                    if var == 'PWD':
                        return str(env_file.parent)
                    return m.group(0)  # Return original if can't resolve
            
            # Replace ${VAR} and $VAR patterns (but not ${shell pwd} which we already handled)
            var_value = re.sub(r'\$\{(\w+)\}|\$(\w+)', substitute_var, var_value)
            
            env_vars[var_name] = var_value
            
    except Exception as e:
        # If parsing fails, return empty dict (will use defaults)
        if os.environ.get('QE_TEST_VERBOSE'):
            print(f"Warning: Could not parse ENVIRONMENT file {env_file}: {e}", file=sys.stderr)
    
    return env_vars


class ExecutableManager:
    """Manages executable paths and DLL dependencies."""
    
    def __init__(self, bin_dir: Optional[Path] = None):
        """Initialize executable manager."""
        self.bin_dir = bin_dir
        if self.bin_dir is None:
            build_dir = os.environ.get('ESPRESSO_BUILD', '')
            if build_dir:
                self.bin_dir = Path(build_dir) / 'bin'
        
        if self.bin_dir:
            self.bin_dir = self.bin_dir.resolve()
    
    def get_executable_path(self, exe_name: str) -> Path:
        """Get full path to executable.
        
        Handles various naming conventions:
        - pw.x -> pw.x.exe (Windows) or pw.x (Unix)
        - pw -> pw.x.exe (Windows) or pw.x (Unix)
        - pw.exe -> pw.exe
        - pw.x.exe -> pw.x.exe
        """
        if not self.bin_dir:
            return Path(exe_name)
        
        # Normalize: remove .x and .exe to get base name
        base = exe_name.replace('.x', '').replace('.exe', '')
        
        # Try different variations
        candidates = []
        if PlatformInfo.is_windows():
            # Windows: try pw.x.exe, pw.exe, pw.x
            candidates = [
                f"{base}.x.exe",
                f"{base}.exe",
                f"{base}.x",
            ]
        else:
            # Unix: try pw.x, pw
            candidates = [
                f"{base}.x",
                base,
            ]
        
        # Also try the original name if it's different
        if exe_name not in candidates:
            candidates.insert(0, exe_name)
        
        for candidate in candidates:
            exe_path = self.bin_dir / candidate
            if exe_path.exists():
                return exe_path
        
        # Fallback: return first candidate (will fail at runtime if not found)
        return self.bin_dir / candidates[0]
    
    def get_environment_with_path(self) -> Dict[str, str]:
        """Get environment with bin_dir added to PATH."""
        env = os.environ.copy()
        
        if self.bin_dir:
            path_sep = ';' if PlatformInfo.is_windows() else ':'
            current_path = env.get('PATH', '')
            if str(self.bin_dir) not in current_path:
                env['PATH'] = f"{self.bin_dir}{path_sep}{current_path}"
        
        return env


class QETestSuite:
    """Main test suite driver."""
    
    def __init__(
        self,
        qe_version: str,
        build_dir: Path,
        test_root: Optional[Path] = None,
        platform: Optional[str] = None,
    ):
        """
        Initialize test suite.
        
        Args:
            qe_version: QE version string (e.g., "qe-7.5")
            build_dir: Directory containing QE binaries
            test_root: Optional path to test-suite root
            platform: Optional platform override ("linux", "macos", "windows")
        """
        self.qe_version = qe_version
        self.build_dir = Path(build_dir).resolve()
        self.test_root = test_root
        self.platform = platform or PlatformInfo.detect()
        
        # Load config for this version
        self.config = self._load_config()
        
        # Executable manager - determine bin directory
        bin_dir = self.build_dir
        # Check if build_dir is already a bin directory (contains executables)
        # Try various naming conventions
        has_executables = False
        for exe_variant in ['pw.x.exe', 'pw.x', 'pw.exe', 'pw']:
            if (self.build_dir / exe_variant).exists():
                has_executables = True
                break
        
        if not has_executables:
            # Try bin subdirectory
            bin_subdir = self.build_dir / 'bin'
            if bin_subdir.exists():
                bin_dir = bin_subdir
        self.exe_mgr = ExecutableManager(bin_dir)
    
    def _load_config(self) -> Dict:
        """Load per-version config file."""
        config_file = Path(__file__).parent / 'configs' / f'{self.qe_version}.yaml'
        if not config_file.exists():
            raise FileNotFoundError(
                f"Config file not found: {config_file}\n"
                f"Please create a config file for version {self.qe_version}"
            )
        
        with open(config_file, 'r') as f:
            return yaml.safe_load(f)
    
    def _find_toolchain_build_dir(self) -> Optional[Path]:
        """Try to find build directory in toolchain structure."""
        # Find repo root (toolchains/quantum_espresso)
        current = Path(__file__).resolve()
        for _ in range(10):
            if current.name == 'quantum_espresso' and current.parent.name == 'toolchains':
                qe_root = current
                break
            if current.parent == current:
                return None
            current = current.parent
        else:
            return None
        
        # Try platform-specific locations
        platform_dirs = {
            'linux': qe_root / 'linux',
            'macos': qe_root / 'macos',
            'windows': qe_root / 'windows',
        }
        
        platform_dir = platform_dirs.get(self.platform)
        if not platform_dir or not platform_dir.exists():
            return None
        
        # Try common build locations
        candidates = []
        
        # For Windows, try both mingw and oneapi
        if self.platform == 'windows':
            for variant in ['mingw', 'oneapi']:
                variant_dir = platform_dir / variant
                if variant_dir.exists():
                    # Try dist directory first (installed binaries)
                    dist_dir = variant_dir / 'dist' / f'win-{variant}'
                    if dist_dir.exists():
                        candidates.append(dist_dir)
                    # Try upstream source bin
                    for qe_dir_name in [f'q-e-{self.qe_version}', 'qe']:
                        upstream_bin = variant_dir / 'upstream' / qe_dir_name / 'bin'
                        if upstream_bin.exists():
                            candidates.append(upstream_bin)
        else:
            # Linux/macOS
            for qe_dir_name in [f'q-e-{self.qe_version}', 'qe']:
                upstream_bin = platform_dir / 'upstream' / qe_dir_name / 'bin'
                if upstream_bin.exists():
                    candidates.append(upstream_bin)
        
        # Return first existing candidate
        for candidate in candidates:
            if candidate.exists():
                return candidate
        
        return None
    
    def ensure_test_suite_present(self) -> Path:
        """
        Ensure QE test-suite is available for this version.
        
        Returns:
            Path to test-suite root
        """
        if self.test_root and self.test_root.exists():
            test_suite_dir = self.test_root
            if (test_suite_dir / 'testcode').exists():
                return test_suite_dir
        
        # Try to find in toolchain structure
        repo_root = Path(__file__).resolve()
        for _ in range(10):
            if repo_root.name == 'quantum_espresso' and repo_root.parent.name == 'toolchains':
                break
            if repo_root.parent == repo_root:
                break
            repo_root = repo_root.parent
        
        qe_root = repo_root.parent / 'quantum_espresso'
        
        # Try existing test-suite location
        existing_test_suite = qe_root / 'tests' / 'test-suite'
        if existing_test_suite.exists() and (existing_test_suite / 'testcode').exists():
            return existing_test_suite
        
        # Try in upstream source trees
        platform_dirs = {
            'linux': qe_root / 'linux',
            'macos': qe_root / 'macos',
            'windows': qe_root / 'windows',
        }
        
        for platform_name, platform_dir in platform_dirs.items():
            if not platform_dir.exists():
                continue
            
            # Try mingw/oneapi for Windows
            if platform_name == 'windows':
                for variant in ['mingw', 'oneapi']:
                    variant_dir = platform_dir / variant / 'upstream'
                    for qe_dir_name in [f'q-e-{self.qe_version}', 'qe']:
                        qe_dir = variant_dir / qe_dir_name
                        test_suite = qe_dir / 'test-suite'
                        if test_suite.exists() and (test_suite / 'testcode').exists():
                            return test_suite
            else:
                # Linux/macOS
                upstream_dir = platform_dir / 'upstream'
                for qe_dir_name in [f'q-e-{self.qe_version}', 'qe']:
                    qe_dir = upstream_dir / qe_dir_name
                    test_suite = qe_dir / 'test-suite'
                    if test_suite.exists() and (test_suite / 'testcode').exists():
                        return test_suite
        
        # If not found, provide helpful error
        raise FileNotFoundError(
            f"QE test-suite not found for version {self.qe_version}.\n"
            f"Please either:\n"
            f"  1. Set test_root parameter to point to test-suite directory\n"
            f"  2. Ensure QE source tree exists at:\n"
            f"     - toolchains/quantum_espresso/tests/test-suite/\n"
            f"     - toolchains/quantum_espresso/<platform>/upstream/q-e-{self.qe_version}/test-suite\n"
            f"  3. Download test-suite manually from QE repository"
        )
    
    def discover_tests(self, groups: Optional[Iterable[str]] = None) -> List[QETestCase]:
        """
        Discover available tests using the per-version config.
        
        Args:
            groups: Optional list of test groups to include (e.g., ["pw", "ph"])
        
        Returns:
            List of QETestCase objects
        """
        test_suite_root = self.ensure_test_suite_present()
        tests = []
        
        config_groups = self.config.get('groups', {})
        if groups:
            config_groups = {k: v for k, v in config_groups.items() if k in groups}
        
        for group_name, group_config in config_groups.items():
            group_root = test_suite_root / group_config.get('root', f'test-suite/{group_name}_*')
            
            # Handle glob patterns in root
            if '*' in str(group_root):
                # Find all matching directories
                parent = group_root.parent
                pattern = group_root.name
                if parent.exists():
                    for item in parent.iterdir():
                        if item.is_dir() and fnmatch.fnmatch(item.name, pattern):
                            self._discover_tests_in_dir(item, group_name, tests, group_config)
            else:
                if group_root.exists():
                    self._discover_tests_in_dir(group_root, group_name, tests, group_config)
        
        return tests
    
    def _discover_tests_in_dir(
        self,
        test_dir: Path,
        group_name: str,
        tests: List[QETestCase],
        group_config: Dict
    ):
        """Discover tests in a specific directory."""
        patterns = group_config.get('patterns', ['*.in'])
        
        for pattern in patterns:
            # Handle directory patterns like "pw_*"
            if '*' in pattern and not pattern.endswith('.in'):
                # Directory pattern
                parent = test_dir
                if not parent.exists():
                    continue
                for item in parent.iterdir():
                    if item.is_dir() and fnmatch.fnmatch(item.name, pattern):
                        # Look for .in files in this directory
                        for in_file in item.glob('*.in'):
                            if in_file.name.startswith('test.') or in_file.name.startswith('benchmark.'):
                                continue
                            tests.append(QETestCase(
                                name=f"{group_name}_{item.name}",
                                group=group_name,
                                workdir=item,
                                input_file=in_file,
                                reference_dir=item
                            ))
            else:
                # File pattern
                for in_file in test_dir.rglob(pattern):
                    if in_file.name.startswith('test.') or in_file.name.startswith('benchmark.'):
                        continue
                    # Use parent directory as workdir
                    tests.append(QETestCase(
                        name=f"{group_name}_{in_file.parent.name}_{in_file.stem}",
                        group=group_name,
                        workdir=in_file.parent,
                        input_file=in_file,
                        reference_dir=in_file.parent
                    ))
    
    def _check_pseudopotentials(self, tests: List[QETestCase]) -> bool:
        """
        Check and download missing pseudopotentials for test cases.
        
        Args:
            tests: List of test cases to check
            
        Returns:
            True if all pseudopotentials are available, False otherwise
        """
        test_suite_root = self.ensure_test_suite_present()
        pseudo_checker_path = test_suite_root / 'qe_pseudo.py'
        
        if not pseudo_checker_path.exists():
            # If checker not available, skip check
            return True
        
        # Import the pseudopotential checker
        try:
            import importlib.util
            spec = importlib.util.spec_from_file_location("qe_pseudo", pseudo_checker_path)
            qe_pseudo = importlib.util.module_from_spec(spec)
            sys.modules["qe_pseudo"] = qe_pseudo
            spec.loader.exec_module(qe_pseudo)
            
            # Set up pseudo directory
            pseudo_dir = test_suite_root.parent / 'pseudo'
            pseudo_dir.mkdir(exist_ok=True)
            
            # Extract all pseudopotential files from test inputs
            pp_files = set()
            for test in tests:
                try:
                    with open(test.input_file, 'r') as f:
                        lines = f.readlines()
                        in_atomic_species = False
                        for line in lines:
                            line_upper = line.upper().strip()
                            if 'ATOMIC_SPECIES' in line_upper:
                                in_atomic_species = True
                                continue
                            # Stop at next section
                            if in_atomic_species and (line.strip().startswith('&') or 
                                                      'ATOMIC_POSITIONS' in line_upper or
                                                      'K_POINTS' in line_upper or
                                                      'CELL_PARAMETERS' in line_upper):
                                in_atomic_species = False
                                continue
                            # Extract pseudopotential filename from ATOMIC_SPECIES section
                            if in_atomic_species and line.strip():
                                parts = line.split()
                                if len(parts) >= 3:
                                    # Last part should be the filename
                                    filename = parts[-1].strip()
                                    if filename.endswith('.UPF') or filename.endswith('.upf'):
                                        pp_files.add(filename)
                except Exception as e:
                    pass  # Skip if can't read file
            
            # Check and download missing pseudopotentials
            checker = qe_pseudo.PseudopotentialChecker(
                pseudo_dir=str(pseudo_dir),
                network_url='https://pseudopotentials.quantum-espresso.org/upf_files/'
            )
            
            all_ok = True
            for pp_file in pp_files:
                local_path = pseudo_dir / pp_file
                if not local_path.exists():
                    if not checker.download_pseudopotential(pp_file):
                        all_ok = False
                        break
            
            return all_ok
        except Exception as e:
            print(f"Warning: Could not check pseudopotentials: {e}", file=sys.stderr)
            import traceback
            if os.environ.get('QE_TEST_VERBOSE'):
                traceback.print_exc()
            return True  # Continue anyway
    
    def run_tests(
        self,
        pattern: Optional[str] = None,
        groups: Optional[Iterable[str]] = None,
        jobs: int = 1,
        env: Optional[Dict] = None,
        limit: Optional[int] = None,
    ) -> List[QETestResult]:
        """
        Run tests serially (parallel support TODO).
        
        Args:
            pattern: Optional glob/substring pattern to match test names
            groups: Optional list of test groups to run
            jobs: Number of parallel jobs (currently ignored, runs serially)
            env: Optional environment variables to set
            limit: Optional limit on number of tests to run
        
        Returns:
            List of QETestResult objects
        """
        tests = self.discover_tests(groups)
        
        # Filter by pattern if provided
        if pattern:
            tests = [t for t in tests if pattern in t.name or fnmatch.fnmatch(t.name, pattern)]
        
        # Apply limit if provided
        if limit is not None and limit > 0:
            tests = tests[:limit]
        
        # Check and download pseudopotentials before running tests
        self._check_pseudopotentials(tests)
        
        results = []
        for test in tests:
            result = self._run_single_test(test, env)
            results.append(result)
        
        return results
    
    def _run_single_test(self, test: QETestCase, env: Optional[Dict] = None) -> QETestResult:
        """Run a single test case."""
        # Get executable for this test group
        executable_map = self.config.get('defaults', {}).get('executable_map', {})
        exe_name = executable_map.get(test.group, f"{test.group}.x")
        
        exe_path = self.exe_mgr.get_executable_path(exe_name)
        if not exe_path.exists():
            return QETestResult(
                case=test,
                status="ERROR",
                details=f"Executable not found: {exe_path}",
                return_code=-1
            )
        
        # Setup output files
        output_file = test.workdir / f"{test.input_file.stem}.out"
        error_file = test.workdir / f"{test.input_file.stem}.err"
        
        # Build command
        # On Windows with Intel Fortran, we need to handle stack overflow issues
        # Use a wrapper script or set stack size via editbin
        if PlatformInfo.is_windows() and str(exe_path).endswith('.exe'):
            # Check if we need to use editbin to increase stack size
            # For now, try running directly - if it fails with stack overflow,
            # we'll need to use editbin or a wrapper
            cmd = [str(exe_path), '-input', str(test.input_file)]
        else:
            cmd = [str(exe_path), '-input', str(test.input_file)]
        
        # Get MPI config from environment
        mpi_prefix = None
        nprocs = os.environ.get('QE_USE_MPI', '')
        if nprocs:
            mpi_prefix = ['mpirun', '-np', nprocs]
            cmd = mpi_prefix + cmd
        
        # Setup environment
        run_env = self.exe_mgr.get_environment_with_path()
        
        # Set up pseudopotential directory
        test_suite_root = self.ensure_test_suite_present()
        pseudo_dir = test_suite_root.parent / 'pseudo'
        pseudo_dir.mkdir(exist_ok=True)
        run_env['ESPRESSO_PSEUDO'] = str(pseudo_dir)
        run_env.setdefault('NETWORK_PSEUDO', 'https://pseudopotentials.quantum-espresso.org/upf_files/')
        
        # Set ESPRESSO_ROOT and ESPRESSO_BUILD if not already set
        if 'ESPRESSO_ROOT' not in run_env:
            run_env['ESPRESSO_ROOT'] = str(test_suite_root.parent)
        if 'ESPRESSO_BUILD' not in run_env:
            run_env['ESPRESSO_BUILD'] = str(self.build_dir)
        
        # Load environment variables from ENVIRONMENT file if it exists
        # This prevents Windows from trying to open .sh files
        env_file = test_suite_root / 'ENVIRONMENT'
        if env_file.exists():
            env_from_file = parse_environment_file(env_file)
            # Merge with existing environment (file takes precedence for specific vars)
            for key, value in env_from_file.items():
                if key not in run_env:  # Don't override already set vars
                    run_env[key] = value
        
        if env:
            run_env.update(env)
        
        # On Windows, prevent .sh files from being executed
        # Ensure shell=False to prevent Windows from trying to open .sh files
        # Also ensure we never pass .sh files to subprocess
        shell = False
        
        # On Windows, ensure executable path is not a .sh file
        if PlatformInfo.is_windows() and str(exe_path).endswith('.sh'):
            return QETestResult(
                case=test,
                status="ERROR",
                details=f"Cannot execute .sh file on Windows: {exe_path}. Use Python wrapper instead.",
                return_code=-1
            )
        
        # On Windows, check all command arguments for .sh files and prevent execution
        if PlatformInfo.is_windows():
            for arg in cmd:
                if str(arg).endswith('.sh'):
                    return QETestResult(
                        case=test,
                        status="ERROR",
                        details=f"Cannot execute .sh file on Windows: {arg}. Use Python wrapper instead.",
                        return_code=-1
                    )
        
        # Run test
        try:
            # On Windows, use CREATE_NO_WINDOW flag to prevent dialogs
            creation_flags = 0
            if PlatformInfo.is_windows():
                # Prevent Windows from showing dialogs when file associations fail
                creation_flags = subprocess.CREATE_NO_WINDOW
            
            with open(output_file, 'w') as out_f, open(error_file, 'w') as err_f:
                result = subprocess.run(
                    cmd,
                    cwd=str(test.workdir),
                    env=run_env,
                    stdout=out_f,
                    stderr=err_f,
                    timeout=300,  # 5 minute timeout per test
                    shell=shell,
                    creationflags=creation_flags
                )
            
            # Check for CRASH file
            crash_file = test.workdir / 'CRASH'
            if crash_file.exists():
                return QETestResult(
                    case=test,
                    status="FAIL",
                    details="CRASH file detected",
                    return_code=result.returncode,
                    output_file=output_file,
                    error_file=error_file
                )
            
            # Simple pass/fail based on return code
            if result.returncode == 0:
                return QETestResult(
                    case=test,
                    status="PASS",
                    details="Test completed successfully",
                    return_code=0,
                    output_file=output_file,
                    error_file=error_file
                )
            else:
                return QETestResult(
                    case=test,
                    status="FAIL",
                    details=f"Executable returned code {result.returncode}",
                    return_code=result.returncode,
                    output_file=output_file,
                    error_file=error_file
                )
        
        except subprocess.TimeoutExpired:
            return QETestResult(
                case=test,
                status="ERROR",
                details="Test timed out after 5 minutes",
                return_code=-1,
                output_file=output_file,
                error_file=error_file
            )
        except Exception as e:
            return QETestResult(
                case=test,
                status="ERROR",
                details=f"Exception: {str(e)}",
                return_code=-1,
                output_file=output_file,
                error_file=error_file
            )

