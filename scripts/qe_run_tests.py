#!/usr/bin/env python3
"""
Thin wrapper CLI for Quantum ESPRESSO test harness.

This script is the main entry point for running QE tests from CI or local development.
"""

import sys
from pathlib import Path

# Add toolchains to path
repo_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(repo_root))

from toolchains.quantum_espresso.testing.qe_tests_cli import main

if __name__ == '__main__':
    sys.exit(main())



