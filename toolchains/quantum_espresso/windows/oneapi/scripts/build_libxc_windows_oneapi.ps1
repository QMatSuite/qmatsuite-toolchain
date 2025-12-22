# scripts/build_libxc_windows_oneapi.ps1
# Windows PowerShell 5.1 compatible
# Build and install libxc (shared/DLL) with Intel oneAPI compilers (icx/ifx),
# using the latest stable GitLab Release (NOT master).

param(
  # GitLab official repo
  [string]$RepoUrl = "https://gitlab.com/libxc/libxc.git",

  # Use latest release tag by default (resolved via GitLab API).
  # You may override with -Ref <tag/commit>.
  [string]$Ref = "",

  [switch]$Clean,
  [switch]$NoUpdate,
  [switch]$EnableTests,  # default false; tests require bzip2 in some libxc versions
  [ValidateSet("Release","RelWithDebInfo","Debug")]
  [string]$Config = "Release"
)

$ErrorActionPreference = "Stop"

function Write-Section([string]$msg) {
  Write-Host ""
  Write-Host ("=== " + $msg + " ===") -ForegroundColor Cyan
}

function Require-Command([string]$name, [string]$hint) {
  $cmd = Get-Command $name -ErrorAction SilentlyContinue
  if (-not $cmd) { throw ("Missing required command: " + $name + ". " + $hint) }
  return $cmd.Source
}

function Import-OneApiEnv([string]$OneApiRoot) {
  $setvars = Join-Path $OneApiRoot "setvars.bat"
  if (-not (Test-Path $setvars)) {
    throw ("oneAPI setvars.bat not found: " + $setvars)
  }

  Write-Host ("Initializing oneAPI env via: " + $setvars) -ForegroundColor Yellow

  # Use cmd.exe to avoid PS native stderr->exception quirks
  $cmdline = "call `"$setvars`" intel64 >nul 2>nul && set"
  $lines = & cmd /c $cmdline

  foreach ($line in $lines) {
    if ($line -match "^(.*?)=(.*)$") {
      [System.Environment]::SetEnvironmentVariable($matches[1], $matches[2], "Process")
    }
  }
}

function Get-LatestLibxcReleaseTag {
  # GitLab Releases API: /projects/:id/releases/permalink/latest
  # Public projects do not require auth.
  $projectPath = "libxc/libxc"
  $encoded = [System.Uri]::EscapeDataString($projectPath)
  $api = "https://gitlab.com/api/v4/projects/$encoded/releases/permalink/latest"

  Write-Host ("Querying GitLab latest release: " + $api) -ForegroundColor DarkGray
  try {
    $resp = Invoke-RestMethod -Method Get -Uri $api -Headers @{ "User-Agent"="qmatsuite-toolchain-libxc" } -TimeoutSec 30
    if ($resp -and $resp.tag_name) { return [string]$resp.tag_name }
  } catch {
    Write-Host ("WARNING: Failed to query latest release from GitLab API: " + $_.Exception.Message) -ForegroundColor Yellow
  }

  # Fallback: known latest on libxc download/release pages (update if needed)
  return "7.0.0"
}

# ----------------------------------------------------------------------
# Resolve repo layout: .../oneapi/scripts/<this file>
# ----------------------------------------------------------------------
$ScriptPath = $MyInvocation.MyCommand.Path
$ScriptsDir = Split-Path -Parent $ScriptPath
$OneApiDir  = Split-Path -Parent $ScriptsDir

$UpstreamDir = Join-Path $OneApiDir "upstream"
$LibxcRoot   = Join-Path $UpstreamDir "libxc"
$SrcDir      = Join-Path $LibxcRoot "src"
$BuildDir    = Join-Path $LibxcRoot "build-win-oneapi"
$InstallDir  = Join-Path $LibxcRoot "install"

Write-Section "libxc Windows oneAPI build (GitLab latest release)"
Write-Host ("OneAPI dir   : " + $OneApiDir)
Write-Host ("RepoUrl      : " + $RepoUrl)
Write-Host ("Config       : " + $Config)
Write-Host ("Clean        : " + $Clean)
Write-Host ("NoUpdate     : " + $NoUpdate)
Write-Host ("EnableTests  : " + $EnableTests)

New-Item -ItemType Directory -Force -Path $UpstreamDir | Out-Null
New-Item -ItemType Directory -Force -Path $LibxcRoot   | Out-Null

# Resolve latest release tag if Ref not provided
if ([string]::IsNullOrWhiteSpace($Ref)) {
  $Ref = Get-LatestLibxcReleaseTag
  Write-Host ("Resolved latest release tag: " + $Ref) -ForegroundColor Green
} else {
  Write-Host ("Using user-provided Ref: " + $Ref) -ForegroundColor Green
}

# ----------------------------------------------------------------------
# oneAPI environment
# ----------------------------------------------------------------------
Write-Section "Discover Intel oneAPI"
$oneApiRoot = $env:ONEAPI_ROOT
if (-not $oneApiRoot) { $oneApiRoot = "C:\Program Files (x86)\Intel\oneAPI" }
if (-not (Test-Path $oneApiRoot)) {
  throw ("Intel oneAPI not found at: " + $oneApiRoot + " (set ONEAPI_ROOT or install oneAPI)")
}
Import-OneApiEnv -OneApiRoot $oneApiRoot

# ----------------------------------------------------------------------
# Tooling
# ----------------------------------------------------------------------
Write-Section "Check tools"
$cmakeExe = Require-Command "cmake" "Install CMake and ensure it is on PATH."
$ninjaExe = Require-Command "ninja" "Install Ninja and ensure it is on PATH."
$gitExe   = Require-Command "git"   "Install Git and ensure it is on PATH."

$icxCmd = Get-Command icx -ErrorAction SilentlyContinue
$ifxCmd = Get-Command ifx -ErrorAction SilentlyContinue
$icxPath = if ($icxCmd) { $icxCmd.Source } else { $null }
$ifxPath = if ($ifxCmd) { $ifxCmd.Source } else { $null }
if (-not $icxPath -or -not $ifxPath) {
  throw "icx/ifx not found in PATH even after setvars.bat. Ensure oneAPI HPC Toolkit is installed."
}
Write-Host ("icx: " + $icxPath) -ForegroundColor Green
Write-Host ("ifx: " + $ifxPath) -ForegroundColor Green

# ----------------------------------------------------------------------
# Fetch/update source
# ----------------------------------------------------------------------
Write-Section "Fetch libxc source"
if (-not (Test-Path $SrcDir)) {
  & $gitExe clone $RepoUrl $SrcDir
}

Push-Location $SrcDir
try {
  if (-not $NoUpdate) {
    & $gitExe fetch --tags --prune
  }

  # Checkout release tag/commit
  & $gitExe checkout $Ref

  # Only pull when on a branch (tags/commits are detached HEAD)
  $head = & $gitExe rev-parse --abbrev-ref HEAD 2>$null
  if ($head -and $head.Trim() -ne "HEAD") {
    & $gitExe pull --ff-only 2>$null | Out-Null
  } else {
    Write-Host "Detached HEAD (tag/commit checkout). Skip git pull." -ForegroundColor DarkGray
  }


  & $gitExe submodule update --init --recursive 2>$null | Out-Null
}
finally {
  Pop-Location
}

# ----------------------------------------------------------------------
# Clean (optional)
# ----------------------------------------------------------------------
if ($Clean) {
  Write-Section "Clean build/install"
  if (Test-Path $BuildDir)   { Remove-Item -Recurse -Force $BuildDir }
  if (Test-Path $InstallDir) { Remove-Item -Recurse -Force $InstallDir }
}

# ----------------------------------------------------------------------
# Configure / Build / Install
# ----------------------------------------------------------------------
Write-Section "CMake configure"
New-Item -ItemType Directory -Force -Path $BuildDir   | Out-Null
New-Item -ItemType Directory -Force -Path $InstallDir | Out-Null

# CMake 4 removed compatibility with <3.5, older projects may need this escape hatch.
# Harmless for modern projects.
$policyMin = "-DCMAKE_POLICY_VERSION_MINIMUM=3.5"

# Tests off by default (avoids bzip2 dependency in some libxc versions)
$testArgs = @("-DBUILD_TESTING=OFF","-DENABLE_TESTS=OFF")
if ($EnableTests) {
  $bzip = Get-Command bzip2 -ErrorAction SilentlyContinue
  if (-not $bzip) { throw "EnableTests requested but 'bzip2' not found on PATH. Install bzip2 or rerun without -EnableTests." }
  $testArgs = @() # let libxc default decide / or user toggles inside source
}

# Force outputs into build root to make artifacts easy to locate
$outArgs = @(
  "-DCMAKE_ARCHIVE_OUTPUT_DIRECTORY=$BuildDir",  # .lib import/archives
  "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=$BuildDir",  # .dll (library outputs)
  "-DCMAKE_RUNTIME_OUTPUT_DIRECTORY=$BuildDir"   # .exe/.dll runtime outputs
)

$cmakeArgs = @(
  "-S", $SrcDir,
  "-B", $BuildDir,
  "-G", "Ninja",
  "-DCMAKE_POLICY_VERSION_MINIMUM=3.5",
  "-DCMAKE_BUILD_TYPE=$Config",
  "-DCMAKE_INSTALL_PREFIX=$InstallDir",
  "-DCMAKE_C_COMPILER=$icxPath",
  "-DCMAKE_Fortran_COMPILER=$ifxPath",
  "-DBUILD_SHARED_LIBS=ON",
  "-DENABLE_FORTRAN=ON",
  "-DBUILD_TESTING=OFF",
  "-DCMAKE_WINDOWS_EXPORT_ALL_SYMBOLS=ON"
  $policyMin
) + $testArgs + $outArgs

Write-Host ("CMake args: " + ($cmakeArgs -join " ")) -ForegroundColor DarkGray
& $cmakeExe @cmakeArgs

Write-Section "Build"
& $cmakeExe --build $BuildDir --config $Config

# Sanity check: xc.lib must exist for normal link usage on Windows
Write-Section "Sanity check (xc.lib)"
$xcLib = Join-Path $BuildDir "xc.lib"
if (-not (Test-Path $xcLib)) {
  Write-Host "ERROR: xc.lib not found in build root." -ForegroundColor Red
  Write-Host "Searching for xc.lib recursively..." -ForegroundColor Yellow
  $hits = Get-ChildItem $BuildDir -Recurse -Filter "xc.lib" -ErrorAction SilentlyContinue | Select-Object -First 10
  if ($hits -and $hits.Count -gt 0) {
    foreach ($h in $hits) { Write-Host ("  found: " + $h.FullName) }
  } else {
    Write-Host "  (no xc.lib found anywhere under build dir)" -ForegroundColor Yellow
  }

  $bn = Join-Path $BuildDir "build.ninja"
  if (Test-Path $bn) {
    Write-Host "Dumping first build.ninja lines containing 'xc.dll' (for /IMPLIB diagnostics):" -ForegroundColor Yellow
    Select-String -Path $bn -Pattern "xc.dll" | Select-Object -First 30 | ForEach-Object { Write-Host $_.Line }
  }
  throw "xc.lib missing: import library required for downstream link on Windows."
}

Write-Section "Install"
& $cmakeExe --install $BuildDir --config $Config

# ----------------------------------------------------------------------
# Report + hints
# ----------------------------------------------------------------------
Write-Section "Result"
$binDir = Join-Path $InstallDir "bin"
$libDir = Join-Path $InstallDir "lib"

Write-Host ("Install prefix : " + $InstallDir)
Write-Host ("DLL dir        : " + $binDir)
Write-Host ("lib dir        : " + $libDir)

$libxcDlls = @()
if (Test-Path $binDir) {
  $libxcDlls = Get-ChildItem $binDir -Filter "*.dll" -ErrorAction SilentlyContinue | Where-Object { $_.Name -match "xc" } | Select-Object -First 20
}
Write-Host ("DLL sample     : " + $libxcDlls.Count)
foreach ($d in $libxcDlls) { Write-Host ("  - " + $d.Name) -ForegroundColor Green }

# Try to locate LibxcConfig.cmake
$libxcDir = $null
$candidateDirs = @(
  (Join-Path $InstallDir "lib\cmake\Libxc"),
  (Join-Path $InstallDir "lib64\cmake\Libxc"),
  (Join-Path $InstallDir "share\cmake\Libxc"),
  (Join-Path $InstallDir "share\cmake\Libxc\Libxc")
)
foreach ($d in $candidateDirs) {
  $cfg = Join-Path $d "LibxcConfig.cmake"
  if (Test-Path $cfg) { $libxcDir = $d; break }
}
if ($libxcDir) {
  Write-Host ("Libxc_DIR      : " + $libxcDir) -ForegroundColor Green
} else {
  Write-Host "LibxcConfig.cmake not found (not fatal). Use CMAKE_PREFIX_PATH=<install> when configuring QE." -ForegroundColor Yellow
}

# Write helper env snippet for QE
$helperLines = @()
$helperLines += "# Auto-generated helper for QE + libxc"
$helperLines += "# Usage:  . `"$($OneApiDir)\upstream\libxc\libxc-env.ps1`""
$helperLines += "`$env:CMAKE_PREFIX_PATH = `"$InstallDir;`$env:CMAKE_PREFIX_PATH`""
if ($libxcDir) { $helperLines += "`$env:Libxc_DIR = `"$libxcDir`"" }

$helperPath = Join-Path $LibxcRoot "libxc-env.ps1"
[System.IO.File]::WriteAllText($helperPath, ($helperLines -join "`r`n"), (New-Object System.Text.UTF8Encoding($false)))
Write-Host ("Wrote helper   : " + $helperPath) -ForegroundColor Green

Write-Host ""
Write-Host ("Pinned libxc release tag: " + $Ref) -ForegroundColor Cyan
Write-Host "Next: configure QE with:" -ForegroundColor Cyan
Write-Host ("  -DQE_ENABLE_LIBXC=ON -DCMAKE_PREFIX_PATH=`"" + $InstallDir + ";...`"") -ForegroundColor Cyan
if ($libxcDir) {
  Write-Host ("  (optional) -DLibxc_DIR=`"" + $libxcDir + "`"") -ForegroundColor Cyan
}
Write-Host ""
Write-Host ("For distribution: bundle libxc DLL(s) from " + $binDir + " into dist\bin, and add libxc license under dist\licenses.") -ForegroundColor Yellow
