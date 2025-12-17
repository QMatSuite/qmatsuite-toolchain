# Wannier90 Source Refresh Script
# Safely cleans Wannier90 upstream sources and clones fresh Wannier90 from git
# Usage: .\scripts\refresh_wannier90_source.ps1 [-DryRun] [-W90Version <branch/tag>]
# Default: clones "develop" branch (CMake support). Use -W90Version v3.1.0 for Make build.

param(
    [switch]$DryRun,           # If specified, only show what would be deleted
    [string]$W90Version = "develop"  # Specific Wannier90 version/branch to use (default: "develop" for CMake support). Use "v3.1.0" for Make build.
)

$ErrorActionPreference = 'Stop'

# Determine script location and repo root
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$RepoRoot = Resolve-Path (Join-Path $ScriptDir "..")
$W90UpstreamDir = Join-Path $RepoRoot "upstream\wannier90"

Write-Host "==========================================" -ForegroundColor Cyan
Write-Host "Wannier90 Source Refresh" -ForegroundColor Cyan
Write-Host "==========================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "Repository root: $RepoRoot"
Write-Host "Wannier90 upstream directory: $W90UpstreamDir"
if ($DryRun) {
    Write-Host "DRY RUN MODE - No files will be deleted" -ForegroundColor Yellow
}
Write-Host ""

# Define allowlist: directories and files that are "our code" and should be preserved
$AllowlistDirs = @(
    "scripts",
    ".github",
    "docs",
    "upstream"  # Preserve upstream/ directory structure if it exists
)

$AllowlistFiles = @(
    ".editorconfig",
    ".gitignore",
    ".gitattributes",
    "LICENSE",
    "LICENSE.txt",
    "README.md"
)

# Define Wannier90 source directories to remove
$W90SourceDirs = @(
    "src", "utility", "config", "cmake", "docs", "pseudo", "pwscf",
    "test-suite", "tutorials", "autodoc"
)

# Define Wannier90 config files to remove
$W90ConfigFiles = @(
    "CMakeLists.txt",
    "Makefile",
    "make.inc",
    "configure",
    "README",
    "README.install",
    "LICENSE",
    "LICENCE",
    "COPYING"
)

# Define build artifact patterns
$BuildArtifactPatterns = @(
    "build*",
    "cmake-build-*",
    "CMakeCache.txt",
    "CMakeFiles",
    "*.obj",
    "*.mod",
    "*.a",
    "*.lib",
    "*.dll",
    "*.exe",
    "*.pdb",
    "*.ilk",
    "wannier90.x",
    "postw90.x",
    "wannier90.exe",
    "postw90.exe"
)

# Step 1: Cleanup Wannier90 sources and build artifacts
Write-Host "Step 1: Cleaning up Wannier90 sources and build artifacts..." -ForegroundColor Yellow
Write-Host ""

$itemsToDelete = @()
$warnings = @()

# Collect Wannier90 source directories
foreach ($dir in $W90SourceDirs) {
    $dirPath = Join-Path $RepoRoot $dir
    if (Test-Path $dirPath) {
        $itemsToDelete += @{Type="Directory"; Path=$dirPath}
    }
}

# Collect Wannier90 config files
foreach ($file in $W90ConfigFiles) {
    $filePath = Join-Path $RepoRoot $file
    if (Test-Path $filePath) {
        $itemsToDelete += @{Type="File"; Path=$filePath}
    }
}

# Collect build artifacts in upstream/wannier90 if it exists
if (Test-Path $W90UpstreamDir) {
    foreach ($pattern in $BuildArtifactPatterns) {
        $matches = Get-ChildItem -Path $W90UpstreamDir -Filter $pattern -Recurse -ErrorAction SilentlyContinue
        foreach ($match in $matches) {
            $itemsToDelete += @{Type="File"; Path=$match.FullName}
        }
    }
    
    # Also check for build directories
    $buildDirs = Get-ChildItem -Path $W90UpstreamDir -Directory -Filter "build*" -ErrorAction SilentlyContinue
    foreach ($buildDir in $buildDirs) {
        $itemsToDelete += @{Type="Directory"; Path=$buildDir.FullName}
    }
}

# Remove duplicates
$itemsToDelete = $itemsToDelete | Sort-Object -Property Path -Unique

if ($itemsToDelete.Count -gt 0) {
    Write-Host "Found $($itemsToDelete.Count) items to delete:" -ForegroundColor Cyan
    foreach ($item in $itemsToDelete) {
        Write-Host "  [$($item.Type)] $($item.Path)"
    }
    Write-Host ""
    
    if ($DryRun) {
        Write-Host "DRY RUN: Would delete $($itemsToDelete.Count) items" -ForegroundColor Cyan
        Write-Host ""
    } else {
        Write-Host "Deleting $($itemsToDelete.Count) items..." -ForegroundColor Yellow
        foreach ($item in $itemsToDelete) {
            try {
                if ($item.Type -eq "Directory") {
                    Remove-Item -Path $item.Path -Recurse -Force -ErrorAction Stop
                } else {
                    Remove-Item -Path $item.Path -Force -ErrorAction Stop
                }
                Write-Host "  Deleted: $($item.Path)" -ForegroundColor Green
            } catch {
                Write-Host "  ERROR deleting $($item.Path): $($_.Exception.Message)" -ForegroundColor Red
                throw
            }
        }
        Write-Host "Cleanup completed." -ForegroundColor Green
        Write-Host ""
    }
} else {
    Write-Host "No Wannier90 sources or build artifacts found to clean." -ForegroundColor Green
    Write-Host ""
}

# Step 2: Determine Wannier90 version to clone
Write-Host "Step 2: Determining Wannier90 version to clone..." -ForegroundColor Yellow

$w90RepoUrl = "https://github.com/wannier-developers/wannier90.git"
$targetVersion = $null

if ($W90Version) {
    $targetVersion = $W90Version
    Write-Host "  Using specified version: $targetVersion" -ForegroundColor Green
} else {
    Write-Host "  Using default version: develop (CMake support)" -ForegroundColor Green
    $targetVersion = "develop"
}
Write-Host ""

# Step 3: Clone or update Wannier90 repository
if (-not $DryRun) {
    Write-Host "Step 3: Updating Wannier90 repository..." -ForegroundColor Yellow
    Write-Host "  Repository: $w90RepoUrl" -ForegroundColor Gray
    Write-Host "  Version: $targetVersion" -ForegroundColor Gray
    Write-Host "  Destination: $W90UpstreamDir" -ForegroundColor Gray
    Write-Host ""
    
    # Check if upstream/wannier90 is already a git repository
    $isGitRepo = Test-Path (Join-Path $W90UpstreamDir ".git")
    
    if ($isGitRepo) {
        # Reuse existing clone: fetch, checkout
        Write-Host "  Reusing existing git repository..." -ForegroundColor Green
        Write-Host ""
        
        Push-Location $W90UpstreamDir
        try {
            # Temporarily disable error action to handle git's stderr output
            $oldErrorAction = $ErrorActionPreference
            $ErrorActionPreference = 'Continue'
            
            # Fetch latest
            Write-Host "  Fetching from remote..." -ForegroundColor Gray
            $null = & git fetch origin 2>&1
            if ($LASTEXITCODE -ne 0) {
                throw "Failed to fetch from remote"
            }
            
            # Checkout target version
            Write-Host "  Checking out $targetVersion..." -ForegroundColor Gray
            $null = & git checkout $targetVersion 2>&1
            if ($LASTEXITCODE -ne 0) {
                throw "Failed to checkout $targetVersion"
            }
            
            # Reset to remote version (discard local changes)
            # Try origin/branch first, then fall back to tag
            $null = & git reset --hard origin/$targetVersion 2>&1
            if ($LASTEXITCODE -ne 0) {
                # If branch doesn't exist, try tag directly
                $null = & git reset --hard $targetVersion 2>&1
            }
            
            # Restore error action
            $ErrorActionPreference = $oldErrorAction
            
            Write-Host "  Repository updated successfully" -ForegroundColor Green
        } catch {
            $ErrorActionPreference = $oldErrorAction
            throw
        } finally {
            Pop-Location
        }
    } else {
        # Fresh clone
        Write-Host "  Cloning fresh repository..." -ForegroundColor Green
        Write-Host ""
        
        # Ensure upstream directory exists
        $upstreamParent = Split-Path -Parent $W90UpstreamDir
        if (-not (Test-Path $upstreamParent)) {
            New-Item -ItemType Directory -Path $upstreamParent -Force | Out-Null
        }
        
        # Remove upstream/wannier90 if it exists but isn't a git repo
        if (Test-Path $W90UpstreamDir) {
            Remove-Item -Path $W90UpstreamDir -Recurse -Force
        }
        
        # Clone repository
        # Temporarily disable error action to handle git's stderr output
        $oldErrorAction = $ErrorActionPreference
        $ErrorActionPreference = 'Continue'
        
        try {
            $null = & git clone --branch $targetVersion --depth 1 $w90RepoUrl $W90UpstreamDir 2>&1
            if ($LASTEXITCODE -ne 0) {
                # If branch doesn't exist, try tag
                Write-Host "  Branch not found, trying tag..." -ForegroundColor Yellow
                $null = & git clone --depth 1 $w90RepoUrl $W90UpstreamDir 2>&1
                if ($LASTEXITCODE -eq 0) {
                    Push-Location $W90UpstreamDir
                    try {
                        $null = & git checkout $targetVersion 2>&1
                        if ($LASTEXITCODE -ne 0) {
                            throw "Failed to checkout $targetVersion"
                        }
                    } finally {
                        Pop-Location
                    }
                } else {
                    throw "Failed to clone Wannier90 repository"
                }
            }
        } finally {
            $ErrorActionPreference = $oldErrorAction
        }
        
        Write-Host "  Repository cloned successfully" -ForegroundColor Green
    }
    
    Write-Host ""
    Write-Host "==========================================" -ForegroundColor Green
    Write-Host "Wannier90 source refresh completed!" -ForegroundColor Green
    Write-Host "==========================================" -ForegroundColor Green
    Write-Host ""
    Write-Host "Source location: $W90UpstreamDir" -ForegroundColor Cyan
    Write-Host ""
} else {
    Write-Host "DRY RUN: Would clone/update Wannier90 repository" -ForegroundColor Cyan
    Write-Host ""
}

