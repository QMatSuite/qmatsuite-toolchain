# Quantum ESPRESSO Source Refresh Script
# Safely cleans QE upstream sources and clones fresh QE from git
# Usage: .\scripts\refresh_qe_source.ps1 [-DryRun] [-QeVersion <tag>]

param(
    [switch]$DryRun,           # If specified, only show what would be deleted
    [string]$QeVersion = $null,  # Specific QE version tag to use (e.g., "qe-7.5"). If not specified, uses newest release
    [switch]$NoPatch           # If specified, skip applying local QE patches
)

$ErrorActionPreference = 'Stop'

# Determine script location and repo root
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$RepoRoot = Resolve-Path (Join-Path $ScriptDir "..")
$QeUpstreamDir = Join-Path $RepoRoot "upstream\qe"

Write-Host "==========================================" -ForegroundColor Cyan
Write-Host "Quantum ESPRESSO Source Refresh" -ForegroundColor Cyan
Write-Host "==========================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "Repository root: $RepoRoot"
Write-Host "QE upstream directory: $QeUpstreamDir"
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
    "LICENSE.txt"
)

# Define QE source directories to remove (known QE components)
$QeSourceDirs = @(
    "PW", "EPW", "PHonon", "Modules", "upflib", "atomic", "CPV", "COUPLE",
    "dft-d3", "FFTXlib", "GUI", "GWW", "HP", "KCW", "KS_Solvers",
    "LAXlib", "LR_Modules", "NEB", "PIOUD", "PP", "PWCOND", "QEHeat",
    "TDDFPT", "UtilXlib", "XClib", "XSpectra", "external", "cmake",
    "install", "include", "Doc", "test-suite", "dev-tools", "pseudo"
)

# Define QE config files to remove
$QeConfigFiles = @(
    "CMakeLists.txt",
    "configure",
    "Makefile",
    "environment_variables",
    "License",
    "logo.jpg",
    "README.md",
    "CONTRIBUTING.md",
    "README_GPU.md"
)

# Define build artifact patterns
$BuildArtifactPatterns = @(
    "build*",
    "cmake-build-*",
    "CMakeCache.txt",
    "CMakeFiles",
    "install",
    "bin",
    "lib",
    "*.obj",
    "*.mod",
    "*.a",
    "*.dll",
    "*.exe",
    "*.pdb",
    "*.ilk",
    "*.ninja*",
    "compile_commands.json",
    "qeConfig.cmake",
    "qeConfigVersion.cmake",
    "quantum_espresso.pc"
)

# Step 1: Cleanup QE sources and build artifacts
Write-Host "Step 1: Cleaning up QE sources and build artifacts..." -ForegroundColor Yellow
Write-Host ""

$itemsToDelete = @()
$warnings = @()

# Collect QE source directories
foreach ($dir in $QeSourceDirs) {
    $dirPath = Join-Path $RepoRoot $dir
    if (Test-Path $dirPath) {
        $itemsToDelete += @{
            Path = $dirPath
            Type = "Directory"
            Reason = "QE source directory"
        }
    }
}

# Collect QE config files
foreach ($file in $QeConfigFiles) {
    $filePath = Join-Path $RepoRoot $file
    if (Test-Path $filePath) {
        # Check if it's in allowlist (e.g., README.md might be ours)
        $isAllowed = $false
        if ($file -eq "README.md" -or $file -eq "LICENSE" -or $file -eq "LICENSE.txt") {
            # Be conservative: if README.md or LICENSE exists, check if it's likely ours
            # For now, we'll delete QE's README.md but warn
            if ($file -eq "README.md") {
                $warnings += "README.md will be deleted. If this is your project README, restore it from git."
            }
        }
        if (-not $isAllowed) {
            $itemsToDelete += @{
                Path = $filePath
                Type = "File"
                Reason = "QE config file"
            }
        }
    }
}

# Collect build artifacts
# Check for build directories
$buildDirs = Get-ChildItem -Path $RepoRoot -Directory -Filter "build*" -ErrorAction SilentlyContinue
foreach ($dir in $buildDirs) {
    if ($dir.Name -notlike "upstream*") {
        $itemsToDelete += @{
            Path = $dir.FullName
            Type = "Directory"
            Reason = "Build directory"
        }
    }
}

# Check for cmake-build directories
$cmakeBuildDirs = Get-ChildItem -Path $RepoRoot -Directory -Filter "cmake-build-*" -ErrorAction SilentlyContinue
foreach ($dir in $cmakeBuildDirs) {
    $itemsToDelete += @{
        Path = $dir.FullName
        Type = "Directory"
        Reason = "CMake build directory"
    }
}

# Check for CMakeCache.txt in root
$cmakeCache = Join-Path $RepoRoot "CMakeCache.txt"
if (Test-Path $cmakeCache) {
    $itemsToDelete += @{
        Path = $cmakeCache
        Type = "File"
        Reason = "CMake cache file"
    }
}

# Check for CMakeFiles in root
$cmakeFiles = Join-Path $RepoRoot "CMakeFiles"
if (Test-Path $cmakeFiles) {
    $itemsToDelete += @{
        Path = $cmakeFiles
        Type = "Directory"
        Reason = "CMake files directory"
    }
}

# Check for other build artifact files in root
$artifactFiles = @("compile_commands.json", "qeConfig.cmake", "qeConfigVersion.cmake", "quantum_espresso.pc")
foreach ($file in $artifactFiles) {
    $filePath = Join-Path $RepoRoot $file
    if (Test-Path $filePath) {
        $itemsToDelete += @{
            Path = $filePath
            Type = "File"
            Reason = "Build artifact"
        }
    }
}

# Display what will be deleted
if ($itemsToDelete.Count -gt 0) {
    Write-Host "Items to be deleted:" -ForegroundColor Yellow
    foreach ($item in $itemsToDelete) {
        Write-Host "  [$($item.Type)] $($item.Path) - $($item.Reason)" -ForegroundColor Gray
    }
    Write-Host ""
    
    if ($warnings.Count -gt 0) {
        Write-Host "Warnings:" -ForegroundColor Yellow
        foreach ($warning in $warnings) {
            Write-Host "  $warning" -ForegroundColor Yellow
        }
        Write-Host ""
    }
    
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
    Write-Host "No QE sources or build artifacts found to clean." -ForegroundColor Green
    Write-Host ""
}

# Step 2: Determine QE version to clone
Write-Host "Step 2: Determining QE version to clone..." -ForegroundColor Yellow

$qeRepoUrl = "https://github.com/QEF/q-e.git"
$targetVersion = $null

if ($QeVersion) {
    $targetVersion = $QeVersion
    Write-Host "  Using specified version: $targetVersion" -ForegroundColor Green
} else {
    Write-Host "  Fetching available release tags from $qeRepoUrl..." -ForegroundColor Gray
    
    try {
        # Get all tags from remote
        $tagsOutput = git ls-remote --tags $qeRepoUrl 2>&1
        if ($LASTEXITCODE -ne 0) {
            throw "Failed to fetch tags: $tagsOutput"
        }
        
        # Parse tags and find newest release (qe-X.Y or qe-X.Y.Z format, ignore RC tags)
        $releaseTags = @()
        foreach ($line in $tagsOutput) {
            if ($line -match 'refs/tags/(qe-\d+\.\d+(?:\.\d+)?)$') {
                $tagName = $matches[1]
                if ($tagName -notmatch 'RC|rc|alpha|beta|dev') {
                    $releaseTags += $tagName
                }
            }
        }
        
        if ($releaseTags.Count -eq 0) {
            throw "No release tags found matching pattern qe-X.Y or qe-X.Y.Z"
        }
        
        # Sort by version (simple string sort should work for semantic versions)
        $releaseTags = $releaseTags | Sort-Object -Descending
        $targetVersion = $releaseTags[0]
        
        Write-Host "  Found newest release: $targetVersion" -ForegroundColor Green
    } catch {
        Write-Host "  ERROR: Failed to determine newest release: $($_.Exception.Message)" -ForegroundColor Red
        Write-Host "  Falling back to 'master' branch" -ForegroundColor Yellow
        $targetVersion = "master"
    }
}
Write-Host ""

# Step 3: Clone or update QE repository
if (-not $DryRun) {
    Write-Host "Step 3: Updating QE repository..." -ForegroundColor Yellow
    Write-Host "  Repository: $qeRepoUrl" -ForegroundColor Gray
    Write-Host "  Version: $targetVersion" -ForegroundColor Gray
    Write-Host "  Destination: $QeUpstreamDir" -ForegroundColor Gray
    Write-Host ""
    
    # Check if upstream/qe is already a git repository
    $isGitRepo = Test-Path (Join-Path $QeUpstreamDir ".git")
    
    if ($isGitRepo) {
        # Reuse existing clone: fetch, checkout, update submodules
        Write-Host "  Reusing existing git repository..." -ForegroundColor Green
        Write-Host ""
        
        Push-Location $QeUpstreamDir
        try {
            # Temporarily disable error action for git commands
            $oldErrorAction = $ErrorActionPreference
            $ErrorActionPreference = 'Continue'
            
            # Check and update remote URL if needed
            Write-Host "  Checking remote URL..." -ForegroundColor Gray
            $currentRemote = git remote get-url origin 2>&1
            if ($LASTEXITCODE -ne 0) {
                Write-Host "    Setting remote URL to $qeRepoUrl" -ForegroundColor Yellow
                git remote set-url origin $qeRepoUrl 2>&1 | Out-Null
            } elseif ($currentRemote -ne $qeRepoUrl) {
                Write-Host "    Remote URL differs. Updating from '$currentRemote' to '$qeRepoUrl'" -ForegroundColor Yellow
                git remote set-url origin $qeRepoUrl 2>&1 | Out-Null
            } else {
                Write-Host "    Remote URL is correct" -ForegroundColor Green
            }
            Write-Host ""
            
            # Fetch tags
            Write-Host "  Fetching tags..." -ForegroundColor Gray
            $fetchOutput = git fetch --tags --prune origin 2>&1
            $fetchOutput | ForEach-Object {
                Write-Host "    $_" -ForegroundColor Gray
            }
            if ($LASTEXITCODE -ne 0) {
                Write-Host "    WARNING: git fetch had issues (exit code $LASTEXITCODE)" -ForegroundColor Yellow
            }
            Write-Host ""
            
            # Determine newest tag from local tags
            Write-Host "  Determining newest release tag from local tags..." -ForegroundColor Gray
            $localTags = git tag --list 'qe-*' 2>&1
            if ($LASTEXITCODE -eq 0 -and $localTags.Count -gt 0) {
                # Filter out RC/pre-release tags and sort
                $releaseTags = @()
                foreach ($tag in $localTags) {
                    if ($tag -match '^qe-\d+\.\d+(?:\.\d+)?$' -and $tag -notmatch 'RC|rc|alpha|beta|dev') {
                        $releaseTags += $tag
                    }
                }
                
                if ($releaseTags.Count -gt 0) {
                    $releaseTags = $releaseTags | Sort-Object -Descending
                    $newestLocalTag = $releaseTags[0]
                    Write-Host "    Found newest local tag: $newestLocalTag" -ForegroundColor Green
                    
                    # Use specified version if provided, otherwise use newest local tag
                    if ($QeVersion) {
                        $checkoutTag = $targetVersion
                    } else {
                        $checkoutTag = $newestLocalTag
                    }
                } else {
                    Write-Host "    No release tags found locally, using specified/default: $targetVersion" -ForegroundColor Yellow
                    $checkoutTag = $targetVersion
                }
            } else {
                Write-Host "    No local tags found, using specified/default: $targetVersion" -ForegroundColor Yellow
                $checkoutTag = $targetVersion
            }
            Write-Host ""
            
            # Checkout the target tag
            Write-Host "  Checking out tag: $checkoutTag" -ForegroundColor Gray
            $checkoutOutput = git checkout $checkoutTag 2>&1
            $checkoutOutput | ForEach-Object {
                Write-Host "    $_" -ForegroundColor Gray
            }
            if ($LASTEXITCODE -ne 0) {
                throw "git checkout failed with exit code $LASTEXITCODE"
            }
            Write-Host ""
            
            # Sync and update submodules
            Write-Host "  Syncing submodules..." -ForegroundColor Gray
            $syncOutput = git submodule sync --recursive 2>&1
            $syncOutput | ForEach-Object {
                Write-Host "    $_" -ForegroundColor Gray
            }
            if ($LASTEXITCODE -ne 0) {
                Write-Host "    WARNING: git submodule sync had issues (exit code $LASTEXITCODE)" -ForegroundColor Yellow
            }
            Write-Host ""
            
            Write-Host "  Updating submodules..." -ForegroundColor Gray
            $submoduleOutput = git submodule update --init --recursive 2>&1
            $submoduleOutput | ForEach-Object {
                Write-Host "    $_" -ForegroundColor Gray
            }
            if ($LASTEXITCODE -ne 0) {
                Write-Host "    WARNING: Submodule update had issues (exit code $LASTEXITCODE)" -ForegroundColor Yellow
            } else {
                Write-Host "    Submodules updated successfully." -ForegroundColor Green
            }
            Write-Host ""
            
            # Restore error action
            $ErrorActionPreference = $oldErrorAction
            
            Write-Host "  Repository update completed successfully." -ForegroundColor Green
            Write-Host ""
            
        } catch {
            # Restore error action before throwing
            $ErrorActionPreference = $oldErrorAction
            Write-Host "  ERROR: Failed to update repository: $($_.Exception.Message)" -ForegroundColor Red
            throw
        } finally {
            Pop-Location
        }
        
    } else {
        # Fresh clone needed
        Write-Host "  No existing git repository found. Cloning fresh..." -ForegroundColor Yellow
        Write-Host ""
        
        # If upstream/qe exists but isn't a git repo, rename it to backup
        if (Test-Path $QeUpstreamDir) {
            $timestamp = Get-Date -Format "yyyyMMdd-HHmmss"
            $backupPath = "$QeUpstreamDir.bak.$timestamp"
            Write-Host "  Existing directory found but not a git repo. Renaming to: $backupPath" -ForegroundColor Yellow
            Move-Item -Path $QeUpstreamDir -Destination $backupPath -Force -ErrorAction Stop
            Write-Host ""
        }
        
        # Create upstream directory
        $upstreamParent = Split-Path -Parent $QeUpstreamDir
        if (-not (Test-Path $upstreamParent)) {
            New-Item -ItemType Directory -Path $upstreamParent -Force | Out-Null
        }
        
        # Clone with recursive submodules
        Write-Host "  Cloning repository (this may take a few minutes)..." -ForegroundColor Yellow
        try {
            # Temporarily disable error action to allow git's informational messages
            $oldErrorAction = $ErrorActionPreference
            $ErrorActionPreference = 'Continue'
            
            # Capture git output (both stdout and stderr)
            $gitOutput = git clone --recursive --branch $targetVersion --depth 1 $qeRepoUrl $QeUpstreamDir 2>&1
            $exitCode = $LASTEXITCODE
            
            # Restore error action
            $ErrorActionPreference = $oldErrorAction
            
            # Display git output (git writes progress to stderr, which PowerShell may interpret as errors)
            $gitOutput | ForEach-Object {
                Write-Host "    $_" -ForegroundColor Gray
            }
            
            if ($exitCode -ne 0) {
                throw "git clone failed with exit code $exitCode"
            }
            
            Write-Host "  Clone completed successfully." -ForegroundColor Green
            Write-Host ""
            
            # Ensure submodules are initialized (in case --recursive didn't work fully)
            Write-Host "  Initializing submodules..." -ForegroundColor Yellow
            Push-Location $QeUpstreamDir
            try {
                # Temporarily disable error action for git submodule output
                $oldErrorAction = $ErrorActionPreference
                $ErrorActionPreference = 'Continue'
                
                $submoduleOutput = git submodule update --init --recursive 2>&1
                $submoduleExitCode = $LASTEXITCODE
                
                # Restore error action
                $ErrorActionPreference = $oldErrorAction
                
                # Display submodule output
                $submoduleOutput | ForEach-Object {
                    Write-Host "    $_" -ForegroundColor Gray
                }
                
                if ($submoduleExitCode -eq 0) {
                    Write-Host "  Submodules initialized successfully." -ForegroundColor Green
                } else {
                    Write-Host "  WARNING: Submodule initialization had issues (exit code $submoduleExitCode)" -ForegroundColor Yellow
                }
            } finally {
                Pop-Location
            }
            Write-Host ""
            
        } catch {
            Write-Host "  ERROR: Failed to clone QE repository: $($_.Exception.Message)" -ForegroundColor Red
            throw
        }
    }
} else {
    # Dry run mode
    $isGitRepo = Test-Path (Join-Path $QeUpstreamDir ".git")
    if ($isGitRepo) {
        Write-Host "Step 3: Would update existing git repository..." -ForegroundColor Cyan
        Write-Host "  Repository: $qeRepoUrl" -ForegroundColor Gray
        Write-Host "  Version: $targetVersion" -ForegroundColor Gray
        Write-Host "  Would run: git fetch --tags, git checkout $targetVersion, git submodule sync/update" -ForegroundColor Gray
    } else {
        Write-Host "Step 3: Would clone QE repository..." -ForegroundColor Cyan
        Write-Host "  Repository: $qeRepoUrl" -ForegroundColor Gray
        Write-Host "  Version: $targetVersion" -ForegroundColor Gray
        Write-Host "  Destination: $QeUpstreamDir" -ForegroundColor Gray
    }
    Write-Host ""
}

# Step 4: Verify submodules
if (-not $DryRun -and (Test-Path $QeUpstreamDir)) {
    Write-Host "Step 4: Verifying submodules..." -ForegroundColor Yellow
    
    $requiredSubmodules = @("external/wannier90", "external/mbd", "external/devxlib")
    $missingSubmodules = @()
    
    foreach ($submodule in $requiredSubmodules) {
        $submodulePath = Join-Path $QeUpstreamDir $submodule
        if (-not (Test-Path $submodulePath)) {
            $missingSubmodules += $submodule
            Write-Host "  MISSING: $submodule" -ForegroundColor Red
        } else {
            # Check if it has content (not just empty directory)
            $hasContent = (Get-ChildItem $submodulePath -Recurse -File -ErrorAction SilentlyContinue | Measure-Object).Count -gt 0
            if (-not $hasContent) {
                $missingSubmodules += $submodule
                Write-Host "  EMPTY: $submodule" -ForegroundColor Yellow
            } else {
                Write-Host "  OK: $submodule" -ForegroundColor Green
            }
        }
    }
    
    if ($missingSubmodules.Count -gt 0) {
        Write-Host ""
        Write-Host "WARNING: Some submodules are missing or empty!" -ForegroundColor Yellow
        Write-Host "  Missing: $($missingSubmodules -join ', ')" -ForegroundColor Yellow
        Write-Host "  Try running: cd upstream/qe && git submodule update --init --recursive" -ForegroundColor Cyan
        Write-Host ""
    } else {
        Write-Host "  All required submodules are present." -ForegroundColor Green
    }
    Write-Host ""
}

# Step 5: Apply local QE patches unless skipped
if (-not $DryRun -and -not $NoPatch) {
    Write-Host "Step 5: Applying local QE patches..." -ForegroundColor Yellow
    $applyScript = Join-Path $RepoRoot "scripts\\apply_qe_patches.ps1"
    $patchDir = Join-Path $RepoRoot "patches"
    if (-not (Test-Path $applyScript)) {
        Write-Host "  ERROR: Patch apply script not found at $applyScript" -ForegroundColor Red
        throw "Missing apply_qe_patches.ps1"
    }
    if (-not (Test-Path $patchDir)) {
        Write-Host "  ERROR: Patch directory not found at $patchDir" -ForegroundColor Red
        throw "Missing patches directory"
    }
    & powershell -ExecutionPolicy Bypass -File $applyScript -QeDir $QeUpstreamDir -PatchDir $patchDir
    if ($LASTEXITCODE -eq 0) {
        Write-Host "  Patches applied successfully." -ForegroundColor Green
    } else {
        throw "Patch application failed (exit $LASTEXITCODE)"
    }
    Write-Host ""
} elseif ($NoPatch) {
    Write-Host "Step 5: Skipped patch application (-NoPatch set)." -ForegroundColor Yellow
    Write-Host ""
}

# Summary
Write-Host "==========================================" -ForegroundColor Green
if ($DryRun) {
    Write-Host "Dry run completed. Use without -DryRun to perform cleanup and clone." -ForegroundColor Cyan
} else {
    Write-Host "QE source refresh completed!" -ForegroundColor Green
    Write-Host ""
    Write-Host "Next steps:" -ForegroundColor Cyan
    Write-Host "  1. QE source is now at: $QeUpstreamDir" -ForegroundColor White
    Write-Host "  2. Run the build script:" -ForegroundColor White
    Write-Host "     .\scripts\build_qe_win_oneapi.ps1 [-NoMpi]" -ForegroundColor Yellow
    Write-Host ""
}
Write-Host "==========================================" -ForegroundColor Green

exit 0

