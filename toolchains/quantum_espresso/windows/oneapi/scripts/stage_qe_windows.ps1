<#
.SYNOPSIS
Stages Windows oneAPI QE executables and their runtime DLLs into a dist folder.

.PARAMETER BuildDir
Path to the QE build directory (default: upstream/qe/build-win-oneapi).

.PARAMETER DistDir
Target dist directory (default: dist/win-oneapi).

.PARAMETER Only
Optional list of executable basenames (with or without .exe) to stage.

.PARAMETER Clean
If set, deletes the dist directory before staging.

.PARAMETER VerboseDeps
Print parsed dependencies for each executable/DLL when using dumpbin.

.PARAMETER WhatIf
Dry run: print actions without copying/deleting.

.PARAMETER ProbeExe
Executable name to probe after staging (default: pw.x.exe if present, else first exe).

.PARAMETER NoStrict
If set, disables strict mode (not recommended). By default, staging fails if smoke test exit code != 0 OR "JOB DONE" not detected.

.PARAMETER NoCleanSmokeTest
If set, skips cleanup of smoke test files (outdir/, *.out, *.in, *.UPF) after smoke test completes. By default, these files are cleaned up.
#>
[CmdletBinding(SupportsShouldProcess = $true)]
param(
    [string]$BuildDir = "upstream/qe/build-win-oneapi",
    [string]$DistDir  = "dist/win-oneapi",
    [string[]]$Only,
    [switch]$Clean,
    [switch]$VerboseDeps,
    [string]$ProbeExe,
    [switch]$NoStrict,  # If set, disables strict mode (not recommended)
    [switch]$NoCleanSmokeTest  # If set, skips cleanup of smoke test files (outdir, *.out, *.in, *.UPF). By default, cleans these files.
)

$ErrorActionPreference = 'Stop'

# Validate parameters
if ([string]::IsNullOrWhiteSpace($BuildDir)) {
    throw "BuildDir parameter cannot be null or empty."
}
if ([string]::IsNullOrWhiteSpace($DistDir)) {
    throw "DistDir parameter cannot be null or empty."
}

#---------------- Helpers ----------------#
function Get-DumpbinPath {
    # Prefer VS toolchain location if present
    if ($env:VCToolsInstallDir) {
        $candidate = Join-Path $env:VCToolsInstallDir "bin\Hostx64\x64\dumpbin.exe"
        if (Test-Path $candidate) { return $candidate }
    }
    $cmd = Get-Command dumpbin.exe -ErrorAction SilentlyContinue
    if ($cmd) { return $cmd.Source }
    return $null
}

# Locate VS2022 dumpbin and add it to PATH for this session if not already present
function Ensure-DumpbinInPath {
    $path = Get-DumpbinPath
    if ($path) { return $path }

    $editions = @("Community","Professional","Enterprise","BuildTools")
    foreach ($ed in $editions) {
        $msvcRoot = "C:\Program Files\Microsoft Visual Studio\2022\$ed\VC\Tools\MSVC"
        if (-not (Test-Path $msvcRoot)) { continue }
        $versions = Get-ChildItem -Path $msvcRoot -Directory -ErrorAction SilentlyContinue | Sort-Object Name -Descending
        foreach ($v in $versions) {
            $candidate = Join-Path $v.FullName "bin\Hostx64\x64\dumpbin.exe"
            if (Test-Path $candidate) {
                # Prepend to PATH for this session so subsequent calls can resolve it
                $dir = Split-Path $candidate -Parent
                $env:PATH = "$dir;$env:PATH"
                # Also set VCToolsInstallDir for compatibility with existing logic
                $env:VCToolsInstallDir = $v.FullName
                return $candidate
            }
        }
    }
    return $null
}

function Get-SearchRoots {
    param(
        [string]$OneApiRoot
    )
    $roots = @()
    if ($OneApiRoot -and (Test-Path $OneApiRoot)) {
        # versioned top-level bins, lib, redist
        $roots += Get-ChildItem -Path $OneApiRoot -Directory -ErrorAction SilentlyContinue | ForEach-Object {
            @(
                (Join-Path $_.FullName "bin"),
                (Join-Path $_.FullName "lib"),
                (Join-Path $_.FullName "redist")
            )
        }
        # compiler tree
        $compilerBase = Join-Path $OneApiRoot "compiler"
        if (Test-Path $compilerBase) {
            $compilerDirs = Get-ChildItem -Path $compilerBase -Directory -ErrorAction SilentlyContinue
            foreach ($dir in $compilerDirs) {
                $roots += @(
                    (Join-Path $dir.FullName "bin"),
                    (Join-Path $dir.FullName "lib"),
                    (Join-Path $dir.FullName "redist\intel64\compiler"),
                    (Join-Path $dir.FullName "redist\intel64"),
                    (Join-Path $dir.FullName "windows\redist\intel64")
                )
            }
        }
        # mkl tree
        $mklBase = Join-Path $OneApiRoot "mkl"
        if (Test-Path $mklBase) {
            $mklDirs = Get-ChildItem -Path $mklBase -Directory -ErrorAction SilentlyContinue
            foreach ($dir in $mklDirs) {
                $roots += @(
                    (Join-Path $dir.FullName "bin"),
                    (Join-Path $dir.FullName "lib"),
                    (Join-Path $dir.FullName "redist\intel64"),
                    (Join-Path $dir.FullName "redist\intel64\mkl")
                )
            }
        }
    }
    $roots = $roots | Where-Object { $_ -and (Test-Path $_) }
    return ($roots | Select-Object -Unique)
}

function Get-LibxcDependentExecutables {
    param(
        [string]$DumpbinExe,
        [string[]]$ExePaths
    )

    $hits = @()

    foreach ($exe in $ExePaths) {
        if (-not (Test-Path $exe)) { continue }

        $out = & $DumpbinExe /DEPENDENTS $exe 2>$null | Out-String
        if ($out -match '(?i)\bxc\.dll\b' -or $out -match '(?i)\bxcf03\.dll\b') {
            $hits += $exe
        }
    }

    return $hits
}
function Get-ExeDependentsDumpbin {
    param(
        [string]$Path,
        [string]$DumpbinPath,
        [switch]$Verbose,
        [int]$TimeoutMs = 5000
    )
    $deps = @()
    $pinfo = New-Object System.Diagnostics.ProcessStartInfo
    $pinfo.FileName = $DumpbinPath
    $pinfo.Arguments = "/DEPENDENTS `"$Path`""
    $pinfo.UseShellExecute = $false
    $pinfo.RedirectStandardOutput = $true
    $pinfo.RedirectStandardError = $true
    $proc = New-Object System.Diagnostics.Process
    $proc.StartInfo = $pinfo
    $proc.Start() | Out-Null
    
    # Wait with timeout (5 seconds default)
    if (-not $proc.WaitForExit($TimeoutMs)) {
        Write-Warning "dumpbin timed out after ${TimeoutMs}ms on $Path"
        $proc.Kill()
        $proc.WaitForExit() | Out-Null
        return $deps
    }
    
    $stdout = $proc.StandardOutput.ReadToEnd()
    $stderr = $proc.StandardError.ReadToEnd()
    if ($proc.ExitCode -ne 0) {
        Write-Warning "dumpbin failed on $Path (exit $($proc.ExitCode)): $stderr"
        return $deps
    }
    $capture = $false
    foreach ($line in $stdout -split "`r?`n") {
        if ($line -match "Image has the following dependencies") { 
            $capture = $true
            continue 
        }
        # Stop capturing when we hit the Summary section or another section header
        if ($capture -and ($line -match "^\s*Summary\s*$" -or $line -match "^\s*[A-Z].*:$")) { 
            $capture = $false
            break
        }
        if ($capture -and $line -match "([A-Za-z0-9._-]+\.dll)") {
            $dll = $Matches[1].ToLowerInvariant()
            $deps += $dll
        }
    }
    if ($Verbose) {
        Write-Host ("Dependencies for {0}:" -f $Path) -ForegroundColor Cyan
        ($deps | Sort-Object -Unique) | ForEach-Object { Write-Host "  $_" }
    }
    return ($deps | Select-Object -Unique)
}

function Is-SystemDll {
    param(
        [string]$Name,
        [string]$ResolvedPath
    )
    $systemNames = @(
        "kernel32.dll","user32.dll","advapi32.dll","gdi32.dll","shell32.dll","ole32.dll",
        "oleaut32.dll","ws2_32.dll","winmm.dll","secur32.dll","bcrypt.dll","ntdll.dll",
        "ucrtbase.dll","msvcrt.dll","vcruntime140.dll","vcruntime140_1.dll","concrt140.dll",
        "api-ms-win-"
    )
    foreach ($s in $systemNames) {
        if ($Name.ToLowerInvariant().StartsWith($s.TrimEnd('*'))) { return $true }
    }
    if ($ResolvedPath) {
        $sysRoots = @("$env:WINDIR\System32","$env:WINDIR\SysWOW64") | Where-Object { $_ }
        foreach ($sr in $sysRoots) {
            if ($ResolvedPath -like "$sr*") { return $true }
        }
    }
    return $false
}

function Is-DenylistedDll {
    param(
        [string]$Name
    )
    # Denylist: these DLLs must NEVER be staged
    $denylistPatterns = @(
        "*_db.dll",      # Debug runtimes (e.g., libiomp5md_db.dll)
        "libimalloc.dll" # Intel scalable allocator (not required)
    )
    foreach ($pattern in $denylistPatterns) {
        if ($Name -like $pattern) {
            return $true
        }
    }
    return $false
}

function Resolve-DllPath {
    param(
        [string]$Name,
        [string[]]$SearchRoots,
        [string[]]$PathDirs  # Kept for compatibility but should be empty for approved roots only
    )
    # Search approved roots first
    foreach ($root in $SearchRoots) {
        if (-not $root -or -not (Test-Path $root)) { continue }
        $found = Get-ChildItem -Path $root -Filter $Name -File -Recurse -ErrorAction SilentlyContinue | Select-Object -First 1
        if ($found) { return $found.FullName }
    }
    # Only search PathDirs if explicitly provided (for backward compatibility, but should be empty)
    foreach ($dir in $PathDirs) {
        if (-not $dir) { continue }
        $candidate = Join-Path $dir $Name
        if (Test-Path $candidate) { return (Resolve-Path $candidate).Path }
    }
    return $null
}

function Stage-Executables {
    param(
        [string]$BinDir,
        [string]$DistBinDir,  # Changed: now expects dist/bin directory
        [string[]]$Only
    )
    $exeFiles = Get-ChildItem -Path $BinDir -Filter *.exe -File
    if ($Only) {
        $normalized = $Only | ForEach-Object { if ($_ -like "*.exe") { $_ } else { "$_.exe" } }
        $exeFiles = $exeFiles | Where-Object { $normalized -contains $_.Name }
        Write-Host "Filtering executables (Only mode): $($normalized -join ', ')" -ForegroundColor Cyan
    } else {
        Write-Host "Staging ALL executables (no filter)" -ForegroundColor Cyan
    }
    if (-not $exeFiles) {
        throw "No executables found to stage. Check build output or -Only filters. BinDir: $BinDir"
    }
    
    $exeCount = ($exeFiles | Measure-Object).Count
    Write-Host "Found $exeCount executable(s) to stage" -ForegroundColor Green
    
    # Ensure dist/bin directory exists
    if (-not $WhatIfPreference -and -not (Test-Path $DistBinDir)) {
        New-Item -ItemType Directory -Force -Path $DistBinDir | Out-Null
    }
    
    # Separate .x.exe files from regular .exe files
    $xExeFiles = $exeFiles | Where-Object { $_.Name -like "*.x.exe" }
    $regularExeFiles = $exeFiles | Where-Object { $_.Name -notlike "*.x.exe" }
    
    # Process .x.exe files: remove existing .exe, copy .x.exe, then rename to .exe
    foreach ($xExe in $xExeFiles) {
        $targetName = $xExe.Name -replace '\.x\.exe$', '.exe'
        $targetPath = Join-Path $DistBinDir $targetName
        $xExeDest = Join-Path $DistBinDir $xExe.Name
        
        if ($PSCmdlet.ShouldProcess($targetPath, "Stage $($xExe.Name) -> $targetName") -and -not $WhatIfPreference) {
            # Remove existing .exe if it exists
            if (Test-Path $targetPath) {
                Remove-Item -Path $targetPath -Force -ErrorAction SilentlyContinue
            }
            # Remove any existing .x.exe in dist (shouldn't exist, but clean up just in case)
            if (Test-Path $xExeDest) {
                Remove-Item -Path $xExeDest -Force -ErrorAction SilentlyContinue
            }
            # Copy .x.exe to dist/bin
            Copy-Item -Force $xExe.FullName $xExeDest
            # Rename to .exe (removes .x.exe from dist)
            Rename-Item -Path $xExeDest -NewName $targetName -Force
        }
        Write-Host "  exe: $($xExe.Name) -> $targetName"
    }
    
    # Process regular .exe files (not .x.exe): just copy them
    foreach ($exe in $regularExeFiles) {
        $dest = Join-Path $DistBinDir $exe.Name
        if ($PSCmdlet.ShouldProcess($dest, "Copy exe $($exe.FullName)") -and -not $WhatIfPreference) {
            Copy-Item -Force $exe.FullName $dest
        }
        Write-Host "  exe: $($exe.Name)"
    }
    
    # Print summary with examples
    $stagedCount = ($xExeFiles | Measure-Object).Count + ($regularExeFiles | Measure-Object).Count
    Write-Host ""
    Write-Host "Staged $stagedCount executable(s) total" -ForegroundColor Green
    if ($stagedCount -gt 0) {
        Write-Host "Examples of staged executables:" -ForegroundColor Cyan
        $examples = $exeFiles | Select-Object -First 5
        foreach ($ex in $examples) {
            $finalName = if ($ex.Name -like "*.x.exe") { $ex.Name -replace '\.x\.exe$', '.exe' } else { $ex.Name }
            Write-Host "  - $finalName" -ForegroundColor Gray
        }
        if ($stagedCount -gt 5) {
            Write-Host "  ... and $($stagedCount - 5) more" -ForegroundColor Gray
        }
    }
    Write-Host ""
    
    return $exeFiles
}

function Copy-DllPatterns {
    param(
        [string[]]$Patterns,
        [string[]]$SearchRoots,
        [string]$DistDir,
        [switch]$Optional
    )
    foreach ($pat in $Patterns) {
        $found = $null
        foreach ($root in $SearchRoots) {
            $found = Get-ChildItem -Path $root -Filter $pat -File -Recurse -ErrorAction SilentlyContinue | Select-Object -First 1
            if ($found) { break }
        }
        if ($found) {
            $dest = Join-Path $DistDir $found.Name
            if ($PSCmdlet.ShouldProcess($dest, "Copy dll $($found.FullName)") -and -not $WhatIfPreference) {
                Copy-Item -Force $found.FullName $dest
            }
            Write-Host "  dll: $($found.Name) (from $($found.DirectoryName))"
        } else {
            if ($Optional) {
                Write-Host "Optional pattern '$pat' not found; will rely on dependency scan / runtime probe."
            } else {
                Write-Warning "Missing DLL pattern '$pat' (not found under oneAPI search roots)."
            }
        }
    }
}

function Get-ValidPathDirs {
    $dirs = @()
    foreach ($p in ($env:PATH -split ';')) {
        if (-not $p) { continue }
        $t = $p.Trim('"')
        try {
            if (Test-Path $t) { $dirs += $t }
        } catch {
            # skip malformed entries
        }
    }
    return $dirs
}

function Find-BuildDirectory {
    param(
        [string]$RepoRoot,
        [string]$DefaultBuildDir
    )
    # If BuildDir was explicitly provided and exists, use it
    $candidate = Join-Path $RepoRoot $DefaultBuildDir
    if (Test-Path $candidate) {
        return $DefaultBuildDir
    }
    
    # If default doesn't exist, search for build* directories in upstream/qe
    $qeDir = Join-Path $RepoRoot "upstream\qe"
    if (Test-Path $qeDir) {
        $buildDirs = Get-ChildItem -Path $qeDir -Directory -Filter "build*" -ErrorAction SilentlyContinue
        if ($buildDirs) {
            # Prefer build-win-oneapi-msmpi if it exists (MPI build), otherwise use first match
            $msmpiDir = $buildDirs | Where-Object { $_.Name -like "*msmpi*" } | Select-Object -First 1
            if ($msmpiDir) {
                $relativePath = $msmpiDir.FullName.Substring($RepoRoot.Length + 1)
                return $relativePath -replace '\\', '/'
            }
            # Otherwise use first build* directory found
            $firstDir = $buildDirs | Select-Object -First 1
            $relativePath = $firstDir.FullName.Substring($RepoRoot.Length + 1)
            return $relativePath -replace '\\', '/'
        }
    }
    
    # Fallback to default (even if it doesn't exist - will error later)
    return $DefaultBuildDir
}

function Get-CMakeCacheBool {
    param(
        [Parameter(Mandatory=$true)][string]$CachePath,
        [Parameter(Mandatory=$true)][string]$Key
    )
    if (-not (Test-Path $CachePath)) { return $null }
    $m = Select-String -Path $CachePath -Pattern ("^\s*" + [regex]::Escape($Key) + ":BOOL=") -ErrorAction SilentlyContinue | Select-Object -First 1
    if (-not $m) { return $null }
    if ((($m.Line -split "=", 2)[1]).Trim() -eq "ON") {return $true}
    return $false
}

function Detect-MPIFromBuild {
    param(
        [string]$BuildDir,
        [string]$RepoRoot
    )
    $buildPath = Join-Path $RepoRoot $BuildDir
    $cmakeCache = Join-Path $buildPath "CMakeCache.txt"
    
    if (-not (Test-Path $cmakeCache)) {
        return $null
    }
    
    # Check for QE_ENABLE_MPI or MPI_FOUND in CMakeCache.txt
    $mpiEnabled = Select-String -Path $cmakeCache -Pattern "^QE_ENABLE_MPI:BOOL=(ON|TRUE)" -CaseSensitive:$false
    $mpiFound = Select-String -Path $cmakeCache -Pattern "^MPI_FOUND:BOOL=TRUE" -CaseSensitive:$false
    
    if ($mpiEnabled -or $mpiFound) {
        # Check if it's MSMPI by looking for MS-MPI references
        $msmpiRef = Select-String -Path $cmakeCache -Pattern "msmpi|MS-MPI|MSMPI" -CaseSensitive:$false
        if ($msmpiRef) {
            return "MSMPI"
        }
        return "MPI"
    }
    
    return "no-MPI"
}

function Enqueue-IfNew {
    param(
        [string]$FilePath,
        [hashtable]$Visited,
        [hashtable]$Enqueued,
        [System.Collections.Queue]$Queue
    )
    $basename = [System.IO.Path]::GetFileName($FilePath).ToLowerInvariant()
    
    # Skip if already visited (scanned)
    if ($Visited.ContainsKey($basename)) {
        return $false
    }
    
    # Skip if already enqueued
    if ($Enqueued.ContainsKey($basename)) {
        return $false
    }
    
    # Add to enqueued and queue
    $Enqueued[$basename] = $true
    $Queue.Enqueue($FilePath)
    return $true
}

function Copy-DllClosure {
    param(
        [string[]]$ExePaths,
        [string]$DumpbinPath,
        [string[]]$SearchRoots,
        [string]$DistBinDir,  # Changed: now expects dist/bin directory
        [string]$BuildBinDir,  # QE build bin directory for approved search roots
        [switch]$VerboseDeps
    )
    $queue = New-Object System.Collections.Queue
    $visited = @{}  # Case-insensitive visited set by basename (already scanned)
    $enqueued = @{}  # Case-insensitive enqueued set by basename (already enqueued)
    $missing = [System.Collections.Generic.List[object]]::new()  # Track unresolved DLLs with parent info
    $copied = [System.Collections.Generic.List[object]]::new()
    $resolved = [System.Collections.Generic.List[object]]::new()  # Track all resolved DLLs for closure file
    $scanCount = 0
    $MAX_SECONDS = 300  # 5 minute wall-clock timeout
    $DUMPBIN_TIMEOUT_MS = 5000
    
    # Build approved search roots: dist/bin, QE build bin, oneAPI redist, MSMPI, VC++ redist
    $approvedRoots = @()
    if (Test-Path $DistBinDir) {
        $approvedRoots += $DistBinDir
    }
    if (Test-Path $BuildBinDir) {
        $approvedRoots += $BuildBinDir
    }
    $approvedRoots += $SearchRoots  # oneAPI redist locations
    # MSMPI and VC++ redist will be added by caller if needed

    Write-Host "Starting recursive dependency closure (max ${MAX_SECONDS}s wall-clock timeout, ${DUMPBIN_TIMEOUT_MS}ms per dumpbin)..."
    Write-Host "  Initial executables to scan: $($ExePaths.Count)"
    Write-Host "  Approved search roots: $($approvedRoots.Count) locations"

    # ---- LIBXC detection (print only) ----
    $cmakeCache = Join-Path $BuildDir "CMakeCache.txt"
    $libxcEnabled = $false

    if (Test-Path $cmakeCache) {
        $line = Select-String -Path $cmakeCache -Pattern '^QE_ENABLE_LIBXC:BOOL=' -ErrorAction SilentlyContinue | Select-Object -First 1
        if ($line -and $line.Line -match 'QE_ENABLE_LIBXC:BOOL=ON') {
            $libxcEnabled = $true
        }
    }

    if ($libxcEnabled) {
        Write-Host ""
        Write-Host "[libxc] QE_ENABLE_LIBXC=ON detected" -ForegroundColor Cyan

        $libxcExes = Get-LibxcDependentExecutables `
            -DumpbinExe $DumpbinPath `
            -ExePaths  $ExePaths
        $libxcExes = $libxcExes | Sort-Object -Unique

        Write-Host "[libxc] Executables linked against xc.dll / xcf03.dll:" -ForegroundColor Cyan
        foreach ($exe in $libxcExes) {
            Write-Host ("[libxc]  " + (Split-Path $exe -Leaf)) -ForegroundColor Yellow
        }

        Write-Host "[libxc] Total: $($libxcExes.Count)" -ForegroundColor Cyan
        Write-Host ""
    }

    # ---- end LIBXC detection ----

    # Initialize BFS by enqueueing all staged exes
    $enqueuedCount = 0
    foreach ($exe in $ExePaths) {
        if (-not (Test-Path $exe)) { 
            Write-Warning "  Executable not found: $exe"
            continue 
        }
        if (Enqueue-IfNew -FilePath $exe -Visited $visited -Enqueued $enqueued -Queue $queue) {
            $enqueuedCount++
        }
    }
    Write-Host "  Enqueued $enqueuedCount unique executable(s) for scanning"

    # Start stopwatch for timeout
    $stopwatch = [System.Diagnostics.Stopwatch]::StartNew()

    while ($queue.Count -gt 0) {
        # Check timeout
        if ($stopwatch.Elapsed.TotalSeconds -gt $MAX_SECONDS) {
            $visitedCount = $visited.Count
            $queueSize = $queue.Count
            $unresolvedSample = $missing | Select-Object -First 10 | ForEach-Object { if ($_ -is [hashtable]) { $_.DllName } else { $_ } }
            throw "Dependency closure exceeded maximum wall-clock time (${MAX_SECONDS}s). Visited: $visitedCount modules, Queue size: $queueSize, Unresolved DLLs (sample): $($unresolvedSample -join ', ')"
        }
        
        $filePath = $queue.Dequeue()
        
        # Get basename for visited check (case-insensitive)
        $basename = [System.IO.Path]::GetFileName($filePath).ToLowerInvariant()
        
        # Skip if already visited (shouldn't happen due to Enqueue-IfNew, but double-check)
        if ($visited.ContainsKey($basename)) { 
            continue 
        }
        
        # Mark as visited and increment scan count (only count unique modules scanned)
        $visited[$basename] = $true
        $scanCount++

        # Get dependencies
        $deps = Get-ExeDependentsDumpbin -Path $filePath -DumpbinPath $DumpbinPath -Verbose:$VerboseDeps -TimeoutMs $DUMPBIN_TIMEOUT_MS
        if ($VerboseDeps) {
            Write-Host "  Scanned: $basename -> $($deps.Count) dependencies"
        }

        foreach ($dllName in $deps) {
            $dllBasename = $dllName.ToLowerInvariant()
            
            # Skip Windows API set pseudo-DLLs (not redistributable files)
            if ($dllBasename -like "api-ms-win-*" -or $dllBasename -like "ext-ms-win-*") {
                continue
            }
            
            # Skip if already visited
            if ($visited.ContainsKey($dllBasename)) { 
                continue 
            }

            # Skip denylisted DLLs
            if (Is-DenylistedDll -Name $dllName) {
                Write-Host "  DENYLIST: Skipping $dllName (prohibited from distribution)"
                continue
            }

            # Resolve DLL path only from approved search roots
            $resolvedPath = Resolve-DllPath -Name $dllName -SearchRoots $approvedRoots -PathDirs @()
            if (-not $resolvedPath) {
                # Record unresolved DLL with parent info
                if (-not (Is-SystemDll -Name $dllName -ResolvedPath $null)) {
                    $missing.Add(@{
                        DllName = $dllName
                        ParentModule = $basename
                    })
                }
                continue
            }
            
            # Skip system DLLs
            if (Is-SystemDll -Name $dllName -ResolvedPath $resolvedPath) { 
                continue 
            }

            # Track all resolved DLLs (for closure file)
            $resolved.Add(@{
                Basename = $dllBasename
                FullPath = $resolvedPath
            })

            # Copy DLL to dist/bin
            $dest = Join-Path $DistBinDir (Split-Path $resolvedPath -Leaf)
            $alreadyExists = (Test-Path $dest)
            if (-not $alreadyExists) {
                if ($PSCmdlet.ShouldProcess($dest, "Copy dll $resolvedPath") -and -not $WhatIfPreference) {
                    Copy-Item -Force $resolvedPath $dest
                }
                $copied.Add(@{
                    Basename = $dllBasename
                    FullPath = $resolvedPath
                })
                Write-Host "  dll: $dllBasename (from $resolvedPath)"
            } else {
                Write-Host "  dll: $dllBasename (already staged)"
            }

            # Enqueue for recursive scanning (using helper to prevent duplicates)
            Enqueue-IfNew -FilePath $resolvedPath -Visited $visited -Enqueued $enqueued -Queue $queue | Out-Null
        }
    }
    
    $stopwatch.Stop()

    # Write deps-closure.txt (include all resolved DLLs, not just newly copied)
    $closureFile = Join-Path (Split-Path $DistBinDir -Parent) "deps-closure.txt"
    if (-not $WhatIfPreference) {
        try {
            if ($resolved.Count -gt 0) {
                $closureLines = @()
                foreach ($item in $resolved) {
                    if ($item -and $item.Basename -and $item.FullPath) {
                        $closureLines += "$($item.Basename)`t$($item.FullPath)"
                    }
                }
                if ($closureLines.Count -gt 0) {
                    $closureLines | Out-File -FilePath $closureFile -Encoding UTF8 -Force
                    Write-Host "  Wrote $($closureLines.Count) entries to deps-closure.txt"
                } else {
                    Write-Warning "resolved list has $($resolved.Count) items but generated no valid lines for deps-closure.txt"
                    "# Dependency closure: resolved list had $($resolved.Count) items but no valid lines generated`n# Scanned $scanCount unique modules in $([math]::Round($stopwatch.Elapsed.TotalSeconds, 2))s." | Out-File -FilePath $closureFile -Encoding UTF8 -Force
                }
            } else {
                # Write file with header comment if no DLLs were resolved
                "# Dependency closure scan completed. No non-system DLLs were resolved via dumpbin.`n# Scanned $scanCount unique modules in $([math]::Round($stopwatch.Elapsed.TotalSeconds, 2))s.`n# Executables scanned: $($ExePaths.Count)" | Out-File -FilePath $closureFile -Encoding UTF8 -Force
                if ($scanCount -eq 0) {
                    Write-Warning "Dependency closure scanned 0 modules. Check if executables were found and dumpbin is working."
                } else {
                    Write-Warning "Dependency closure found no DLLs to resolve (scanned $scanCount modules). This may indicate all dependencies are system DLLs."
                }
            }
        } catch {
            Write-Error "Failed to write deps-closure.txt: $_"
            throw
        }
    }

    # Print summary
    $elapsedSeconds = [math]::Round($stopwatch.Elapsed.TotalSeconds, 2)
    Write-Host "Dependency closure complete: $($copied.Count) DLLs copied, $($resolved.Count) DLLs resolved, $scanCount unique modules scanned in ${elapsedSeconds}s"
    if ($missing.Count -gt 0) {
        $missingDlls = $missing | ForEach-Object { if ($_ -is [hashtable]) { $_.DllName } else { $_ } } | Sort-Object -Unique
        $missingSample = $missingDlls | Select-Object -First 10
        # Note: Most unresolved DLLs are Windows system DLLs (e.g., api-ms-win-*, ext-ms-win-*, kernel32.dll, etc.)
        # that are provided by the OS and not redistributable. This is expected and not an error.
        Write-Warning "Unresolved non-system DLLs ($($missingDlls.Count) total, mostly Windows system DLLs): $($missingSample -join ', ')"
        if ($missingDlls.Count -gt 10) {
            Write-Warning "  ... and $($missingDlls.Count - 10) more unresolved DLLs (mostly system DLLs)"
        }
    }
    
    # Return summary info
    return @{
        ScannedCount = $scanCount
        CopiedCount = $copied.Count
        ResolvedCount = $resolved.Count
        UnresolvedCount = $missing.Count
        UnresolvedDlls = $missing
        ElapsedSeconds = $elapsedSeconds
    }
}

function Copy-MustHaveDllPatterns {
    param(
        [string[]]$SearchRoots,
        [string]$DistDir,
        [string]$OneApiRoot,
        [string]$MklRoot
    )
    Write-Host "Copying must-have DLL patterns (runtime safety net)..."

    # Build extended search roots from environment variables
    $extendedRoots = $SearchRoots | ForEach-Object { $_ }
    if ($env:ONEAPI_ROOT -and (Test-Path $env:ONEAPI_ROOT)) {
        $extendedRoots += Get-SearchRoots -OneApiRoot $env:ONEAPI_ROOT
    }
    if ($env:MKLROOT -and (Test-Path $env:MKLROOT)) {
        $mklBin = Join-Path $env:MKLROOT "bin"
        $mklRedist = Join-Path $env:MKLROOT "redist\intel64"
        if (Test-Path $mklBin) { $extendedRoots += $mklBin }
        if (Test-Path $mklRedist) { $extendedRoots += $mklRedist }
    }
    if ($env:I_MPI_ROOT -and (Test-Path $env:I_MPI_ROOT)) {
        $mpiBin = Join-Path $env:I_MPI_ROOT "bin"
        if (Test-Path $mpiBin) { $extendedRoots += $mpiBin }
    }
    $extendedRoots = $extendedRoots | Select-Object -Unique | Where-Object { $_ -and (Test-Path $_) }

    # Define explicit allowlist: only these DLLs should be staged
    # This policy is based on multiple validated builds (local + CI artifacts)
    $categories = @{
        "Intel Fortran Runtime" = @{
            Patterns = @("libifcoremd*.dll", "libifportmd*.dll", "libmmd*.dll")
            Required = $true
            Found = 0
            Description = "Required for any ifx-compiled Fortran code"
        }
        "Intel SVML" = @{
            Patterns = @("svml_dispmd*.dll")
            Required = $true
            Found = 0
            Description = "Vector math functions used by ifx-compiled code"
        }
        "Intel OpenMP Runtime" = @{
            Patterns = @("libiomp5md*.dll")
            Required = $true
            Found = 0
            PreferRelease = $true  # Prefer release DLL, avoid _db.dll
            Description = "Required by QE and MKL threading"
        }
        "MKL Core" = @{
            Patterns = @("mkl_core*.dll", "mkl_def*.dll", "mkl_rt*.dll", "mkl_intel_thread*.dll")
            Required = $true
            Found = 0
            Description = "BLAS / LAPACK / FFT backend"
        }
        "MKL CPU Features" = @{
            Patterns = @("mkl_avx*.dll", "mkl_vml*.dll", "mkl_mc*.dll")
            Required = $false
            Found = 0
            Comment = "Optional: Enable runtime CPU dispatch (AVX2, etc.) for better performance. QE will run without them."
        }
    }
    
    # Define explicit denylist: these must NEVER be staged
    $denylistPatterns = @(
        "*_db.dll",      # Debug runtimes (e.g., libiomp5md_db.dll)
        "libimalloc.dll" # Intel scalable allocator (not required)
    )

    foreach ($catName in $categories.Keys) {
        $cat = $categories[$catName]
        Write-Host "  Category: $catName"
        if ($cat.Comment) {
            Write-Host "    ($($cat.Comment))"
        }
        
        foreach ($pattern in $cat.Patterns) {
            $found = $null
            $allMatches = @()
            
            # Collect all matches first (for OpenMP release preference logic)
            foreach ($root in $extendedRoots) {
                $matches = Get-ChildItem -Path $root -Filter $pattern -File -Recurse -ErrorAction SilentlyContinue
                if ($matches) {
                    $allMatches += $matches
                }
            }
            
            # Filter all matches against denylist first
            $filteredMatches = $allMatches | Where-Object { 
                $name = $_.Name
                $denied = $false
                foreach ($denyPattern in $denylistPatterns) {
                    if ($name -like $denyPattern) {
                        $denied = $true
                        break
                    }
                }
                -not $denied
            }
            
            # For OpenMP with PreferRelease, prioritize exact libiomp5md.dll
            if ($cat.PreferRelease -and $pattern -eq "libiomp5md*.dll") {
                # First try exact match: libiomp5md.dll (already filtered by denylist)
                $exactMatch = $filteredMatches | Where-Object { $_.Name -eq "libiomp5md.dll" } | Select-Object -First 1
                if ($exactMatch) {
                    $found = $exactMatch
                    Write-Host "    [OpenMP] Selected exact match: libiomp5md.dll (release)"
                } else {
                    # Fallback: any non-_db libiomp5md*.dll (denylist already filtered _db.dll)
                    $found = $filteredMatches | Select-Object -First 1
                    if ($found) {
                        Write-Host "    [OpenMP] Selected: $($found.Name) (release)"
                    }
                }
            } else {
                # Default behavior: first match from filtered list
                $found = $filteredMatches | Select-Object -First 1
            }
            
            if ($found) {
                $dest = Join-Path $DistDir (Split-Path $found.FullName -Leaf)
                $alreadyExists = (Test-Path $dest)
                if (-not $alreadyExists) {
                    if ($PSCmdlet.ShouldProcess($dest, "Copy pattern $pattern") -and -not $WhatIfPreference) {
                        Copy-Item -Force $found.FullName $dest
                    }
                    Write-Host "    + $pattern -> $($found.Name) (from $($found.DirectoryName))"
                    $cat.Found++
                } else {
                    Write-Host "    = $pattern (already staged)"
                    $cat.Found++
                }
            } else {
                Write-Host "    - $pattern (not found)"
            }
        }
        
        if ($cat.Required -and $cat.Found -eq 0) {
            throw "Required category '$catName' found no DLLs. Staging failed. Check oneAPI installation and environment variables (ONEAPI_ROOT, MKLROOT)."
        } elseif (-not $cat.Required -and $cat.Found -eq 0) {
            Write-Warning "Optional category '$catName' found no DLLs. Consider installing VC++ Redistributable for MSVC runtime."
        }
    }

    Write-Host "Must-have patterns complete. Summary:"
    foreach ($catName in $categories.Keys) {
        $cat = $categories[$catName]
        $status = if ($cat.Required) { "REQUIRED" } else { "optional" }
        Write-Host "  $catName ($status): $($cat.Found) DLL(s) found"
    }
    
    # Post-staging validation: enforce denylist
    Write-Host ""
    Write-Host "Validating staged DLLs against denylist..." -ForegroundColor Cyan
    $violations = @()
    Get-ChildItem -Path $DistDir -Filter "*.dll" -File | ForEach-Object {
        if (Is-DenylistedDll -Name $_.Name) {
            $violations += $_.Name
        }
    }
    
    if ($violations.Count -gt 0) {
        Write-Error "DENYLIST VIOLATION: The following DLLs must not be staged:"
        foreach ($violation in $violations) {
            Write-Error "  - $violation"
        }
        Write-Error "Staging failed. These DLLs should not be included in the distribution."
        throw "Denylist validation failed. Found $($violations.Count) prohibited DLL(s)."
    } else {
        Write-Host "  Validation passed: No denylisted DLLs found" -ForegroundColor Green
    }
}

function Stage-MSMPIRuntime {
    param(
        [string]$DistDir
    )
    
    Write-Host "=== Staging MS-MPI Runtime ===" -ForegroundColor Cyan
    
    # Primary location: MS-MPI runtime install directory
    $msmpiBin = "C:\Program Files\Microsoft MPI\Bin"
    
    # Files to stage
    $filesToStage = @(
        @{ Name = "msmpi.dll"; Required = $true; Source = $null }
        @{ Name = "mpiexec.exe"; Required = $true; Source = $null }
        @{ Name = "smpd.exe"; Required = $true; Source = $null }
    )
    
    # Find files in primary location
    foreach ($file in $filesToStage) {
        $candidate = Join-Path $msmpiBin $file.Name
        if (Test-Path $candidate) {
            $file.Source = $candidate
        }
    }
    
    # Fallback for msmpi.dll: check System32
    if (-not $filesToStage[0].Source) {
        $fallback = Join-Path $env:WINDIR "System32\msmpi.dll"
        if (Test-Path $fallback) {
            $filesToStage[0].Source = $fallback
        }
    }
    
    # Stage each file
    foreach ($file in $filesToStage) {
        if ($file.Source -and (Test-Path $file.Source)) {
            $dest = Join-Path $DistDir $file.Name
            if ($PSCmdlet.ShouldProcess($dest, "Copy MS-MPI $($file.Name)") -and -not $WhatIfPreference) {
                Copy-Item -Force $file.Source $dest
            }
            Write-Host "  + $($file.Name) (from $($file.Source))" -ForegroundColor Green
        } else {
            if ($file.Required) {
                throw "Required MS-MPI file '$($file.Name)' not found. Expected at: $msmpiBin\$($file.Name) or (for msmpi.dll only) $env:WINDIR\System32\msmpi.dll"
            } else {
                Write-Host "  - $($file.Name) (optional, not found)" -ForegroundColor Yellow
            }
        }
    }
    
    Write-Host "MS-MPI runtime staging complete." -ForegroundColor Green
}


function Run-SmokeTest {
    param(
        [string]$DistDir,
        [string]$ResourcesDir,
        [switch]$Strict,  # Strict mode: fail if exit code != 0 OR "JOB DONE" not detected
        [switch]$UseMpi,   # If set, run under mpiexec
        [string]$MpiLauncher,  # Path to mpiexec.exe
        [int]$MpiNP = 2   # Number of MPI processes
    )
    
    if (-not $DistDir) {
        throw "DistDir parameter is required for smoke test"
    }
    
    if (-not $ResourcesDir) {
        throw "ResourcesDir parameter is required for smoke test"
    }
    
    Write-Host "Running clean-room smoke test (15s timeout)..."
    
    if (-not (Test-Path -Path $ResourcesDir)) {
        throw "Shared resources folder not found at: $ResourcesDir"
    }
    
    $smokeDir = $DistDir  # Run test in dist/bin (where pw.exe is)
    $smokeInput = Join-Path $smokeDir "scf-cg.in"
    $smokeOutput = Join-Path $smokeDir "scf-cg.out"
    $pseudoFile = Join-Path $smokeDir "Si.pz-vbc.UPF"
    $binDir = $DistDir  # Executables are in dist/bin
    $pwExe = Join-Path $binDir "pw.exe"
    if (-not (Test-Path $pwExe)) {
        $pwExe = Join-Path $binDir "pw.x.exe"
    }
    if (-not (Test-Path $pwExe)) {
        Write-Warning "pw.exe not found for smoke test, skipping"
        return
    }

    # Copy shared resources to smoke directory
    $resourceInput = Join-Path $resourcesDir "scf-cg.in"
    $resourcePseudo = Join-Path $resourcesDir "Si.pz-vbc.UPF"
    
    if (-not (Test-Path $resourceInput)) {
        throw "Shared input file not found: $resourceInput"
    }
    if (-not (Test-Path $resourcePseudo)) {
        throw "Shared pseudopotential file not found: $resourcePseudo"
    }
    
    if (-not $WhatIfPreference) {
        Copy-Item -Force $resourceInput $smokeInput
        Copy-Item -Force $resourcePseudo $pseudoFile
        Write-Host "  Copied test resources: scf-cg.in, Si.pz-vbc.UPF"
    }

    # Run smoke test with clean-room PATH (only dist/bin prepended)
    # Use shell redirection to match manual execution: pw.exe -i scf-cg.in > scf-cg.out
    $cleanPath = "$binDir;$env:PATH"
    $SMOKE_TIMEOUT_MS = 15000  # 15 seconds (calculation takes ~1s, this is abundant)
    
    # Build command line
    if ($UseMpi) {
        if (-not $MpiLauncher -or -not (Test-Path $MpiLauncher)) {
            throw "MPI launcher not found at: $MpiLauncher. Ensure mpiexec.exe is staged in dist root."
        }
        Write-Host "  Running under MPI: mpiexec -n $MpiNP"
        # Use relative paths with pushd/popd in a single cmd.exe /c string
        $cmdLine = "pushd `"$DistDir`" && `".\mpiexec.exe`" -n $MpiNP `".\pw.exe`" -in `"scf-cg.in`" > `"scf-cg.out`" 2>&1 && popd"
        Write-Host "  MPI command: $cmdLine" -ForegroundColor Cyan
    } else {
        $cmdLine = "`"$pwExe`" -i scf-cg.in > scf-cg.out 2>&1"
    }
    
    if (-not $WhatIfPreference) {
        # Use cmd.exe to run with shell redirection (avoids deadlock issues with large output)
        # This matches the manual command: pw.exe -i scf-cg.in > scf-cg.out
        $psi = New-Object System.Diagnostics.ProcessStartInfo
        $psi.FileName = "cmd.exe"
        $psi.Arguments = "/c $cmdLine"
        $psi.UseShellExecute = $false
        $psi.WorkingDirectory = $smokeDir
        $psi.Environment["PATH"] = $cleanPath
        $psi.RedirectStandardOutput = $true
        $psi.RedirectStandardError = $true
        $proc = New-Object System.Diagnostics.Process
        $proc.StartInfo = $psi
        
        $proc.Start() | Out-Null
        if (-not $proc.WaitForExit($SMOKE_TIMEOUT_MS)) {
            Write-Warning "Smoke test timed out after ${SMOKE_TIMEOUT_MS}ms"
            $proc.Kill()
            $proc.WaitForExit() | Out-Null
            $exitCode = -1
        } else {
            $exitCode = $proc.ExitCode
        }
        
        # Read the output file that was created by shell redirection
        if (Test-Path -Path $smokeOutput) {
            $combined = Get-Content -Path $smokeOutput -Raw -ErrorAction SilentlyContinue
            if (-not $combined) {
                $combined = ""
            }
            # Show last 5 lines of output for debugging
            $outputLines = Get-Content -Path $smokeOutput -ErrorAction SilentlyContinue
            if ($outputLines) {
                $lastLines = $outputLines | Select-Object -Last 5
                Write-Host "  Last 5 lines of ${smokeOutput}:"
                foreach ($line in $lastLines) {
                    Write-Host "    $line"
                }
            }
        } else {
            $combined = ""
            Write-Warning "  Output file not found: $smokeOutput"
        }
        
        # Check for "JOB DONE" in the .out file (case-insensitive)
        $jobDone = ($combined -match "JOB\s+DONE" -or $combined -match "job\s+done")
        
        # Lightweight MPI process count verification (informational only, does not affect pass/fail)
        if ($UseMpi -and $combined) {
            $foundNP = $null
            $mpiVersionFound = $false
            
            # Check for "Parallel version (MPI), running on" (OpenMP disabled)
            if ($combined -match "Parallel version \(MPI\), running on") {
                $mpiVersionFound = $true
            }
            # Check for "Parallel version (MPI & OpenMP), running on" (OpenMP enabled)
            elseif ($combined -match "Parallel version \(MPI & OpenMP\), running on") {
                $mpiVersionFound = $true
            }
            
            if ($mpiVersionFound) {
                # Try to extract MPI process count from "Number of MPI processes:" line
                if ($combined -match "Number of MPI processes:\s+(\d+)") {
                    $foundNP = $Matches[1]
                    if ($foundNP -eq $MpiNP) {
                        Write-Host "  MPI process count verified: found '$MpiNP' processes" -ForegroundColor Green
                    } else {
                        Write-Warning "  MPI process count mismatch: expected '$MpiNP' but found '$foundNP'"
                    }
                } else {
                    Write-Warning "  Parallel version detected but could not extract MPI process count from 'Number of MPI processes:' line"
                }
            } else {
                Write-Warning "  MPI process count line not found in output (expected 'Parallel version (MPI)' or 'Parallel version (MPI & OpenMP)')"
            }
        }
        
        # Format exit code
        $exitCodeHex = if ($exitCode -ge 0) { "0x{0:X8}" -f $exitCode } else { "N/A" }
        
        Write-Host "Smoke test results:"
        Write-Host "  ExitCode: $exitCode ($exitCodeHex)"
        Write-Host "  JobDone: $jobDone"
        Write-Host "  StrictMode: $Strict"
        Write-Host "  Output: $smokeOutput"
        
        # Check for loader errors
        if ($exitCode -eq 0xC0000135 -or $exitCode -eq 0xC000007B) {
            Write-Warning "DLL missing or incompatible image (exit code indicates loader failure)"
        }
        
        # Acceptance criteria (strict mode is default)
        $passed = $false
        if ($Strict) {
            if ($exitCode -eq 0 -and $jobDone) {
                $passed = $true
            } else {
                Write-Error "Smoke test FAILED (ExitCode=$exitCode, JobDone=$jobDone)"
                throw "Smoke test failed: exit code must be 0 and 'JOB DONE' must be detected in output"
            }
        } else {
            # Non-strict mode (not recommended, only for debugging)
            if ($exitCode -eq 0) {
                $passed = $true
                if (-not $jobDone) {
                    Write-Warning "Smoke test passed (exit 0) but 'JOB DONE' not detected. Check $smokeOutput"
                }
            } else {
                Write-Warning "Smoke test failed (exit $exitCode). Check $smokeOutput"
            }
        }
        
        if ($passed) {
            Write-Host "Smoke test PASSED" -ForegroundColor Green
        }
    } else {
        if ($UseMpi) {
            Write-Host "  (WhatIf: would run pushd `"$DistDir`" && .\mpiexec.exe -n $MpiNP .\pw.exe -in scf-cg.in > scf-cg.out 2>&1 && popd)"
        } else {
            Write-Host "  (WhatIf: would run $pwExe -i scf-cg.in > scf-cg.out with clean PATH)"
        }
    }
}

function Get-StackReserveLine {
    param(
        [string]$DumpbinPath,
        [string]$PwExe
    )
    
    # Use PowerShell's native filtering to avoid stdout deadlock
    # Run dumpbin and pipe directly to Select-String for filtering
    # This avoids cmd.exe escaping issues and works reliably
    try {
        # Run dumpbin with stderr redirected, pipe to Select-String for filtering
        $output = & $DumpbinPath /headers $PwExe 2>&1 | Select-String -Pattern "size of stack reserve" -CaseSensitive:$false
        
        if ($output) {
            # Select-String returns MatchInfo objects
            if ($output -is [array]) {
                # Multiple matches - take the first one (should be only one anyway)
                $line = $output[0].Line
            } else {
                $line = $output.Line
            }
            return $line.Trim()
        }
        return $null
    } catch {
        return $null
    }
}

function Get-QeVersion {
    param(
        [string]$QeSourceDir
    )
    # Try to get version from git
    $gitDir = Join-Path $QeSourceDir ".git"
    if (Test-Path $gitDir) {
        try {
            Push-Location $QeSourceDir
            $gitSha = git rev-parse --short HEAD 2>$null
            $gitTag = git describe --tags --exact-match 2>$null
            if ($gitTag) {
                return $gitTag
            } elseif ($gitSha) {
                return "git-$gitSha"
            }
        } catch {
            # Ignore git errors
        } finally {
            Pop-Location
        }
    }
    # Try to read from version file or CMakeLists.txt
    $versionFile = Join-Path $QeSourceDir "VERSION"
    if (Test-Path $versionFile) {
        $version = Get-Content $versionFile -First 1 -ErrorAction SilentlyContinue
        if ($version) { return $version.Trim() }
    }
    return "unknown"
}

function Get-OneApiVersion {
    param(
        [string]$OneApiRoot,
        [string]$DistBinDir  # Optional: check bundled DLLs if licensing path fails
    )
    # Strategy 1: Try to extract version from licensing directory path
    # e.g., "licensing\2025.3\licensing\2025.3\license.htm" -> "2025.3"
    $licensingBase = Join-Path $OneApiRoot "licensing"
    if (Test-Path $licensingBase) {
        # Check latest\licensing\* first
        $latestLicensePath = Join-Path $licensingBase "latest\licensing"
        if (Test-Path $latestLicensePath) {
            $versionDirs = Get-ChildItem -Path $latestLicensePath -Directory -ErrorAction SilentlyContinue | 
                Where-Object { $_.Name -match '^\d+\.\d+' }
            if ($versionDirs) {
                $version = ($versionDirs | Sort-Object { [Version]$_.Name } -Descending)[0].Name
                Write-Host "  Extracted oneAPI version from licensing path: $version" -ForegroundColor Gray
                return $version
            }
        }
        
        # Check versioned directories (e.g., licensing\2025.3\licensing\2025.3\)
        $versionDirs = Get-ChildItem -Path $licensingBase -Directory -ErrorAction SilentlyContinue | 
            Where-Object { $_.Name -match '^\d+\.\d+' } | 
            Sort-Object { [Version]$_.Name } -Descending
        
        foreach ($versionDir in $versionDirs) {
            $versionLicensePath = Join-Path $versionDir.FullName "licensing"
            if (Test-Path $versionLicensePath) {
                $subDirs = Get-ChildItem -Path $versionLicensePath -Directory -ErrorAction SilentlyContinue | 
                    Where-Object { $_.Name -match '^\d+\.\d+' }
                if ($subDirs) {
                    $version = ($subDirs | Sort-Object { [Version]$_.Name } -Descending)[0].Name
                    Write-Host "  Extracted oneAPI version from licensing path: $version" -ForegroundColor Gray
                    return $version
                }
                # Fallback: use parent directory version
                $version = $versionDir.Name
                Write-Host "  Extracted oneAPI version from licensing directory: $version" -ForegroundColor Gray
                return $version
            }
        }
    }
    
    # Strategy 2: Try to get version from bundled Intel DLL (if dist/bin provided)
    if ($DistBinDir -and (Test-Path $DistBinDir)) {
        # Try libifcoremd.dll first (Intel Fortran runtime, always present)
        $intelDlls = @("libifcoremd.dll", "mkl_rt*.dll", "libiomp5md.dll")
        foreach ($dllPattern in $intelDlls) {
            $dllFiles = Get-ChildItem -Path $DistBinDir -Filter $dllPattern -File -ErrorAction SilentlyContinue | Select-Object -First 1
            if ($dllFiles) {
                try {
                    $v = [System.Diagnostics.FileVersionInfo]::GetVersionInfo($dllFiles.FullName)
                    $version = if ($v.ProductVersion) { $v.ProductVersion } else { $v.FileVersion }
                    if ($version -and $version -ne "0.0.0.0") {
                        Write-Host "  Extracted oneAPI version from bundled DLL ($($dllFiles.Name)): $version (runtime file version)" -ForegroundColor Gray
                        return $version
                    }
                } catch {
                    # Ignore errors
                }
            }
        }
    }
    
    # Strategy 3: Try environment variable
    if ($env:ONEAPI_VERSION) {
        return $env:ONEAPI_VERSION
    }
    
    Write-Warning "Could not determine Intel oneAPI version"
    return "unknown"
}

function Get-MsmpiVersion {
    param(
        [string]$DistBinDir  # Check bundled msmpi.dll first
    )
    # Strategy 1: Check bundled msmpi.dll in dist/bin (most reliable)
    if ($DistBinDir -and (Test-Path $DistBinDir)) {
        $msmpiDll = Join-Path $DistBinDir "msmpi.dll"
        if (Test-Path $msmpiDll) {
            try {
                $v = [System.Diagnostics.FileVersionInfo]::GetVersionInfo($msmpiDll)
                $version = if ($v.ProductVersion) { $v.ProductVersion } else { $v.FileVersion }
                if ($version -and $version -ne "0.0.0.0") {
                    Write-Host "  Extracted MS-MPI version from bundled msmpi.dll: $version" -ForegroundColor Gray
                    return $version
                }
            } catch {
                # Ignore errors, try fallback
            }
        }
    }
    
    # Strategy 2: Try installed MS-MPI
    $msmpiBin = "C:\Program Files\Microsoft MPI\Bin"
    if (Test-Path $msmpiBin) {
        $msmpiDll = Join-Path $msmpiBin "msmpi.dll"
        if (Test-Path $msmpiDll) {
            try {
                $v = [System.Diagnostics.FileVersionInfo]::GetVersionInfo($msmpiDll)
                $version = if ($v.ProductVersion) { $v.ProductVersion } else { $v.FileVersion }
                if ($version -and $version -ne "0.0.0.0") {
                    Write-Host "  Extracted MS-MPI version from installed msmpi.dll: $version" -ForegroundColor Gray
                    return $version
                }
            } catch {
                # Ignore errors
            }
        }
    }
    
    Write-Warning "Could not determine MS-MPI version"
    return "unknown"
}

function Copy-QeLicenses {
    param(
        [string]$QeSourceDir,
        [string]$LicensesDir
    )
    $qeLicenseDir = Join-Path $LicensesDir "quantum-espresso"
    if (-not (Test-Path $qeLicenseDir)) {
        New-Item -ItemType Directory -Force -Path $qeLicenseDir | Out-Null
    }
    
    $copied = @()
    $missing = @()
    
    # Main QE license (required)
    $mainLicense = Join-Path $QeSourceDir "License"
    if (Test-Path $mainLicense) {
        $dst = Join-Path $qeLicenseDir "License"
        Copy-Item -Force $mainLicense $dst
        $copied += "License"
        Write-Host "  Copied QE main license: License" -ForegroundColor Green
    } else {
        throw "Required QE main license file not found: $mainLicense"
    }
    
    # External library licenses (including EPW for consistency)
    $externalLicenses = @(
        @{ Path = "EPW\License"; Name = "external/EPW/License" },
        @{ Path = "external\wannier90\LICENSE"; Name = "external/wannier90/LICENSE" },
        @{ Path = "external\lapack\LICENSE"; Name = "external/lapack/LICENSE" },
        @{ Path = "external\fox\LICENSE"; Name = "external/fox/LICENSE" },
        @{ Path = "external\mbd\LICENSE"; Name = "external/mbd/LICENSE" },
        @{ Path = "external\qe-gipaw\LICENSE"; Name = "external/qe-gipaw/LICENSE" },
        @{ Path = "external\d3q\License"; Name = "external/d3q/License" }
    )
    
    foreach ($lic in $externalLicenses) {
        $src = Join-Path $QeSourceDir $lic.Path
        if (Test-Path $src) {
            # Create directory structure using the target name path (e.g., external/EPW/License -> external/EPW/)
            # This ensures EPW goes to external/EPW/ not just EPW/
            $targetRelativePath = Split-Path $lic.Name -Parent
            $targetDir = Join-Path $qeLicenseDir $targetRelativePath
            if (-not (Test-Path $targetDir)) {
                New-Item -ItemType Directory -Force -Path $targetDir | Out-Null
            }
            $fileName = Split-Path $lic.Path -Leaf
            $dst = Join-Path $targetDir $fileName
            Copy-Item -Force $src $dst
            $copied += $lic.Name
            Write-Host "  Copied QE license: $($lic.Name)" -ForegroundColor Green
        } else {
            Write-Warning "  QE license not found: $($lic.Path)"
        }
    }
    
    return @{ Copied = $copied; Missing = $missing }
}

function Copy-OneApiLicenses {
    param(
        [string]$OneApiRoot,
        [string]$LicensesDir,
        [string]$RepoRootPath
    )
    $oneApiLicenseDir = Join-Path $LicensesDir "intel-oneapi"
    if (-not (Test-Path $oneApiLicenseDir)) {
        New-Item -ItemType Directory -Force -Path $oneApiLicenseDir | Out-Null
    }
    
    $copied = @()
    
    # Strategy 1: Try latest\licensing\*\license.htm (robust for CI)
    $latestLicensePath = Join-Path $OneApiRoot "licensing\latest\licensing"
    if (Test-Path $latestLicensePath) {
        $licenseFiles = Get-ChildItem -Path $latestLicensePath -Filter "license.htm" -Recurse -ErrorAction SilentlyContinue
        foreach ($file in $licenseFiles) {
            $dst = Join-Path $oneApiLicenseDir "license.htm"
            Copy-Item -Force $file.FullName $dst
            $copied += "license.htm"
            Write-Host "  Copied oneAPI license: license.htm (from latest)" -ForegroundColor Green
            break  # Take first match
        }
    }
    
    # Strategy 2.1: Try versioned path (e.g., 2025.3\licensing\2025.3\license.htm)
    if ($copied.Count -eq 0) {
        $licensingBase = $OneApiRoot
        if (Test-Path $licensingBase) {
            # Find versioned directories
            $versionDirs = Get-ChildItem -Path $licensingBase -Directory -ErrorAction SilentlyContinue | 
                Where-Object { $_.Name -match '^\d+\.\d+' } | 
                Sort-Object { [Version]$_.Name } -Descending
            
            foreach ($versionDir in $versionDirs) {
                $versionLicensePath = Join-Path $versionDir.FullName "licensing"
                if (Test-Path $versionLicensePath) {
                    $versionSubDirs = Get-ChildItem -Path $versionLicensePath -Directory -ErrorAction SilentlyContinue
                    foreach ($subDir in $versionSubDirs) {
                        $licenseFile = Join-Path $subDir.FullName "license.htm"
                        if (Test-Path $licenseFile) {
                            $dst = Join-Path $oneApiLicenseDir "license.htm"
                            Copy-Item -Force $licenseFile $dst
                            $copied += "license.htm"
                            Write-Host "  Copied oneAPI license: license.htm (from $($versionDir.Name))" -ForegroundColor Green
                            break
                        }
                    }
                    if ($copied.Count -gt 0) { break }
                }
            }
        }
    }
    
    # Strategy 2.2: Try versioned path (e.g., licensing\2025.3\licensing\2025.3\license.htm)
    if ($copied.Count -eq 0) {
        $licensingBase = Join-Path $OneApiRoot "licensing"
        if (Test-Path $licensingBase) {
            # Find versioned directories
            $versionDirs = Get-ChildItem -Path $licensingBase -Directory -ErrorAction SilentlyContinue | 
                Where-Object { $_.Name -match '^\d+\.\d+' } | 
                Sort-Object { [Version]$_.Name } -Descending
            
            foreach ($versionDir in $versionDirs) {
                $versionLicensePath = Join-Path $versionDir.FullName "licensing"
                if (Test-Path $versionLicensePath) {
                    $versionSubDirs = Get-ChildItem -Path $versionLicensePath -Directory -ErrorAction SilentlyContinue
                    foreach ($subDir in $versionSubDirs) {
                        $licenseFile = Join-Path $subDir.FullName "license.htm"
                        if (Test-Path $licenseFile) {
                            $dst = Join-Path $oneApiLicenseDir "license.htm"
                            Copy-Item -Force $licenseFile $dst
                            $copied += "license.htm"
                            Write-Host "  Copied oneAPI license: license.htm (from $($versionDir.Name))" -ForegroundColor Green
                            break
                        }
                    }
                    if ($copied.Count -gt 0) { break }
                }
            }
        }
    }

    # Strategy 3: Fallback - search for any license.htm in licensing directory
    if ($copied.Count -eq 0) {
        $licensingBase = Join-Path $OneApiRoot "licensing"
        if (Test-Path $licensingBase) {
            $licenseFiles = Get-ChildItem -Path $licensingBase -Filter "license.htm" -Recurse -ErrorAction SilentlyContinue | Select-Object -First 1
            if ($licenseFiles) {
                $dst = Join-Path $oneApiLicenseDir "license.htm"
                Copy-Item -Force $licenseFiles.FullName $dst
                $copied += "license.htm"
                Write-Host "  Copied oneAPI license: license.htm (fallback search)" -ForegroundColor Green
            }
        }
    }
    if ($copied.Count -eq 0) {
        $LicenseFallbackDir = Join-Path $RepoRootPath "license_fallbacks\intel-oneapi"
        if (Test-Path $LicenseFallbackDir) {
            $licenseFiles = Get-ChildItem -Path $LicenseFallbackDir -Filter "license.htm" -Recurse -ErrorAction SilentlyContinue | Select-Object -First 1
            if ($licenseFiles) {
                $dst = Join-Path $oneApiLicenseDir "license.htm"
                Copy-Item -Force $licenseFiles.FullName $dst
                $copied += "license.htm"
                Write-Host "  Copied oneAPI license: license.htm (from repo fallback)" -ForegroundColor Yellow
            }
        }
    }
    if ($copied.Count -eq 0) {
        Write-Warning "No Intel oneAPI license file (license.htm) found under $OneApiRoot\licensing"
    }
    
    return @{ Copied = $copied }
}

function Copy-MsmpiLicense {
    param(
        [string]$LicensesDir
    )
    $msmpiLicenseDir = Join-Path $LicensesDir "microsoft-msmpi"
    if (-not (Test-Path $msmpiLicenseDir)) {
        New-Item -ItemType Directory -Force -Path $msmpiLicenseDir | Out-Null
    }
    
    # Only use Runtime License, not SDK
    # Note: License is a directory, not a file - need to copy files from inside it
    $msmpiRuntimeLicenseDir = Join-Path $env:ProgramFiles "Microsoft MPI\License"
    if (-not (Test-Path $msmpiRuntimeLicenseDir)) {
        # Try alternative path
        $msmpiRuntimeLicenseDir = "C:\Program Files\Microsoft MPI\License"
    }
    
    $copied = @()
    if (Test-Path $msmpiRuntimeLicenseDir) {
        # Copy all license files from the License directory
        $licenseFiles = Get-ChildItem -Path $msmpiRuntimeLicenseDir -File -ErrorAction SilentlyContinue
        foreach ($file in $licenseFiles) {
            $dst = Join-Path $msmpiLicenseDir $file.Name
            Copy-Item -Force $file.FullName $dst
            $copied += $file.Name
            Write-Host "  Copied MS-MPI Runtime license: $($file.Name)" -ForegroundColor Green
        }
        
        if ($copied.Count -gt 0) {
            return @{ Copied = $copied }
        } else {
            Write-Warning "MS-MPI Runtime license directory found but contains no files: $msmpiRuntimeLicenseDir"
        }
    } else {
        Write-Warning "MS-MPI Runtime license directory not found at: $msmpiRuntimeLicenseDir"
    }
    
    return @{ Copied = $copied }
}


function Generate-ThirdPartyNotices {
    param(
        [string]$LicensesDir,
        [string]$QeVersion,
        [string]$OneApiVersion,
        [string]$MsmpiVersion,
        [string]$BuildDateUtc
    )
    $noticesFile = Join-Path $LicensesDir "THIRD_PARTY_NOTICES.txt"
    $content = @"
QMatSuite Toolchain - Quantum ESPRESSO Windows (oneAPI + MS-MPI) Binary Distribution

Build date (UTC): $BuildDateUtc

This package redistributes precompiled binaries of Quantum ESPRESSO and bundles
third-party runtime components needed to run on Windows without requiring a
separate oneAPI/MS-MPI installation. Each component remains licensed under its
own license terms. See the subfolders in this directory for the full license
texts.

Included components

1) Quantum ESPRESSO ($QeVersion)
   Copyright: Quantum ESPRESSO Foundation and contributors
   License: GNU General Public License v2 or later (GPL-2.0-or-later)
   License text: licenses/quantum-espresso/

2) Intel oneAPI Runtime Libraries ($OneApiVersion)
   Examples: Intel OpenMP (libiomp5md), Intel MKL (mkl_rt), compiler runtimes
   License: Intel terms for redistributable runtime components
   License text: licenses/intel-oneapi/

3) Microsoft MPI (MS-MPI) Runtime ($MsmpiVersion)
   License: Microsoft Software License Terms
   License text: licenses/microsoft-msmpi/

4) External libraries linked by Quantum ESPRESSO
   License: Check individual library licenses
   License text: licenses/quantum-espresso/external/

Notes

- Windows system components (e.g., api-ms-win-* and other OS-provided DLLs) are not
  bundled.

- If any license text could not be found during staging, the build log will
  include a prominent warning. Please refer to the vendor documentation for the
  corresponding version.
"@
    
    $content | Out-File -FilePath $noticesFile -Encoding UTF8 -NoNewline
    Write-Host "  Generated THIRD_PARTY_NOTICES.txt" -ForegroundColor Green
}

function Generate-Readme {
    param(
        [string]$DistRoot,
        [string]$QeVersion,
        [string]$BuildDateUtc,
        [int]$ExeCount,
        [int]$DllCount,
        [bool]$HasMsmpiLauncher
    )

    $readmeFile = Join-Path $DistRoot "README.txt"

    $distLayoutSummary = "Contains $ExeCount executable(s) in bin/ and $DllCount runtime DLL(s) in bin/."

    # --- code snippets (single-quoted here-strings => NO variable expansion) ---
    $pathHint = @'
Optional: add bin/ to PATH for this terminal session only:

  set "PATH=%CD%;%PATH%"               (CMD, when you are already inside bin\)
  $env:PATH = "$PWD;$env:PATH"         (PowerShell, when you are already inside bin\)

If you are NOT inside bin\ yet, use:

  set "PATH=%CD%\bin;%PATH%"           (CMD)
  $env:PATH = "$PWD\bin;$env:PATH"     (PowerShell)
'@

    $ompHint = @'
Set OMP_NUM_THREADS to control OpenMP threading:

  set OMP_NUM_THREADS=8                (CMD)
  $env:OMP_NUM_THREADS=8               (PowerShell)
'@

    $mklHint = @'
Set MKL_NUM_THREADS to control Intel MKL threading (recommended when using OpenMP):

  set MKL_NUM_THREADS=8                (CMD)
  $env:MKL_NUM_THREADS=8               (PowerShell)
'@

    $msmpiHint = if ($HasMsmpiLauncher) {
@'
MS-MPI launcher is included (mpiexec/msmpiexec) in bin/.

Tip:
- If you run from inside bin\, mpiexec should be found automatically.
- Otherwise, add bin/ to PATH for the current terminal session (see the PATH hint above).
'@
    } else {
@'
MS-MPI launcher is NOT bundled.

To run MPI jobs, install Microsoft MPI and ensure mpiexec is available in PATH.
(You can still run serial/OpenMP calculations without MS-MPI.)
'@
    }

    # --- main README content (double-quoted here-string => expands $variables we want) ---
    $content = @"
Quantum ESPRESSO for Windows (Intel oneAPI + MS-MPI) — Precompiled Binaries

Build date (UTC): $BuildDateUtc
Quantum ESPRESSO version/tag: $QeVersion

Why this exists (and why it matters)

Building Quantum ESPRESSO on Windows is significantly harder than on Linux/macOS:
toolchains, Fortran runtime DLLs, MPI, and dependency closure can be painful even for experts.

This distribution is a "just run it" package:
- All executables and required runtime DLLs are colocated in bin/
- No extra compiler/runtime installation needed on a typical Windows machine
- Signed executables (when built via CI release workflow) for safer downloading and sharing

Contents

$distLayoutSummary

Directory layout

- bin/
  Quantum ESPRESSO executables (*.exe), bundled runtime DLL dependencies, and (optionally) MS-MPI launcher.

- licenses/
  License texts for Quantum ESPRESSO and bundled third-party components.
  See licenses/THIRD_PARTY_NOTICES.txt for a summary.

- VERSION.txt
  Build metadata for traceability (exact build flags, toolchain versions, etc.)

Quick start (serial)

1) Open PowerShell or CMD.
2) cd into the bin\ folder:
     cd .\bin
3) Run an executable, e.g.:
     .\pw.exe -i scf-cg.in > scf-cg.out

$pathHint

MPI usage

$msmpiHint

Example (MPI): run from inside bin\:

  mpiexec -n 4 .\pw.exe -i scf-cg.in > scf-cg.out

OpenMP / threading

$ompHint

$mklHint

Tip: For hybrid MPI+OpenMP runs, avoid oversubscription.
Example: on a 16-core CPU, try mpiexec -n 4 with OMP_NUM_THREADS=4 as a starting point. Tune the above parameters to find the best performance for your calculations.

About QMatSuite (GUI)

If you want a more user-friendly workflow (modern user interface, project management, input generation, job management, result browsing),
check out QMatSuite — a modern GUI that can drive Quantum ESPRESSO and other engines. It is free and open source.

- Project home / downloads: www.qmatsuite.com
- Source code / releases: github.com/QMatSuite

This repo/toolchain not only produces ready-to-run executables, but also serves as the "high-performance, reproducible build backend" for QMatSuite,
and as a reference compilation recipe for anyone exploring QE/Wannier/etc. on Windows.

Licensing (important)

- Quantum ESPRESSO is licensed under GPL v2 or later.
- Bundled third-party runtime components (e.g., Intel oneAPI runtimes, Microsoft MPI, MSVC runtime DLLs)
  are redistributed under their respective licenses.
  See the licenses/ folder for the exact texts.

Disclaimer

These binaries are provided "as is" without warranty.

If you report an issue, please include:
- Windows version
- CPU model
- Whether you used MPI and/or OpenMP (and the values of OMP_NUM_THREADS / MKL_NUM_THREADS)
- The command line used (including mpiexec arguments)
- The contents of VERSION.txt
"@

    $content | Out-File -FilePath $readmeFile -Encoding UTF8 -NoNewline
    Write-Host "  Generated README.txt" -ForegroundColor Green
}

function Generate-Version {
    param(
        [string]$DistRoot,
        [string]$QeVersion,
        [string]$BuildDateUtc,
        [string]$BuildMode,
        [int]$ExeCount,
        [int]$DllCount,
        [string]$OneApiVersion,
        [string]$MsmpiVersion,

        # --- new optional fields (safe defaults) ---
        [bool]$LibxcEnabled = $false,
        [bool]$MpiEnabled   = $true,
        [bool]$OpenmpEnabled = $true,

        # Human-ish but still machine-readable backends
        [string]$FftBackend   = "Intel DFTI (MKL)",
        [string]$BlasBackend  = "Intel MKL",
        [string]$LapackBackend= "Intel MKL",

        # Tool versions (best-effort)
        [string]$CmakeVersion = $null,
        [string]$NinjaVersion = $null,

        # QE git metadata (best-effort)
        [string]$QeGitDescribe = $null,
        [string]$QeCommit = $null,

        # Traceability: toolchain repo commit (CI)
        [string]$ToolchainCommit = $env:GITHUB_SHA
    )

    $versionFile = Join-Path $DistRoot "VERSION.txt"

    function _norm([string]$s, [string]$fallback="unknown") {
        if ([string]::IsNullOrWhiteSpace($s)) { return $fallback }
        return $s.Trim()
    }

    # Best-effort discover versions if not provided
    if ([string]::IsNullOrWhiteSpace($CmakeVersion)) {
        try {
            $line = (& cmake --version 2>$null | Select-Object -First 1)
            if ($line -match 'cmake version\s+(.+)$') { $CmakeVersion = $matches[1].Trim() }
        } catch {}
    }
    if ([string]::IsNullOrWhiteSpace($NinjaVersion)) {
        try {
            $NinjaVersion = (& ninja --version 2>$null | Select-Object -First 1)
        } catch {}
    }

    $content = @"
build_date_utc=$BuildDateUtc
toolchain_repo=github.com/QMatSuite/QMatSuite
toolchain_commit=$(_norm $ToolchainCommit)

workflow=qe-windows-oneapi-msmpi-release
qe_version_or_commit=$(_norm $QeVersion)
qe_git_describe=$(_norm $QeGitDescribe)
qe_commit=$(_norm $QeCommit)

build_mode=$(_norm $BuildMode)
exe_count=$ExeCount
dll_count=$DllCount

oneapi_version=$(_norm $OneApiVersion)
msmpi_version=$(_norm $MsmpiVersion)

libxc=$([string]::new(@("OFF","ON")[[int]$LibxcEnabled]))
mpi=$([string]::new(@("OFF","ON")[[int]$MpiEnabled]))
openmp=$([string]::new(@("OFF","ON")[[int]$OpenmpEnabled]))

fft_backend=$($FftBackend)
blas_backend=$($BlasBackend)
lapack_backend=$($LapackBackend)

cmake_version=$(_norm $CmakeVersion)
ninja_version=$(_norm $NinjaVersion)
"@

    $content | Out-File -FilePath $versionFile -Encoding UTF8 -NoNewline
    Write-Host "  Generated VERSION.txt" -ForegroundColor Green
}



function Clean-SmokeTestFiles {
    param(
        [string]$DistRoot,
        [string]$DistBinDir  # Also check dist/bin for outdir
    )
    Write-Host "Cleaning up smoke test files..." -ForegroundColor Cyan
    
    $cleaned = @()
    
    # Remove outdir folder - check both dist/bin and dist root
    $outdirLocations = @(
        (Join-Path $DistBinDir "outdir"),
        (Join-Path $DistRoot "outdir")
    )
    
    foreach ($outdir in $outdirLocations) {
        if (Test-Path $outdir) {
            try {
                Remove-Item -Recurse -Force $outdir -ErrorAction Stop
                $cleaned += "outdir/ ($outdir)"
                Write-Host "  Removed: outdir/ ($outdir)" -ForegroundColor Gray
            } catch {
                Write-Warning "  Failed to remove outdir at $outdir : $_"
            }
        }
    }
    
    # Also search recursively for any outdir folders in dist/bin
    if (Test-Path $DistBinDir) {
        $allOutdirs = Get-ChildItem -Path $DistBinDir -Filter "outdir" -Directory -Recurse -ErrorAction SilentlyContinue
        foreach ($outdir in $allOutdirs) {
            if (Test-Path $outdir.FullName) {
                try {
                    Remove-Item -Recurse -Force $outdir.FullName -ErrorAction Stop
                    $cleaned += "outdir/ ($($outdir.FullName))"
                    Write-Host "  Removed: outdir/ ($($outdir.FullName))" -ForegroundColor Gray
                } catch {
                    Write-Warning "  Failed to remove outdir at $($outdir.FullName) : $_"
                }
            }
        }
    }
    
    # Remove .out files
    $outFiles = Get-ChildItem -Path $DistBinDir -Filter "*.out" -File -ErrorAction SilentlyContinue
    foreach ($file in $outFiles) {
        try {
            Remove-Item -Force $file.FullName -ErrorAction Stop
            $cleaned += $file.Name
            Write-Host "  Removed: $($file.Name)" -ForegroundColor Gray
        } catch {
            Write-Warning "  Failed to remove $($file.Name) : $_"
        }
    }
    
    # Remove .in files
    $inFiles = Get-ChildItem -Path $DistBinDir -Filter "*.in" -File -ErrorAction SilentlyContinue
    foreach ($file in $inFiles) {
        try {
            Remove-Item -Force $file.FullName -ErrorAction Stop
            $cleaned += $file.Name
            Write-Host "  Removed: $($file.Name)" -ForegroundColor Gray
        } catch {
            Write-Warning "  Failed to remove $($file.Name) : $_"
        }
    }
    
    # Remove .UPF files
    $upfFiles = Get-ChildItem -Path $DistBinDir -Filter "*.UPF" -File -ErrorAction SilentlyContinue
    foreach ($file in $upfFiles) {
        try {
            Remove-Item -Force $file.FullName -ErrorAction Stop
            $cleaned += $file.Name
            Write-Host "  Removed: $($file.Name)" -ForegroundColor Gray
        } catch {
            Write-Warning "  Failed to remove $($file.Name) : $_"
        }
    }
    
    # Remove deps-closure.txt from dist root
    $depsClosureFile = Join-Path $DistRoot "deps-closure.txt"
    if (Test-Path $depsClosureFile) {
        try {
            Remove-Item -Force $depsClosureFile -ErrorAction Stop
            $cleaned += "deps-closure.txt"
            Write-Host "  Removed: deps-closure.txt" -ForegroundColor Gray
        } catch {
            Write-Warning "  Failed to remove deps-closure.txt : $_"
        }
    }
    
    if ($cleaned.Count -gt 0) {
        Write-Host "  Cleaned $($cleaned.Count) smoke test file(s)/folder(s)" -ForegroundColor Green
    } else {
        Write-Host "  No smoke test files to clean" -ForegroundColor Gray
    }
}

function Verify-DistContents {
    param(
        [string]$DistRoot,
        [string]$DistBinDir
    )
    Write-Host "Verifying distribution contents..." -ForegroundColor Cyan
    
    # Check bin directory exists and has files
    if (-not (Test-Path $DistBinDir)) {
        throw "Distribution bin directory not found: $DistBinDir"
    }
    
    $binFiles = Get-ChildItem -Path $DistBinDir -File
    $exeFiles = $binFiles | Where-Object { $_.Extension -eq ".exe" }
    $dllFiles = $binFiles | Where-Object { $_.Extension -eq ".dll" }
    $otherFiles = $binFiles | Where-Object { $_.Extension -ne ".exe" -and $_.Extension -ne ".dll" }
    
    Write-Host "  bin/ contains: $($exeFiles.Count) .exe, $($dllFiles.Count) .dll" -ForegroundColor Green
    
    if ($exeFiles.Count -eq 0) {
        throw "No executables found in bin/ directory. Staging failed."
    }
    
    if ($otherFiles.Count -gt 0) {
        Write-Warning "bin/ contains non-exe/dll files:"
        foreach ($file in $otherFiles) {
            Write-Warning "  - $($file.Name)"
        }
    }
    
    # Check dist root for unexpected files (should only have bin/, licenses/, README.txt, VERSION.txt, deps-closure.txt)
    $rootFiles = Get-ChildItem -Path $DistRoot -File
    $expectedRootFiles = @("README.txt", "VERSION.txt", "deps-closure.txt")
    $unexpectedRootFiles = $rootFiles | Where-Object { $expectedRootFiles -notcontains $_.Name }
    
    if ($unexpectedRootFiles.Count -gt 0) {
        Write-Warning "Dist root contains unexpected files:"
        foreach ($file in $unexpectedRootFiles) {
            Write-Warning "  - $($file.Name)"
        }
    }
    
    # Check licenses directory exists
    $licensesDir = Join-Path $DistRoot "licenses"
    if (-not (Test-Path $licensesDir)) {
        Write-Warning "licenses/ directory not found"
    }
    
    Write-Host "  Verification complete" -ForegroundColor Green
    return @{
        ExeCount = $exeFiles.Count
        DllCount = $dllFiles.Count
        OtherFilesInBin = $otherFiles.Count
        UnexpectedRootFiles = $unexpectedRootFiles.Count
    }
}

function Show-PostBuildDiagnostics {
    param(
        [string]$BuildBin,
        [string]$BuildDir
    )
    
    Write-Host "=== Post-build diagnostics (from staging) ===" -ForegroundColor Cyan
    
    # Find pw.exe
    $pwExe = Join-Path $BuildBin "pw.exe"
    if (-not (Test-Path $pwExe)) {
        $pwExe = Join-Path $BuildBin "pw.x.exe"
    }
    
    if (-not (Test-Path $pwExe)) {
        Write-Warning "pw.exe not found for stack reserve check (searched: $BuildBin\pw.exe, $BuildBin\pw.x.exe)"
    } else {
        Write-Host "pw.exe path: $pwExe"
        
        # Resolve dumpbin using existing logic
        $dumpbinPath = Ensure-DumpbinInPath
        if (-not $dumpbinPath) {
            Write-Warning "dumpbin.exe not found, skipping stack reserve check"
        } else {
            Write-Host "Stack reserve size:"
            try {
                $line = Get-StackReserveLine -DumpbinPath $dumpbinPath -PwExe $pwExe
                if ($line) {
                    Write-Host "  $line"
                    
                    # Parse hex value and convert to MB
                    # Pattern matches: optional whitespace, hex digits, whitespace, "size of stack reserve"
                    if ($line -match "\s*([0-9A-Fa-f]+)\s+size of stack reserve") {
                        $hexValue = $Matches[1]
                        try {
                            $bytes = [Convert]::ToInt64($hexValue, 16)
                            $mb = [math]::Round($bytes / 1MB, 2)
                            Write-Host "  Interpreted: $mb MB ($hexValue hex = $bytes bytes)"
                        } catch {
                            Write-Host "  (Could not parse hex value: $hexValue)"
                        }
                    }
                } else {
                    Write-Warning "  Stack reserve line not found in dumpbin output"
                }
            } catch {
                Write-Warning "  Failed to run dumpbin: $_"
            }
        }
    }
    
    # Check CMakeCache.txt - build directory is parent of bin directory
    $buildDirPath = Split-Path $BuildBin -Parent
    $cmakeCache = Join-Path $buildDirPath "CMakeCache.txt"
    
    if (-not (Test-Path $cmakeCache)) {
        Write-Warning "CMakeCache.txt not found at: $cmakeCache"
    } else {
        Write-Host "Relevant CMakeCache flags:"
        $flags = Select-String -Path $cmakeCache -Pattern "CMAKE_Fortran_FLAGS_RELEASE:STRING=|CMAKE_EXE_LINKER_FLAGS:STRING=" | ForEach-Object { $_.Line }
        if ($flags) {
            foreach ($flag in $flags) {
                Write-Host "  $flag"
            }
        } else {
            Write-Warning "  No matching flags found in CMakeCache.txt"
        }
    }
    
    Write-Host ""
}

function Run-Probe {
    param(
        [string]$ExePath
    )
    if ($WhatIfPreference) { Write-Host "Probe skipped (WhatIf)."; return }
    if (-not (Test-Path $ExePath)) {
        Write-Warning "Probe exe not found: $ExePath"
        return
    }
    $tmpOut = [System.IO.Path]::GetTempFileName()
    $tmpErr = [System.IO.Path]::GetTempFileName()
    $tmpInput = $null
    try {
        $exeName = [System.IO.Path]::GetFileName($ExePath).ToLowerInvariant()
        $isPwExe = $exeName -eq "pw.exe" -or $exeName -eq "pw.x.exe"
        
        if ($isPwExe) {
            # For pw.exe, use a real input file to test runtime dependencies
            $tmpInput = [System.IO.Path]::GetTempFileName()
            $inputContent = @"
&control
   calculation = 'scf'
   pseudo_dir = '.'
   outdir = './outdir'
/
&system
   ibrav=2, celldm(1) =10.20, 
   nat=2, ntyp=1,
   ecutwfc=12.0
/
&electrons
   diago_thr_init = 1.D-5
   diagonalization='cg'
/
ATOMIC_SPECIES
 Si  28.086  Si.pz-vbc.UPF
ATOMIC_POSITIONS (alat)
 Si 0.00 0.00 0.00
 Si 0.25 0.25 0.25
K_POINTS
  2 
   0.250000  0.250000  0.250000   1.00
   0.250000  0.250000  0.750000   3.00
"@
            $inputContent | Out-File -FilePath $tmpInput -Encoding ASCII
            $psi = New-Object System.Diagnostics.ProcessStartInfo
            $psi.FileName = $ExePath
            $psi.Arguments = "-i `"$tmpInput`""
        } else {
            # For other executables, use -help flag
            $psi = New-Object System.Diagnostics.ProcessStartInfo
            $psi.FileName = $ExePath
            $psi.Arguments = "-help"
        }
        
        $psi.RedirectStandardOutput = $true
        $psi.RedirectStandardError = $true
        $psi.RedirectStandardInput = $true  # Redirect stdin so we can close it to prevent hanging on interactive prompts
        $psi.UseShellExecute = $false
        $proc = New-Object System.Diagnostics.Process
        $proc.StartInfo = $psi
        $proc.Start() | Out-Null
        # Close stdin immediately to send EOF - prevents executables that wait for input from hanging
        $proc.StandardInput.Close()
        
        # Wait for process with timeout (10 seconds) - some QE utilities may not exit on -help or may take time to process input
        $timeoutMs = 10000
        if (-not $proc.WaitForExit($timeoutMs)) {
            Write-Warning "Probe timed out after ${timeoutMs}ms for $ExePath (process may be waiting for input or taking too long)"
            $proc.Kill()
            $proc.WaitForExit() | Out-Null  # Ensure cleanup after kill
            return
        }
        
        # Read output after process has exited (safe since stdin is closed and process terminated)
        $outText = $proc.StandardOutput.ReadToEnd()
        $errText = $proc.StandardError.ReadToEnd()
        $rc = $proc.ExitCode
        $errCombined = "$outText`n$errText"
        if ($rc -ne 0 -and ($errCombined -match "(was not found|cannot proceed because|missing)")) {
            Write-Warning "Probe indicates missing dependency for $ExePath (exit $rc):"
            Write-Warning $errCombined
        } else {
            Write-Host "Probe finished for $ExePath (exit $rc)."
        }
    } finally {
        Remove-Item $tmpOut,$tmpErr -ErrorAction SilentlyContinue
        if ($tmpInput) {
            Remove-Item $tmpInput -ErrorAction SilentlyContinue
        }
    }
}

#---------------- Main ----------------#
# Get script directory - handle both direct execution and dot-sourcing
$scriptPath = $null
if ($MyInvocation.MyCommand.Path) {
    $scriptPath = $MyInvocation.MyCommand.Path
} elseif ($PSCommandPath) {
    $scriptPath = $PSCommandPath
} elseif ($MyInvocation.PSCommandPath) {
    $scriptPath = $MyInvocation.PSCommandPath
}

if (-not $scriptPath) {
    throw "Cannot determine script path. Please run this script directly, not via dot-sourcing or other indirect methods."
}
if (-not (Test-Path $scriptPath)) {
    Write-Warning "Script path determined but file does not exist: $scriptPath. Continuing anyway..."
}

$ScriptDir = Split-Path -Parent $scriptPath
if (-not $ScriptDir) {
    throw "Cannot determine script directory from path: $scriptPath"
}

# Build paths with null checks
$RepoRootPath = Join-Path $ScriptDir ".."
if (-not $RepoRootPath) {
    throw "Cannot build RepoRoot path from ScriptDir: $ScriptDir"
}
$RepoRoot = Resolve-Path $RepoRootPath -ErrorAction SilentlyContinue
if (-not $RepoRoot) { 
    $RepoRoot = $RepoRootPath 
}

# Auto-detect build directory if using default
$actualBuildDir = Find-BuildDirectory -RepoRoot $RepoRoot -DefaultBuildDir $BuildDir
if ($actualBuildDir -ne $BuildDir) {
    Write-Host "Auto-detected build directory: $actualBuildDir (was: $BuildDir)" -ForegroundColor Cyan
    $BuildDir = $actualBuildDir
}


# Detect MPI configuration from build
$mpiConfig = Detect-MPIFromBuild -BuildDir $BuildDir -RepoRoot $RepoRoot
if ($mpiConfig) {
    Write-Host "=== MPI Detection ===" -ForegroundColor Cyan
    Write-Host "Detected: $mpiConfig" -ForegroundColor $(if ($mpiConfig -eq "no-MPI") { "Yellow" } else { "Green" })
    Write-Host ""
} else {
    Write-Host "=== MPI Detection ===" -ForegroundColor Cyan
    Write-Host "Could not detect MPI configuration (CMakeCache.txt not found or invalid)" -ForegroundColor Yellow
    Write-Host ""
}

# Build BuildBin path
$BuildBinPath = Join-Path $RepoRoot (Join-Path $BuildDir "bin")
if (-not $BuildBinPath) {
    throw "Cannot build BuildBin path from RepoRoot: $RepoRoot, BuildDir: $BuildDir"
}
$BuildBin = Resolve-Path $BuildBinPath -ErrorAction SilentlyContinue
if (-not $BuildBin) {
    $BuildBin = $BuildBinPath
}

# Build DistRoot path
$DistRootPath = Join-Path $RepoRoot $DistDir
if (-not $DistRootPath) {
    throw "Cannot build DistRoot path from RepoRoot: $RepoRoot, DistDir: $DistDir"
}
$DistRoot = Resolve-Path $DistRootPath -ErrorAction SilentlyContinue
if (-not $DistRoot) { 
    $DistRoot = $DistRootPath 
}

# oneAPI root discovery
$oneApiRoot = $env:ONEAPI_ROOT
if (-not $oneApiRoot -or -not (Test-Path $oneApiRoot)) {
    $oneApiRoot = "C:\Program Files (x86)\Intel\oneAPI"
}
if (-not (Test-Path $oneApiRoot)) {
    throw "ONEAPI_ROOT not found. Set ONEAPI_ROOT or install Intel oneAPI (expected at '$oneApiRoot')."
}

if (-not (Test-Path $BuildBin)) {
    throw "QE build bin not found at '$BuildBin'. Build first or set -BuildDir."
}

Write-Host "QE build bin : $BuildBin"
Write-Host "Dist target  : $DistRoot"
Write-Host "ONEAPI_ROOT  : $oneApiRoot"

# Show post-build diagnostics
Show-PostBuildDiagnostics -BuildBin $BuildBin -BuildDir $BuildDir

if ($Clean -and (Test-Path $DistRoot)) {
    if ($PSCmdlet.ShouldProcess($DistRoot, "Clean dist folder") -and -not $WhatIf) {
        Remove-Item -Recurse -Force $DistRoot
    }
}
if (-not (Test-Path $DistRoot) -and -not $WhatIf) {
    New-Item -ItemType Directory -Force -Path $DistRoot | Out-Null
}

# Create dist/bin directory structure
$DistBinDir = Join-Path $DistRoot "bin"
if (-not (Test-Path $DistBinDir) -and -not $WhatIf) {
    New-Item -ItemType Directory -Force -Path $DistBinDir | Out-Null
}


# Stage executables to dist/bin
Write-Host "=== Staging executables ===" -ForegroundColor Cyan
$stagedExeFiles = Stage-Executables -BinDir $BuildBin -DistBinDir $DistBinDir -Only $Only -WhatIf:$WhatIf
if (-not $DistBinDir) {
    throw "DistBinDir is null. Cannot stage executables."
}
# Construct paths using final staged names (not original build names)
# For .x.exe files, they are renamed to .exe in dist/bin, so we need to use the final name
$stagedExePaths = $stagedExeFiles | ForEach-Object { 
    if ($_ -and $_.Name) {
        # If it's a .x.exe file, use the renamed .exe version in dist/bin
        if ($_.Name -like "*.x.exe") {
            $finalName = $_.Name -replace '\.x\.exe$', '.exe'
            Join-Path $DistBinDir $finalName
        } else {
            # Regular .exe files keep their name
            Join-Path $DistBinDir $_.Name
        }
    }
} | Where-Object { $_ -ne $null }

# --- libxc auto-stage (based on QE_ENABLE_LIBXC in CMakeCache.txt) ---
$buildPath = Join-Path $RepoRoot $BuildDir
$cmakeCache = Join-Path $buildPath "CMakeCache.txt"
$libxcEnabled = Get-CMakeCacheBool -CachePath $cmakeCache -Key "QE_ENABLE_LIBXC"

if ($libxcEnabled) {
    Write-Host "libxc on (QE_ENABLE_LIBXC=ON) -> will stage libxc runtime bin/*" -ForegroundColor Green

    $libxcBin = Join-Path $RepoRoot "upstream\libxc\install\bin"
    if (-not (Test-Path $libxcBin)) {
        throw "QE_ENABLE_LIBXC=ON but libxc bin dir not found: $libxcBin"
    }

    # Dist bin dir (use your existing dist root; keep consistent with your script)
    New-Item -ItemType Directory -Force -Path $DistBinDir | Out-Null

    Copy-Item -Force -Path (Join-Path $libxcBin "*") -Destination $DistBinDir


    # Also stage libxc license (prefer COPYING, fallback to LICENSE*)
    $libxcLicenseCandidates = @(
        (Join-Path $RepoRoot "upstream\libxc\src\COPYING"),
        (Join-Path $RepoRoot "upstream\libxc\src\COPYING.txt"),
        (Join-Path $RepoRoot "upstream\libxc\src\LICENSE"),
        (Join-Path $RepoRoot "upstream\libxc\src\LICENSE.txt"),
        (Join-Path $RepoRoot "upstream\libxc\src\LICENSE.md"),
        (Join-Path $RepoRoot "upstream\libxc\src\LICENCE"),       # UK spelling (rare)
        (Join-Path $RepoRoot "upstream\libxc\src\LICENCE.txt")
    )

    $libxcLicenseSrc = $null
    foreach ($p in $libxcLicenseCandidates) {
        if (Test-Path $p) { $libxcLicenseSrc = $p; break }
    }

    # If not found at repo root, try anywhere under upstream/libxc/src as a last resort
    if (-not $libxcLicenseSrc) {
        $srcRoot = Join-Path $RepoRoot "upstream\libxc\src"
        if (Test-Path $srcRoot) {
            $hit = Get-ChildItem $srcRoot -Recurse -File -ErrorAction SilentlyContinue |
                Where-Object { $_.Name -match '^(COPYING|LICENSE|LICENCE)(\..*)?$' } |
                Select-Object -First 1
            if ($hit) { $libxcLicenseSrc = $hit.FullName }
        }
    }

    if ($libxcLicenseSrc) {
        $libxcLicenseDstDir = Join-Path $distRoot "licenses\quantum-espresso\external\libxc"
        New-Item -ItemType Directory -Force -Path $libxcLicenseDstDir | Out-Null

        # Keep original filename to avoid ambiguity
        $dstName = Split-Path -Leaf $libxcLicenseSrc
        Copy-Item -Force -Path $libxcLicenseSrc -Destination (Join-Path $libxcLicenseDstDir $dstName)

        Write-Host "Staged libxc license: $dstName" -ForegroundColor Green
    } else {
        Write-Host "WARNING: QE_ENABLE_LIBXC=ON but no libxc license file found under upstream\libxc\src (COPYING/LICENSE/LICENCE). Skipping." -ForegroundColor Yellow
        # If you prefer hard fail instead of warning, replace the line above with:
        # throw "QE_ENABLE_LIBXC=ON but no libxc license found under upstream\libxc\src"
    }

}
elseif ($libxcEnabled -eq "OFF") {
    Write-Host "libxc off (QE_ENABLE_LIBXC=OFF) -> skip libxc staging" -ForegroundColor DarkGray
}
else {
    Write-Host "QE_ENABLE_LIBXC not found in CMakeCache.txt -> skip libxc staging" -ForegroundColor Yellow
}
# --- end libxc auto-stage ---

# Layer 1: Recursive dumpbin dependency closure (static transitive closure)
$dumpbinPath = Ensure-DumpbinInPath
$searchRoots = Get-SearchRoots -OneApiRoot $oneApiRoot
$closureFile = Join-Path $DistRoot "deps-closure.txt"
$closureSummary = $null
if ($dumpbinPath) {
    Write-Host "=== Layer 1: Recursive dependency closure ===" -ForegroundColor Cyan
    Write-Host "Using dumpbin at: $dumpbinPath"
    $closureSummary = Copy-DllClosure -ExePaths $stagedExePaths -DumpbinPath $dumpbinPath -SearchRoots $searchRoots -DistBinDir $DistBinDir -BuildBinDir $BuildBin -VerboseDeps:$VerboseDeps -WhatIf:$WhatIf
} else {
    Write-Warning "dumpbin.exe not found; skipping recursive closure. Staging may be incomplete."
    # Create empty closure file with note
    if (-not $WhatIfPreference) {
        "# Dependency closure skipped: dumpbin.exe not found" | Out-File -FilePath $closureFile -Encoding UTF8
    }
    $closureSummary = @{
        ScannedCount = 0
        CopiedCount = 0
        ResolvedCount = 0
        UnresolvedCount = 0
        UnresolvedDlls = @()
        ElapsedSeconds = 0
    }
}

# Layer 2: Must-have DLL patterns (runtime safety net)
Write-Host "=== Layer 2: Must-have DLL patterns ===" -ForegroundColor Cyan
$mklRoot = $env:MKLROOT
if (-not $mklRoot -or -not (Test-Path $mklRoot)) {
    $mklRoot = Join-Path $oneApiRoot "mkl\latest"
}
Copy-MustHaveDllPatterns -SearchRoots $searchRoots -DistDir $DistBinDir -OneApiRoot $oneApiRoot -MklRoot $mklRoot

# Stage MS-MPI runtime if MSMPI is detected
$hasMsmpi = $false
if ($mpiConfig -eq "MSMPI") {
    Stage-MSMPIRuntime -DistDir $DistBinDir
    $hasMsmpi = $true
}

# Layer 3: Clean-room smoke test (final acceptance)
Write-Host "=== Layer 3: Clean-room smoke test ===" -ForegroundColor Cyan
# Calculate resources directory from script directory (go up 3 levels: scripts -> oneapi -> windows -> quantum_espresso)
$quantumEspressoRoot = Resolve-Path (Join-Path $ScriptDir "..\..\..") -ErrorAction Stop
$resourcesDir = Join-Path $quantumEspressoRoot "tests\resources"
# Run smoke test in strict mode by default (require exit code 0 AND "JOB DONE")
$strictMode = -not $NoStrict
if ($mpiConfig -eq "MSMPI") {
    $mpiexecPath = Join-Path $DistBinDir "mpiexec.exe"
    Run-SmokeTest -DistDir $DistBinDir -ResourcesDir $resourcesDir -Strict:$strictMode -UseMpi -MpiLauncher $mpiexecPath -MpiNP 2
} else {
    Run-SmokeTest -DistDir $DistBinDir -ResourcesDir $resourcesDir -Strict:$strictMode
}

# Clean up smoke test files (optional, enabled by default)
if (-not $NoCleanSmokeTest) {
    Clean-SmokeTestFiles -DistRoot $DistRoot -DistBinDir $DistBinDir
} else {
    Write-Host "Skipping smoke test file cleanup (NoCleanSmokeTest specified)" -ForegroundColor Yellow
}

# Package licenses and notices
Write-Host "=== Packaging licenses and notices ===" -ForegroundColor Cyan
$licensesDir = Join-Path $DistRoot "licenses"
if (-not (Test-Path $licensesDir) -and -not $WhatIf) {
    New-Item -ItemType Directory -Force -Path $licensesDir | Out-Null
}

# Get QE source directory
# BuildBin is typically: .../upstream/qe/build-win-oneapi/bin
# So we need to go up 2 levels to get to .../upstream/qe
$buildDir = Split-Path $BuildBin -Parent  # .../upstream/qe/build-win-oneapi
$qeSourceDir = Split-Path $buildDir -Parent  # .../upstream/qe
if (-not (Test-Path $qeSourceDir)) {
    # Try alternative path
    $qeSourceDir = Join-Path $RepoRoot "upstream\qe"
}
if (-not (Test-Path $qeSourceDir)) {
    Write-Warning "QE source directory not found at: $qeSourceDir. License copying may fail."
}

# Get versions (pass DistBinDir for best-effort version extraction from bundled DLLs)
$buildDateUtc = (Get-Date).ToUniversalTime().ToString("yyyy-MM-dd HH:mm:ss UTC")
$qeVersion = Get-QeVersion -QeSourceDir $qeSourceDir
$oneApiVersion = Get-OneApiVersion -OneApiRoot $oneApiRoot -DistBinDir $DistBinDir
$msmpiVersion = Get-MsmpiVersion -DistBinDir $DistBinDir

# Copy license files
$qeLicenseResult = Copy-QeLicenses -QeSourceDir $qeSourceDir -LicensesDir $licensesDir
$oneApiLicenseResult = Copy-OneApiLicenses -OneApiRoot $oneApiRoot -LicensesDir $licensesDir -RepoRootPath $RepoRootPath
$msmpiLicenseResult = Copy-MsmpiLicense -LicensesDir $licensesDir

# Generate THIRD_PARTY_NOTICES.txt
Generate-ThirdPartyNotices -LicensesDir $licensesDir -QeVersion $qeVersion -OneApiVersion $oneApiVersion -MsmpiVersion $msmpiVersion -BuildDateUtc $buildDateUtc

# Detect build mode (try to infer from exe count or environment)
$buildMode = if ($env:BUILD_MODE) { $env:BUILD_MODE } else { if ($stagedExeFiles.Count -gt 50) { "all" } else { "pw" } }

# Get file counts for verification (must match what's actually in dist/bin)
$binFiles = if (Test-Path $DistBinDir) { Get-ChildItem -Path $DistBinDir -File } else { @() }
$exeCount = ($binFiles | Where-Object { $_.Extension -eq ".exe" }).Count
$dllCount = ($binFiles | Where-Object { $_.Extension -eq ".dll" }).Count

# Check if MS-MPI launcher is actually bundled
$hasMsmpiLauncher = (Test-Path (Join-Path $DistBinDir "mpiexec.exe")) -or (Test-Path (Join-Path $DistBinDir "msmpiexec.exe"))

# Generate README.txt and VERSION.txt (using actual counts from dist/bin)
$libxcEnabled = Get-CMakeCacheBool -CachePath $cmakeCache -Key "QE_ENABLE_LIBXC"
$mpiEnabled = Get-CMakeCacheBool -CachePath $cmakeCache -Key "QE_ENABLE_MPI"
$openmpEnabled = Get-CMakeCacheBool -CachePath $cmakeCache -Key "QE_ENABLE_OPENMP"
$qeGitDescribe = ""
$qeCommit = ""
try {
  $qeCommit = (& git -C $QeSourceRoot rev-parse HEAD 2>$null).Trim()
  $qeGitDescribe = (& git -C $QeSourceRoot describe --tags --always --dirty 2>$null).Trim()
} catch {}

Generate-Readme -DistRoot $DistRoot -QeVersion $qeVersion -BuildDateUtc $buildDateUtc -ExeCount $exeCount -DllCount $dllCount -HasMsmpiLauncher $hasMsmpiLauncher
Generate-Version -DistRoot $DistRoot -QeVersion $qeVersion -BuildDateUtc $buildDateUtc `
    -BuildMode $buildMode -ExeCount $exeCount -DllCount $dllCount -OneApiVersion $oneApiVersion `
    -MsmpiVersion $msmpiVersion -LibxcEnabled $libxcEnabled -MpiEnabled $mpiEnabled `
    -OpenmpEnabled $openmpEnabled -QeGitDescribe $qeGitDescribe -QeCommit $qeCommit

# Verify distribution contents
$verifyResult = Verify-DistContents -DistRoot $DistRoot -DistBinDir $DistBinDir

# Print final summary
Write-Host ""
Write-Host "=== Staging Summary ===" -ForegroundColor Green
Write-Host "  Executables staged: $exeCount" -ForegroundColor Green
Write-Host "  Runtime DLLs bundled: $dllCount" -ForegroundColor Green
if ($closureSummary) {
    Write-Host "  Unique modules scanned: $($closureSummary.ScannedCount)" -ForegroundColor Green
    # Note: Unresolved DLLs are mostly Windows system DLLs (not redistributable), which is expected
    Write-Host "  Unresolved DLLs: $($closureSummary.UnresolvedCount) (mostly Windows system DLLs, expected)" -ForegroundColor $(if ($closureSummary.UnresolvedCount -eq 0) { "Green" } else { "Yellow" })
}
Write-Host "  License directories:" -ForegroundColor Green
Write-Host "    - quantum-espresso: $(if ($qeLicenseResult.Copied.Count -gt 0) { 'OK' } else { 'MISSING' })" -ForegroundColor $(if ($qeLicenseResult.Copied.Count -gt 0) { "Green" } else { "Red" })
Write-Host "    - intel-oneapi: $(if ($oneApiLicenseResult.Copied.Count -gt 0) { 'OK' } else { 'WARNING' })" -ForegroundColor $(if ($oneApiLicenseResult.Copied.Count -gt 0) { "Green" } else { "Yellow" })
Write-Host "    - microsoft-msmpi: $(if ($msmpiLicenseResult.Copied.Count -gt 0) { 'OK' } else { 'WARNING' })" -ForegroundColor $(if ($msmpiLicenseResult.Copied.Count -gt 0) { "Green" } else { "Yellow" })
Write-Host ""
Write-Host "Staging complete. Distribution in: $DistRoot" -ForegroundColor Green


