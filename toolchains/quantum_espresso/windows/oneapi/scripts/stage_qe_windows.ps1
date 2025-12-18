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
#>
[CmdletBinding(SupportsShouldProcess = $true)]
param(
    [string]$BuildDir = "upstream/qe/build-win-oneapi",
    [string]$DistDir  = "dist/win-oneapi",
    [string[]]$Only,
    [switch]$Clean,
    [switch]$VerboseDeps,
    [string]$ProbeExe,
    [switch]$NoStrict  # If set, disables strict mode (not recommended)
)

$ErrorActionPreference = 'Stop'

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
        if ($line -match "Image has the following dependencies") { $capture = $true; continue }
        if ($capture -and [string]::IsNullOrWhiteSpace($line)) { $capture = $false }
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
        [string[]]$PathDirs
    )
    foreach ($root in $SearchRoots) {
        $found = Get-ChildItem -Path $root -Filter $Name -File -Recurse -ErrorAction SilentlyContinue | Select-Object -First 1
        if ($found) { return $found.FullName }
    }
    foreach ($dir in $PathDirs) {
        $candidate = Join-Path $dir $Name
        if (Test-Path $candidate) { return (Resolve-Path $candidate).Path }
    }
    return $null
}

function Stage-Executables {
    param(
        [string]$BinDir,
        [string]$DistDir,
        [string[]]$Only
    )
    $exeFiles = Get-ChildItem -Path $BinDir -Filter *.exe -File
    if ($Only) {
        $normalized = $Only | ForEach-Object { if ($_ -like "*.exe") { $_ } else { "$_.exe" } }
        $exeFiles = $exeFiles | Where-Object { $normalized -contains $_.Name }
    }
    if (-not $exeFiles) {
        throw "No executables found to stage. Check build output or -Only filters."
    }
    
    # Separate .x.exe files from regular .exe files
    $xExeFiles = $exeFiles | Where-Object { $_.Name -like "*.x.exe" }
    $regularExeFiles = $exeFiles | Where-Object { $_.Name -notlike "*.x.exe" }
    
    # Process .x.exe files: remove existing .exe, copy .x.exe, then rename to .exe
    foreach ($xExe in $xExeFiles) {
        $targetName = $xExe.Name -replace '\.x\.exe$', '.exe'
        $targetPath = Join-Path $DistDir $targetName
        $xExeDest = Join-Path $DistDir $xExe.Name
        
        if ($PSCmdlet.ShouldProcess($targetPath, "Stage $($xExe.Name) -> $targetName") -and -not $WhatIfPreference) {
            # Remove existing .exe if it exists
            if (Test-Path $targetPath) {
                Remove-Item -Path $targetPath -Force -ErrorAction SilentlyContinue
            }
            # Remove any existing .x.exe in dist (shouldn't exist, but clean up just in case)
            if (Test-Path $xExeDest) {
                Remove-Item -Path $xExeDest -Force -ErrorAction SilentlyContinue
            }
            # Copy .x.exe to dist
            Copy-Item -Force $xExe.FullName $xExeDest
            # Rename to .exe (removes .x.exe from dist)
            Rename-Item -Path $xExeDest -NewName $targetName -Force
        }
        Write-Host "  exe: $($xExe.Name) -> $targetName"
    }
    
    # Process regular .exe files (not .x.exe): just copy them
    foreach ($exe in $regularExeFiles) {
        $dest = Join-Path $DistDir $exe.Name
        if ($PSCmdlet.ShouldProcess($dest, "Copy exe $($exe.FullName)") -and -not $WhatIfPreference) {
            Copy-Item -Force $exe.FullName $dest
        }
        Write-Host "  exe: $($exe.Name)"
    }
    
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

function Copy-DllClosure {
    param(
        [string[]]$ExePaths,
        [string]$DumpbinPath,
        [string[]]$SearchRoots,
        [string]$DistDir,
        [switch]$VerboseDeps
    )
    $pathDirs = Get-ValidPathDirs
    $queue = New-Object System.Collections.Queue
    $visited = @{}  # Case-insensitive visited set by basename
    $missing = [System.Collections.Generic.List[string]]::new()
    $copied = [System.Collections.Generic.List[object]]::new()
    $resolved = [System.Collections.Generic.List[object]]::new()  # Track all resolved DLLs for closure file
    $scanCount = 0
    $MAX_SCAN_COUNT = 500
    $DUMPBIN_TIMEOUT_MS = 5000

    Write-Host "Starting recursive dependency closure (max ${MAX_SCAN_COUNT} modules, ${DUMPBIN_TIMEOUT_MS}ms per dumpbin)..."
    Write-Host "  Initial executables to scan: $($ExePaths.Count)"

    # Enqueue initial executables
    $enqueuedCount = 0
    foreach ($exe in $ExePaths) {
        if (-not (Test-Path $exe)) { 
            Write-Warning "  Executable not found: $exe"
            continue 
        }
        $queue.Enqueue($exe)
        $enqueuedCount++
    }
    Write-Host "  Enqueued $enqueuedCount executable(s) for scanning"

    while ($queue.Count -gt 0) {
        if ($scanCount -ge $MAX_SCAN_COUNT) {
            throw "Dependency closure exceeded maximum scan count (${MAX_SCAN_COUNT}). Possible resolution loop. Check for circular dependencies or increase MAX_SCAN_COUNT."
        }
        
        $filePath = $queue.Dequeue()
        $scanCount++
        
        # Get basename for visited check (case-insensitive)
        $basename = [System.IO.Path]::GetFileName($filePath).ToLowerInvariant()
        if ($visited.ContainsKey($basename)) { continue }
        $visited[$basename] = $true

        # Get dependencies
        $deps = Get-ExeDependentsDumpbin -Path $filePath -DumpbinPath $DumpbinPath -Verbose:$VerboseDeps -TimeoutMs $DUMPBIN_TIMEOUT_MS
        if ($VerboseDeps) {
            Write-Host "  Scanned: $basename -> $($deps.Count) dependencies"
        }

        foreach ($dllName in $deps) {
            $dllBasename = $dllName.ToLowerInvariant()
            if ($visited.ContainsKey($dllBasename)) { continue }

            # Skip denylisted DLLs
            if (Is-DenylistedDll -Name $dllName) {
                Write-Host "  DENYLIST: Skipping $dllName (prohibited from distribution)"
                continue
            }

            $resolvedPath = Resolve-DllPath -Name $dllName -SearchRoots $SearchRoots -PathDirs $pathDirs
            if (-not $resolvedPath) {
                if (-not (Is-SystemDll -Name $dllName -ResolvedPath $null)) {
                    $missing.Add($dllName)
                }
                continue
            }
            if (Is-SystemDll -Name $dllName -ResolvedPath $resolvedPath) { continue }

            # Track all resolved DLLs (for closure file)
            $resolved.Add(@{
                Basename = $dllBasename
                FullPath = $resolvedPath
            })

            # Copy DLL
            $dest = Join-Path $DistDir (Split-Path $resolvedPath -Leaf)
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

            # Enqueue for recursive scanning
            $queue.Enqueue($resolvedPath)
        }
    }

    # Write deps-closure.txt (include all resolved DLLs, not just newly copied)
    $closureFile = Join-Path $DistDir "deps-closure.txt"
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
                    "# Dependency closure: resolved list had $($resolved.Count) items but no valid lines generated`n# Scanned $scanCount modules." | Out-File -FilePath $closureFile -Encoding UTF8 -Force
                }
            } else {
                # Write file with header comment if no DLLs were resolved
                "# Dependency closure scan completed. No non-system DLLs were resolved via dumpbin.`n# Scanned $scanCount modules.`n# Executables scanned: $($ExePaths.Count)" | Out-File -FilePath $closureFile -Encoding UTF8 -Force
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

    Write-Host "Dependency closure complete: $($copied.Count) DLLs copied, $($resolved.Count) DLLs resolved, $scanCount modules scanned"
    if ($missing.Count -gt 0) {
        Write-Warning ("Missing non-system DLLs: " + ($missing | Sort-Object -Unique -join ", "))
        Write-Warning "Ensure oneAPI redistributables are installed and ONEAPI_ROOT is correct; re-run setvars if needed."
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

function Run-SmokeTest {
    param(
        [string]$DistDir,
        [string]$ResourcesDir,
        [switch]$Strict  # Strict mode: fail if exit code != 0 OR "JOB DONE" not detected
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
    
    $smokeDir = $DistDir  # Run test in dist root (where pw.exe is)
    $smokeInput = Join-Path $smokeDir "scf-cg.in"
    $smokeOutput = Join-Path $smokeDir "scf-cg.out"
    $pseudoFile = Join-Path $smokeDir "Si.pz-vbc.UPF"
    $binDir = $DistDir  # Executables are in dist root
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
    
    if (-not $WhatIfPreference) {
        # Use cmd.exe to run with shell redirection (avoids deadlock issues with large output)
        # This matches the manual command: pw.exe -i scf-cg.in > scf-cg.out
        $cmdLine = "`"$pwExe`" -i scf-cg.in > scf-cg.out 2>&1"
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
        Write-Host "  (WhatIf: would run $pwExe -i scf-cg.in > scf-cg.out with clean PATH)"
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
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$RepoRoot  = Resolve-Path (Join-Path $ScriptDir "..")
$BuildBin  = Resolve-Path (Join-Path $RepoRoot (Join-Path $BuildDir "bin")) -ErrorAction SilentlyContinue
$DistRoot  = Resolve-Path (Join-Path $RepoRoot $DistDir) -ErrorAction SilentlyContinue
if (-not $DistRoot) { $DistRoot = Join-Path $RepoRoot $DistDir }

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

# Stage executables
Write-Host "Copying executables..."
$stagedExeFiles = Stage-Executables -BinDir $BuildBin -DistDir $DistRoot -Only $Only -WhatIf:$WhatIf
$stagedExePaths = $stagedExeFiles | ForEach-Object { Join-Path $DistRoot $_.Name }

# Layer 1: Recursive dumpbin dependency closure (static transitive closure)
$dumpbinPath = Ensure-DumpbinInPath
$searchRoots = Get-SearchRoots -OneApiRoot $oneApiRoot
$closureFile = Join-Path $DistRoot "deps-closure.txt"
if ($dumpbinPath) {
    Write-Host "=== Layer 1: Recursive dependency closure ===" -ForegroundColor Cyan
    Write-Host "Using dumpbin at: $dumpbinPath"
    Copy-DllClosure -ExePaths $stagedExePaths -DumpbinPath $dumpbinPath -SearchRoots $searchRoots -DistDir $DistRoot -VerboseDeps:$VerboseDeps -WhatIf:$WhatIf
} else {
    Write-Warning "dumpbin.exe not found; skipping recursive closure. Staging may be incomplete."
    # Create empty closure file with note
    if (-not $WhatIfPreference) {
        "# Dependency closure skipped: dumpbin.exe not found" | Out-File -FilePath $closureFile -Encoding UTF8
    }
}

# Layer 2: Must-have DLL patterns (runtime safety net)
Write-Host "=== Layer 2: Must-have DLL patterns ===" -ForegroundColor Cyan
$mklRoot = $env:MKLROOT
if (-not $mklRoot -or -not (Test-Path $mklRoot)) {
    $mklRoot = Join-Path $oneApiRoot "mkl\latest"
}
Copy-MustHaveDllPatterns -SearchRoots $searchRoots -DistDir $DistRoot -OneApiRoot $oneApiRoot -MklRoot $mklRoot

# Layer 3: Clean-room smoke test (final acceptance)
Write-Host "=== Layer 3: Clean-room smoke test ===" -ForegroundColor Cyan
# Calculate resources directory from script directory (go up 3 levels: scripts -> oneapi -> windows -> quantum_espresso)
$quantumEspressoRoot = Resolve-Path (Join-Path $ScriptDir "..\..\..") -ErrorAction Stop
$resourcesDir = Join-Path $quantumEspressoRoot "tests\resources"
# Run smoke test in strict mode by default (require exit code 0 AND "JOB DONE")
$strictMode = -not $NoStrict
Run-SmokeTest -DistDir $DistRoot -ResourcesDir $resourcesDir -Strict:$strictMode

Write-Host ""
Write-Host "Staging complete. Files in: $DistRoot" -ForegroundColor Green

