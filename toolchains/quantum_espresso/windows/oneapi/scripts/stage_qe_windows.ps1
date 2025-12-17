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
#>
[CmdletBinding(SupportsShouldProcess = $true)]
param(
    [string]$BuildDir = "upstream/qe/build-win-oneapi",
    [string]$DistDir  = "dist/win-oneapi",
    [string[]]$Only,
    [switch]$Clean,
    [switch]$VerboseDeps,
    [string]$ProbeExe
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
        [switch]$Verbose
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
    $stdout = $proc.StandardOutput.ReadToEnd()
    $stderr = $proc.StandardError.ReadToEnd()
    $proc.WaitForExit()
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
    foreach ($exe in $exeFiles) {
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

function Copy-DllClosure {
    param(
        [string[]]$ExePaths,
        [string]$DumpbinPath,
        [string[]]$SearchRoots,
        [string]$DistDir,
        [switch]$VerboseDeps,
        [int]$MaxDepth = 5
    )
    $pathDirs = ($env:PATH -split ';') | Where-Object { $_ -and (Test-Path $_) }
    $queue = New-Object System.Collections.Queue
    $processed = @{}
    $missing = [System.Collections.Generic.List[string]]::new()

    foreach ($exe in $ExePaths) {
        $deps = Get-ExeDependentsDumpbin -Path $exe -DumpbinPath $DumpbinPath -Verbose:$VerboseDeps
        foreach ($d in $deps) { $queue.Enqueue(@($d,0)) }
    }

    while ($queue.Count -gt 0) {
        $item = $queue.Dequeue()
        $name = $item[0]
        $depth = $item[1]
        if ($processed.ContainsKey($name)) { continue }
        $processed[$name] = $true

        $resolved = Resolve-DllPath -Name $name -SearchRoots $SearchRoots -PathDirs $pathDirs
        if (-not $resolved) {
            $missing.Add($name)
            continue
        }
        if (Is-SystemDll -Name $name -ResolvedPath $resolved) { continue }

        $dest = Join-Path $DistDir (Split-Path $resolved -Leaf)
        if ($PSCmdlet.ShouldProcess($dest, "Copy dll $resolved") -and -not $WhatIfPreference) {
            Copy-Item -Force $resolved $dest
        }
        if ($VerboseDeps) { Write-Host "  dll: $name (from $resolved)" }

        if ($depth -lt $MaxDepth) {
            $dllDeps = Get-ExeDependentsDumpbin -Path $resolved -DumpbinPath $DumpbinPath -Verbose:$VerboseDeps
            foreach ($d in $dllDeps) { $queue.Enqueue(@($d, $depth + 1)) }
        }
    }

    if ($missing.Count -gt 0) {
        Write-Warning ("Missing non-system DLLs: " + ($missing | Sort-Object -Unique -join ", "))
        Write-Warning "Ensure oneAPI redistributables are installed and ONEAPI_ROOT is correct; re-run setvars if needed."
    }
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
    try {
        $psi = New-Object System.Diagnostics.ProcessStartInfo
        $psi.FileName = $ExePath
        $psi.Arguments = "-help"
        $psi.RedirectStandardOutput = $true
        $psi.RedirectStandardError = $true
        $psi.UseShellExecute = $false
        $proc = New-Object System.Diagnostics.Process
        $proc.StartInfo = $psi
        $proc.Start() | Out-Null
        $outText = $proc.StandardOutput.ReadToEnd()
        $errText = $proc.StandardError.ReadToEnd()
        $proc.WaitForExit()
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

# Dependency-based copy using dumpbin if available
$dumpbinPath = Get-DumpbinPath
if ($dumpbinPath) {
    Write-Host "Using dumpbin at: $dumpbinPath"
    $searchRoots = Get-SearchRoots -OneApiRoot $oneApiRoot
    Copy-DllClosure -ExePaths $stagedExePaths -DumpbinPath $dumpbinPath -SearchRoots $searchRoots -DistDir $DistRoot -VerboseDeps:$VerboseDeps -WhatIf:$WhatIf
} else {
    Write-Warning "dumpbin.exe not found; falling back to pattern-based DLL copy."
    $searchRoots = Get-SearchRoots -OneApiRoot $oneApiRoot
    Write-Host "Copying required DLLs (patterns)..."
    # OpenMP
    Copy-DllPatterns -Patterns @("libiomp5md.dll") -SearchRoots $searchRoots -DistDir $DistRoot -WhatIf:$WhatIf
    # MKL
    Copy-DllPatterns -Patterns @("mkl_core*.dll","mkl_intel_lp64*.dll","mkl_intel_thread*.dll","mkl_rt*.dll") -SearchRoots $searchRoots -DistDir $DistRoot -Optional -WhatIf:$WhatIf
    # Intel Fortran runtime
    Copy-DllPatterns -Patterns @("libifcoremd*.dll","libifportmd*.dll","libmmd*.dll","svml_dispmd*.dll","libintlc*.dll","libimalloc*.dll") -SearchRoots $searchRoots -DistDir $DistRoot -Optional -WhatIf:$WhatIf
}

# Optional extra pattern copy even if dumpbin succeeded (to catch loose files)
if ($dumpbinPath) {
    $searchRoots = Get-SearchRoots -OneApiRoot $oneApiRoot
    Write-Host "Optional pattern copy (supplemental)..."
    Copy-DllPatterns -Patterns @("mkl_intel_lp64*.dll","libintlc*.dll") -SearchRoots $searchRoots -DistDir $DistRoot -Optional -WhatIf:$WhatIf
}

# Probe
$probeTarget = $ProbeExe
if (-not $probeTarget -and $stagedExeFiles) {
    $probeTarget = ($stagedExeFiles | Where-Object { $_.Name -ieq "pw.x.exe" }).Name
    if (-not $probeTarget) { $probeTarget = $stagedExeFiles[0].Name }
}
if ($probeTarget) {
    $probePath = Join-Path $DistRoot $probeTarget
    Write-Host "Probing '$probePath' for missing DLLs..."
    Run-Probe -ExePath $probePath -WhatIf:$WhatIf
}

Write-Host "Staging complete. Files in: $DistRoot"

