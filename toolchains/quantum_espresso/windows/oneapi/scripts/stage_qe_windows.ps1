<#
.SYNOPSIS
Stages Windows oneAPI QE executables plus runtime DLLs into dist/win-oneapi.

.PARAMETER BuildDir
QE build dir containing bin (default: ../upstream/qe/build-win-oneapi)

.PARAMETER DistDir
Destination dist directory (default: dist/win-oneapi under oneapi root)

.PARAMETER Only
Optional exe basenames to stage.

.PARAMETER Clean
Wipe dist dir before staging.

.PARAMETER VerboseDeps
Show parsed dependencies when using dumpbin.
#>
[CmdletBinding(SupportsShouldProcess=$true)]
param(
    [string]$BuildDir = $null,
    [string]$DistDir  = $null,
    [string[]]$Only,
    [switch]$Clean,
    [switch]$VerboseDeps
)

$ErrorActionPreference = 'Stop'
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$OneApiRoot = Split-Path $ScriptDir -Parent

if (-not $BuildDir) { $BuildDir = Join-Path $OneApiRoot "upstream/qe/build-win-oneapi" }
if (-not $DistDir)  { $DistDir  = Join-Path $OneApiRoot "dist/win-oneapi" }

function Resolve-PathSafe {
    param([string]$Path)
    if ([System.IO.Path]::IsPathRooted($Path)) {
        if (Test-Path $Path) { return (Resolve-Path $Path -ErrorAction Stop).Path }
        return $Path
    } else {
        $joined = Join-Path $OneApiRoot $Path
        if (Test-Path $joined) { return (Resolve-Path $joined -ErrorAction Stop).Path }
        return $joined
    }
}

$BuildDir = Resolve-PathSafe $BuildDir
$DistDir  = Resolve-PathSafe $DistDir

Write-Host "QE build bin : $BuildDir\\bin"
Write-Host "Dist target  : $DistDir"

$BinDir = Join-Path $BuildDir "bin"
if (-not (Test-Path $BinDir)) { throw "Bin directory not found: $BinDir" }

if ($Clean -and (Test-Path $DistDir)) {
    if ($PSCmdlet.ShouldProcess($DistDir,"Clean dist")) { Remove-Item -Recurse -Force $DistDir }
}
if (-not (Test-Path $DistDir)) { New-Item -ItemType Directory -Force -Path $DistDir | Out-Null }

# dumpbin detection
function Get-Dumpbin {
    if ($env:VCToolsInstallDir) {
        $cand = Join-Path $env:VCToolsInstallDir "bin\\Hostx64\\x64\\dumpbin.exe"
        if (Test-Path $cand) { return $cand }
    }
    $cmd = Get-Command dumpbin.exe -ErrorAction SilentlyContinue
    if ($cmd) { return $cmd.Source }
    return $null
}
$DumpbinPath = Get-Dumpbin

# stage exes
$exeFiles = Get-ChildItem -Path $BinDir -Filter *.exe -File
if ($Only) {
    $norm = $Only | ForEach-Object { if ($_ -like "*.exe") { $_ } else { "$_.exe" } }
    $exeFiles = $exeFiles | Where-Object { $norm -contains $_.Name }
}
if (-not $exeFiles) { throw "No executables to stage (check build or -Only filters)." }

Write-Host "Copying executables..."
foreach ($exe in $exeFiles) {
    $dest = Join-Path $DistDir $exe.Name
    if ($PSCmdlet.ShouldProcess($dest,"Copy exe")) { Copy-Item -Force $exe.FullName $dest }
    Write-Host "  exe: $($exe.Name)"
}
$stagedExePaths = $exeFiles | ForEach-Object { Join-Path $DistDir $_.Name }

# search roots for DLLs
$searchRoots = @()
if ($env:ONEAPI_ROOT) { $oneApi = $env:ONEAPI_ROOT } else { $oneApi = "C:\Program Files (x86)\Intel\oneAPI" }
if (Test-Path $oneApi) {
    $dirs = Get-ChildItem -Path $oneApi -Directory -ErrorAction SilentlyContinue
    foreach ($d in $dirs) {
        $searchRoots += @(
            (Join-Path $d.FullName "bin"),
            (Join-Path $d.FullName "lib"),
            (Join-Path $d.FullName "redist")
        )
    }
    $compilerBase = Join-Path $oneApi "compiler"
    if (Test-Path $compilerBase) {
        foreach ($d in Get-ChildItem $compilerBase -Directory -ErrorAction SilentlyContinue) {
            $searchRoots += @(
                (Join-Path $d.FullName "bin"),
                (Join-Path $d.FullName "lib"),
                (Join-Path $d.FullName "redist\\intel64\\compiler"),
                (Join-Path $d.FullName "redist\\intel64"),
                (Join-Path $d.FullName "windows\\redist\\intel64")
            )
        }
    }
    $mklBase = Join-Path $oneApi "mkl"
    if (Test-Path $mklBase) {
        foreach ($d in Get-ChildItem $mklBase -Directory -ErrorAction SilentlyContinue) {
            $searchRoots += @(
                (Join-Path $d.FullName "bin"),
                (Join-Path $d.FullName "lib"),
                (Join-Path $d.FullName "redist\\intel64"),
                (Join-Path $d.FullName "redist\\intel64\\mkl")
            )
        }
    }
}
$searchRoots = $searchRoots | Where-Object { $_ -and (Test-Path $_) } | Select-Object -Unique
$pathDirs = ($env:PATH -split ';') | Where-Object { $_ -and (Test-Path $_) }

function Is-SystemDll {
    param($Name,$Resolved)
    $sysRoots = @("$env:WINDIR\\System32","$env:WINDIR\\SysWOW64")
    foreach ($sr in $sysRoots) { if ($Resolved -like "$sr*") { return $true } }
    $prefix = @("kernel32","user32","advapi32","gdi32","shell32","ole32","oleaut32","ws2_32","secur32","bcrypt","ntdll","ucrtbase","msvcrt","vcruntime")
    foreach ($p in $prefix) { if ($Name.ToLower().StartsWith($p)) { return $true } }
    return $false
}

function Resolve-DllPath {
    param($Name)
    foreach ($r in $searchRoots) {
        $f = Get-ChildItem -Path $r -Filter $Name -File -Recurse -ErrorAction SilentlyContinue | Select-Object -First 1
        if ($f) { return $f.FullName }
    }
    foreach ($d in $pathDirs) {
        $cand = Join-Path $d $Name
        if (Test-Path $cand) { return (Resolve-Path $cand).Path }
    }
    return $null
}

function Get-Dependents {
    param($Path)
    $deps=@()
    if (-not $DumpbinPath) { return $deps }
    $out = & $DumpbinPath /DEPENDENTS "$Path" 2>$null
    $capture=$false
    foreach ($line in $out) {
        if ($line -match "Image has the following dependencies") { $capture=$true; continue }
        if ($capture -and [string]::IsNullOrWhiteSpace($line)) { $capture=$false }
        if ($capture -and $line -match "([A-Za-z0-9._-]+\\.dll)") {
            $deps += $matches[1].ToLower()
        }
    }
    return ($deps | Select-Object -Unique)
}

function Copy-DllClosure {
    param([string[]]$Roots)
    if (-not $DumpbinPath) { return }
    $queue = New-Object System.Collections.Queue
    $seen = @{}
    foreach ($exe in $stagedExePaths) {
        foreach ($d in Get-Dependents -Path $exe) { $queue.Enqueue($d) }
    }
    while ($queue.Count -gt 0) {
        $name = $queue.Dequeue()
        if ($seen.ContainsKey($name)) { continue }
        $seen[$name] = $true
        $resolved = Resolve-DllPath $name
        if (-not $resolved) { Write-Warning "Missing DLL: $name"; continue }
        if (Is-SystemDll -Name $name -Resolved $resolved) { continue }
        $dest = Join-Path $DistDir (Split-Path $resolved -Leaf)
        if ($PSCmdlet.ShouldProcess($dest,"Copy dll $resolved")) { Copy-Item -Force $resolved $dest }
        if ($VerboseDeps) { Write-Host "  dll: $name (from $resolved)" }
        foreach ($d in Get-Dependents -Path $resolved) { $queue.Enqueue($d) }
    }
}

if ($DumpbinPath) {
    Write-Host "Using dumpbin at: $DumpbinPath"
    Copy-DllClosure -Roots $searchRoots
} else {
    Write-Warning "dumpbin.exe not found; dependency-based copy disabled."
}

# fallback patterns (optional)
function Copy-Patterns($patterns,$optional=$true) {
    foreach ($p in $patterns) {
        $found=$null
        foreach ($r in $searchRoots) {
            $found = Get-ChildItem -Path $r -Filter $p -File -Recurse -ErrorAction SilentlyContinue | Select-Object -First 1
            if ($found) { break }
        }
        if ($found) {
            $dest = Join-Path $DistDir $found.Name
            if ($PSCmdlet.ShouldProcess($dest,"Copy dll $($found.FullName)")) { Copy-Item -Force $found.FullName $dest }
            Write-Host "  dll: $($found.Name) (from $($found.DirectoryName))"
        } elseif (-not $optional) {
            Write-Warning "Missing DLL pattern '$p'"
        }
    }
}
Copy-Patterns @("libiomp5md.dll","mkl_core*.dll","mkl_intel_thread*.dll","mkl_rt*.dll","libifcoremd*.dll","libifportmd*.dll","libmmd*.dll","svml_dispmd*.dll","libintlc*.dll","libimalloc*.dll")

# probe
$probeTarget = ($exeFiles | Where-Object { $_.Name -ieq "pw.x.exe" } | Select-Object -First 1).FullName
if (-not $probeTarget) { $probeTarget = $exeFiles[0].FullName }
$probeExe = Join-Path $DistDir (Split-Path $probeTarget -Leaf)
Write-Host "Probing '$probeExe'..."
$psi = New-Object System.Diagnostics.ProcessStartInfo
$psi.FileName = $probeExe
$psi.Arguments = "-help"
$psi.UseShellExecute = $false
$psi.RedirectStandardError = $true
$psi.RedirectStandardOutput = $true
$proc = New-Object System.Diagnostics.Process
$proc.StartInfo = $psi
$proc.Start() | Out-Null
$out = $proc.StandardOutput.ReadToEnd()
$err = $proc.StandardError.ReadToEnd()
$proc.WaitForExit()
if ($proc.ExitCode -ne 0 -and ($err + $out) -match "(was not found|missing|cannot proceed)") {
    Write-Warning "Probe indicates missing dependency (exit $($proc.ExitCode)):`n$err"
} else {
    Write-Host "Probe finished (exit $($proc.ExitCode))."
}

Write-Host "Staging complete. Files in: $DistDir"

