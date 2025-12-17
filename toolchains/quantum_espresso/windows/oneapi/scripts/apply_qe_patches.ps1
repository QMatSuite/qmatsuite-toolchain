param(
    [string]$QeDir   = "upstream/qe",
    [string]$PatchDir = "patches"
)

$ErrorActionPreference = 'Stop'

function Resolve-RepoPath {
    param([string]$Path)
    $resolved = Resolve-Path $Path -ErrorAction SilentlyContinue
    if (-not $resolved) { throw "Path not found: $Path" }
    return $resolved.Path
}

$QePath = Resolve-RepoPath $QeDir
$PatchPath = Resolve-RepoPath $PatchDir

if (-not (Test-Path (Join-Path $QePath ".git"))) {
    throw "QE directory '$QePath' is not a git repository (.git missing). Clone QE first."
}

$patches = @(
    "qe-win-cmake-generation.patch",
    "qe-win-c-portability.patch"
)

function Apply-OnePatch {
    param([string]$PatchFile)
    $full = Join-Path $PatchPath $PatchFile
    if (-not (Test-Path $full)) { throw "Patch not found: $full" }

    Write-Host "Checking patch: $PatchFile"
    & git -C $QePath apply --check $full
    $rc = $LASTEXITCODE
    if ($rc -eq 0) {
        Write-Host "Applying patch: $PatchFile"
        & git -C $QePath apply $full
        if ($LASTEXITCODE -ne 0) { throw "Failed to apply $PatchFile (exit $LASTEXITCODE)" }
        return
    }

    # Check if already applied
    & git -C $QePath apply --reverse --check $full
    if ($LASTEXITCODE -eq 0) {
        Write-Host "Skipping (already applied): $PatchFile"
        return
    }

    throw "Patch $PatchFile failed --check and is not already applied. Resolve manually."
}

foreach ($p in $patches) {
    Apply-OnePatch -PatchFile $p
}

# Replace devxlib timer.c directly (avoids patch application issues on CI)
$srcTimer = Join-Path $PatchPath "devxlib-timer.c"
$dstTimer = Join-Path $QePath "external/devxlib/src/timer.c"
if (-not (Test-Path $srcTimer)) { throw "Patched timer not found: $srcTimer" }
Write-Host "Overwriting devxlib timer.c with patched version"
Copy-Item -Path $srcTimer -Destination $dstTimer -Force

Write-Host "All patches processed successfully."

