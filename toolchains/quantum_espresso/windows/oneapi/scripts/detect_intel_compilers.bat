@echo off
REM Script to detect Intel compilers (ifort, ifx, icx) after sourcing setvars.bat
REM Usage: detect_intel_compilers.bat

setlocal enabledelayedexpansion

set "SETVARS_PATH=C:\Program Files (x86)\Intel\oneAPI\setvars.bat"

echo Checking for Intel oneAPI compilers...
echo Setvars path: %SETVARS_PATH%

REM Check if setvars.bat exists
if not exist "!SETVARS_PATH!" (
    echo ERROR: setvars.bat not found at: !SETVARS_PATH!
    exit /b 1
)

REM Source setvars.bat
call "!SETVARS_PATH!" >nul 2>&1
if errorlevel 1 (
    echo WARNING: Failed to source setvars.bat, but continuing...
)

echo.
echo === Intel Compiler Detection Results ===
echo.

set FOUND_COUNT=0

REM Check for ifort
where ifort >nul 2>&1
if !ERRORLEVEL! == 0 (
    echo [OK] ifort: FOUND
    for /f "delims=" %%i in ('where ifort 2^>nul') do (
        set IFORT_PATH=%%i
        echo     Path: %%i
    )
    echo     Checking version...
    ifort --version >nul 2>&1
    set /a FOUND_COUNT=!FOUND_COUNT!+1
) else (
    echo [FAIL] ifort: NOT FOUND
)

echo.

REM Check for ifx
where ifx >nul 2>&1
if !ERRORLEVEL! == 0 (
    echo [OK] ifx: FOUND
    for /f "delims=" %%i in ('where ifx 2^>nul') do (
        set IFX_PATH=%%i
        echo     Path: %%i
    )
    echo     Checking version...
    ifx --version >nul 2>&1
    set /a FOUND_COUNT=!FOUND_COUNT!+1
) else (
    echo [FAIL] ifx: NOT FOUND
)

echo.

REM Check for icx
where icx >nul 2>&1
if !ERRORLEVEL! == 0 (
    echo [OK] icx: FOUND
    for /f "delims=" %%i in ('where icx 2^>nul') do (
        set ICX_PATH=%%i
        echo     Path: %%i
    )
    echo     Checking version...
    icx --version >nul 2>&1
    set /a FOUND_COUNT=!FOUND_COUNT!+1
) else (
    echo [FAIL] icx: NOT FOUND
)

echo.

if !FOUND_COUNT! EQU 0 (
    echo No Intel compilers found. Make sure Intel oneAPI is properly installed.
    exit /b 1
) else (
    call echo Found !FOUND_COUNT! compiler^(s^) available.
    exit /b 0
)

