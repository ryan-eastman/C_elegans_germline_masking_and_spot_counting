@echo off
REM ============================================================================
REM  quantify.bat - count RAD-51 spots in ONE .nd2 image. No coding needed.
REM
REM  HOW TO USE (two ways):
REM    1) Drag your .nd2 file ONTO this file's icon, OR
REM    2) Double-click this file, then drag your .nd2 into the black window and
REM       press Enter.
REM  Your results appear in a new folder next to your image (named *_results).
REM ============================================================================
setlocal
cd /d "%~dp0"
set "PY=.venv\Scripts\python.exe"

if not exist "%PY%" (
  echo.
  echo ERROR: the software environment is not set up ^(the .venv folder is missing^).
  echo Ask whoever set up this computer, or see the README section
  echo   "Setting up a new computer (one time)".
  echo.
  pause
  exit /b 1
)

set "IMG=%~1"
if "%IMG%"=="" set /p IMG="Drag your .nd2 image into this window, then press Enter: "
set IMG=%IMG:"=%
if not exist "%IMG%" (
  echo.
  echo ERROR: that file was not found:
  echo   %IMG%
  echo Make sure you dragged a real .nd2 file in.
  echo.
  pause
  exit /b 1
)

for %%F in ("%IMG%") do set "OUT=%%~dpnF_results"
echo.
echo ============================================================
echo  Quantifying RAD-51 spots
echo  Image:    %IMG%
echo  Results:  %OUT%
echo  This takes about 15-25 minutes. Leave this window open.
echo ============================================================
echo.

"%PY%" -m germquant.cli run "%IMG%" --config "config\config.yaml" --out "%OUT%"

echo.
if errorlevel 1 (
  echo *** Something went wrong - read the messages above, then check the
  echo     README section "If something breaks". ***
) else (
  echo *** DONE! Your results are in: ***
  echo     %OUT%
  echo.
  echo Open the file that ends in "__nuclei.csv" in Excel.
  echo The "n_spots" column is the RAD-51 count for each nucleus.
)
echo.
pause
