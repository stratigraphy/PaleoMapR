@echo off
title PaleoMapR setup
rem Uses the newest released R in C:\Program Files\R (R-4.x.y), not R-devel
set RS=
for /d %%D in ("C:\Program Files\R\R-?.*") do (
  if exist "%%D\bin\x64\Rscript.exe" set RS=%%D\bin\x64\Rscript.exe
)
if "%RS%"=="" (
  for /d %%D in ("%LOCALAPPDATA%\Programs\R\R-?.*") do (
    if exist "%%D\bin\x64\Rscript.exe" set RS=%%D\bin\x64\Rscript.exe
  )
)
if "%RS%"=="" (
  echo R was not found. Install R from https://cloud.r-project.org/bin/windows/base/ and run SETUP.bat again.
  pause
  exit /b 1
)
echo Using %RS%
"%RS%" "%~dp0setup.R" %*
echo.
pause
