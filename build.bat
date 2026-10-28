@echo off
REM Compilation script for Windows (MSYS2 UCRT64 / MinGW GFortran)
echo Compiling modernized QUADESTH Fortran 90 solver...

if exist C:\msys64\ucrt64\bin (
    set "PATH=C:\msys64\ucrt64\bin;%PATH%"
)

gfortran -O3 -fopenmp -fdefault-real-8 -c data_mod.f90
gfortran -O3 -fopenmp -fdefault-real-8 -c solver_mod.f90
gfortran -O3 -fopenmp -fdefault-real-8 -c assembly_mod.f90
gfortran -O3 -fopenmp -fdefault-real-8 -c stress_mod.f90
gfortran -O3 -fopenmp -fdefault-real-8 -c io_mod.f90
gfortran -O3 -fopenmp -fdefault-real-8 -c main.f90

gfortran -O3 -fopenmp -fdefault-real-8 data_mod.o solver_mod.o assembly_mod.o stress_mod.o io_mod.o main.o -static -llapack -lblas -o quadesth_modern.exe

if %ERRORLEVEL% EQU 0 (
    echo Compilation successful: quadesth_modern.exe
    del /f /q *.o *.mod 2>nul
) else (
    echo Compilation failed!
)
