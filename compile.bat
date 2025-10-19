@echo off
setlocal

echo Building xar_g.exe with gfortran...
gfortran -O0 -Wall -Werror=unused-parameter -Werror=unused-variable -Werror=unused-function -Wno-maybe-uninitialized -Wno-surprising -fbounds-check -static -g -fmodule-private m.f90 xar.f90 -o xar_g.exe
if errorlevel 1 goto end

echo Building xar_ifx.exe with ifx...
ifx /nologo /traceback /check:bounds /check:uninit /warn:all /warn:unused /gen-interfaces /warn:interfaces /F512000000 m.f90 xar.f90 /exe:xar_ifx.exe

:end
endlocal
