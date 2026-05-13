@echo off
:: ============================================================
::  MalwareODE Simulator — Script de empaquetado para Windows
::  Genera un .exe standalone con PyInstaller
:: ============================================================

echo [1/3] Instalando PyInstaller...
pip install pyinstaller --quiet

echo [2/3] Generando ejecutable...
pyinstaller --onefile ^
            --windowed ^
            --name "MalwareODE_Simulator" ^
            --add-data "src;src" ^
            src\malware_ode_simulator.py

echo [3/3] Listo!
echo El ejecutable se encuentra en: dist\MalwareODE_Simulator.exe
pause
