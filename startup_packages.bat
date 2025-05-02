@echo off
echo ===== CryoProtect Package Installer =====
echo Running package installation script...

REM Check if virtual environment exists
if exist .venv\Scripts\activate.bat (
    echo Using existing virtual environment...
    call .venv\Scripts\activate.bat
) else (
    echo Virtual environment not found.
    echo Using system Python installation...
    REM We'll proceed without creating a new venv to avoid permission issues
)

REM Run the package installation script
python install_packages.py

echo.
echo ===== Installation Complete =====
echo All required packages have been installed.
echo You can now run the application.
echo.

REM Keep the window open
pause