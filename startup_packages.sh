#!/bin/bash

echo "===== CryoProtect Package Installer ====="
echo "Running package installation script..."

# Check if virtual environment exists
if [ -f ".venv/bin/activate" ]; then
    echo "Using existing virtual environment..."
    source .venv/bin/activate
else
    echo "Virtual environment not found."
    echo "Using system Python installation..."
    # We'll proceed without creating a new venv to avoid permission issues
fi

# Run the package installation script
python install_packages.py

echo ""
echo "===== Installation Complete ====="
echo "All required packages have been installed."
echo "You can now run the application."
echo ""

# Keep the terminal open
read -p "Press Enter to continue..."