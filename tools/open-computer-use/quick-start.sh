#!/bin/bash

# Quick start script for Open Computer Use
# This script provides a simplified way to access the E2B Desktop sandbox

cd "$(dirname "$0")"

# Clear the terminal
clear

echo "======================================================"
echo "                Open Computer Use"
echo "======================================================"
echo ""
echo "This tool lets you interact with a virtual Linux desktop"
echo "controlled by AI. You can use it to:"
echo ""
echo "  1. Test AI-controlled browser automation"
echo "  2. Try out Linux commands in a safe sandbox"
echo "  3. Experiment with desktop interactions"
echo ""
echo "The tool will open a sandbox and give you a URL."
echo "Copy and paste that URL into your browser to see the desktop."
echo ""
echo "======================================================"
echo ""
echo "Starting the sandbox..."

# Run the simple Python script
python run_simple.py "echo 'Welcome to the sandbox! Try using the terminal.'"

# End message
echo ""
echo "======================================================"
echo "If you want to launch the sandbox again, run:"
echo "cd $(pwd) && ./quick-start.sh"
echo "======================================================"