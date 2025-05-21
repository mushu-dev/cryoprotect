#!/bin/bash

# Simple wrapper script to run the Open Computer Use tool
# with proper environment setup

cd "$(dirname "$0")"

# Default prompt if none provided
PROMPT="${1:-"Open a browser and go to GitHub"}"

# Run the Python script with poetry
poetry run python run_with_prompt.py "$PROMPT"

# Print a message about how to use the VNC URL
echo ""
echo "===================================================="
echo "If you see a VNC URL above, copy and paste it into your browser"
echo "to interact with the virtual computer controlled by AI."
echo "===================================================="