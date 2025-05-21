#!/bin/bash

# Headless script to run the Open Computer Use tool in a simplified way
# This script will run the tool without using the browser component
# It will still provide a VNC URL that can be opened in a browser

cd "$(dirname "$0")"

# Default prompt if none provided
PROMPT="${1:-"Open firefox and search for Claude AI"}"

# Clear the terminal
clear

echo "===================================================="
echo "Starting Open Computer Use in headless mode"
echo "===================================================="

# Print instructions
echo "When you see a VNC URL, copy and paste it into your browser"
echo "to view the virtual computer controlled by AI."
echo ""
echo "Browser component is disabled, but you will still be able to"
echo "view the virtual computer by opening the VNC URL manually."
echo ""
echo "The AI will automatically try to complete the task: $PROMPT"
echo ""
echo "Press Ctrl+C to exit when done."
echo "===================================================="
echo ""

# Run Python directly (not through Poetry) to avoid terminal interaction issues
cd "$(dirname "$0")"
python -c "
import os
import asyncio
from main import start, initialize_output_directory
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Configure proper values
output_dir = initialize_output_directory(lambda id: f'./output/run_{id}')

# Create event loop
loop = asyncio.new_event_loop()
asyncio.set_event_loop(loop)

# Run with the prompt
prompt = '$PROMPT'
loop.run_until_complete(start(user_input=prompt, output_dir=output_dir))
"

echo ""
echo "===================================================="
echo "Open Computer Use has been stopped."
echo "===================================================="