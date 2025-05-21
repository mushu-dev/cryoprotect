#!/usr/bin/env python3
"""
A simple script to run the Open Computer Use tool with a prompt.
This allows running the tool from a script rather than requiring
an interactive terminal.
"""

import asyncio
import os
import sys
import importlib
from dotenv import load_dotenv
from main import start, initialize_output_directory

# Load environment variables from .env file
load_dotenv()

# Configure E2B
os.environ["E2B_API_KEY"] = os.getenv("E2B_API_KEY")

# Set default prompt if not provided
prompt = sys.argv[1] if len(sys.argv) > 1 else "Please open a web browser and search for the current weather"

# Fix config to use Anthropic provider since we have that key
def update_config():
    # Create a proper config.py file with a known-good configuration
    with open('config.py', 'w') as f:
        f.write("""# Using Anthropic for both vision and action models
from os_computer_use import providers

grounding_model = providers.OSAtlasProvider()
vision_model = providers.AnthropicProvider("claude-3.5-sonnet")
action_model = providers.AnthropicProvider("claude-3.5-sonnet")
""")
    
    # Force reload config module
    import os_computer_use.config
    importlib.reload(os_computer_use.config)

def main():
    update_config()
    output_dir = initialize_output_directory(lambda id: f"./output/run_{id}")
    loop = asyncio.new_event_loop()
    asyncio.set_event_loop(loop)
    loop.run_until_complete(start(user_input=prompt, output_dir=output_dir))

if __name__ == "__main__":
    main()