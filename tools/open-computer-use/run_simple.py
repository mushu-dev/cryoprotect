#!/usr/bin/env python3
"""
Simple script to run the Open Computer Use tool without browser dependencies.
This script focuses on getting a VNC URL that can be opened manually.
"""

import os
import sys
import asyncio
from dotenv import load_dotenv
from os_computer_use.streaming import Sandbox

# Load environment variables
load_dotenv()

# Configure E2B
os.environ["E2B_API_KEY"] = os.getenv("E2B_API_KEY")

# Command to run (default or from command line)
command = sys.argv[1] if len(sys.argv) > 1 else "firefox"

async def main():
    print("Starting the E2B Desktop sandbox...")
    
    # Create the sandbox
    sandbox = Sandbox()
    
    try:
        # Start the VNC stream
        print("Starting the VNC server...")
        sandbox.stream.start()
        
        # Get and display the VNC URL
        vnc_url = sandbox.stream.get_url()
        print(f"VNC URL: {vnc_url}")
        print("Copy and paste this URL into your browser to view the desktop")
        print("")
        print(f"Running command: {command}")
        
        # Run the command in the sandbox
        result = sandbox.commands.run(command)
        print(f"Command result: {result.stdout}")
        
        # Keep the sandbox alive for 2 minutes
        print("\nSandbox will remain active for 2 minutes.")
        print("Press Ctrl+C to exit sooner.")
        await asyncio.sleep(120)
        
    except asyncio.CancelledError:
        print("\nOperation cancelled by user")
    except Exception as e:
        print(f"\nError: {str(e)}")
    finally:
        # Clean up
        print("Stopping the sandbox...")
        sandbox.kill()
        print("Sandbox stopped.")

if __name__ == "__main__":
    try:
        asyncio.run(main())
    except KeyboardInterrupt:
        print("\nOperation cancelled by user")
        sys.exit(0)