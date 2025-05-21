#!/usr/bin/env python3
"""
Vercel Build Script for CryoProtect v2

This script optimizes the build process for Vercel deployment by:
1. Consolidating routes to reduce serverless function count
2. Pre-generating static assets where possible
3. Setting up proper dependencies
"""

import os
import sys
import shutil
import subprocess
import json
from pathlib import Path

def main():
    """Main build function for Vercel deployment."""
    print("Starting Vercel build optimization...")
    
    # Define paths
    project_root = Path(__file__).parent.resolve()
    static_dir = project_root / "static"
    templates_dir = project_root / "templates"
    
    # Ensure directories exist
    static_dir.mkdir(exist_ok=True)
    (static_dir / "build").mkdir(exist_ok=True)
    
    # Copy static files
    copy_static_files(project_root)
    
    # Pre-generate static pages where possible
    pregenerate_static_pages(templates_dir, static_dir)
    
    # Create .vercel directory if it doesn't exist
    vercel_dir = project_root / ".vercel"
    vercel_dir.mkdir(exist_ok=True)
    
    # Create output directory
    output_dir = vercel_dir / "output"
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Copy static files to output directory
    copy_dir_contents(static_dir, output_dir / "static")
    
    # Create Vercel project settings file if it doesn't exist
    create_vercel_project_settings(project_root)
    
    print("Vercel build optimization completed successfully!")
    return 0

def copy_static_files(project_root):
    """Copy static files to the build directory."""
    source_dirs = [
        project_root / "static/js",
        project_root / "static/css",
        project_root / "static/img"
    ]
    
    target_dir = project_root / "static/build"
    target_dir.mkdir(exist_ok=True)
    
    for source_dir in source_dirs:
        if source_dir.exists():
            target_subdir = target_dir / source_dir.name
            target_subdir.mkdir(exist_ok=True)
            
            for file in source_dir.glob("*"):
                if file.is_file():
                    shutil.copy2(file, target_subdir)

def pregenerate_static_pages(templates_dir, static_dir):
    """Pre-generate static HTML for pages that don't require dynamic data."""
    static_pages = [
        "index.html",
        "login.html",
        "register.html",
        "error.html"
    ]
    
    for page in static_pages:
        template_file = templates_dir / page
        if template_file.exists():
            # For simple static generation, just copy the template
            # In a real scenario, you'd render the template with its context
            output_file = static_dir / page
            shutil.copy2(template_file, output_file)
            print(f"Pre-generated static file: {output_file}")

def copy_dir_contents(source_dir, target_dir):
    """Copy directory contents recursively."""
    if not source_dir.exists():
        return
    
    target_dir.mkdir(exist_ok=True)
    
    for item in source_dir.glob("*"):
        if item.is_file():
            shutil.copy2(item, target_dir / item.name)
        elif item.is_dir():
            copy_dir_contents(item, target_dir / item.name)

def create_vercel_project_settings(project_root):
    """Create Vercel project settings file."""
    vercel_dir = project_root / ".vercel"
    vercel_dir.mkdir(exist_ok=True)
    
    project_json = vercel_dir / "project.json"
    if not project_json.exists():
        project_data = {
            "projectId": "prj_cryoprotect_v2",
            "orgId": "your-org-id"
        }
        
        with open(project_json, 'w') as f:
            json.dump(project_data, f, indent=2)

if __name__ == "__main__":
    sys.exit(main())