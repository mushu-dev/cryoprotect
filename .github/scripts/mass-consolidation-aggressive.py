#!/usr/bin/env python3
"""
Aggressive Mass Issue Consolidation Script

This script aggressively consolidates GitHub issues by:
1. Categorizing all issues into strategic epics
2. Creating new well-structured epic issues
3. Extracting key information from existing issues
4. Closing almost all existing issues except for a very small set

Goal: Reduce issue count from 240+ to under 50

Usage:
    python mass-consolidation-aggressive.py [--execute] [--quiet]
"""

import argparse
import json
import os
import re
import subprocess
import sys
import time
from collections import defaultdict
from datetime import datetime
from typing import Dict, List, Set, Tuple, Optional

# Epic categories to consolidate issues into
EPIC_CATEGORIES = {
    "database": {
        "title": "Database Implementation and Optimization",
        "description": "Database schema, migrations, optimization, and data population tasks",
        "patterns": ["database", "schema", "migration", "sql", "table", "supabase", "populate", "data"]
    },
    "chembl_integration": {
        "title": "ChEMBL Integration and Data Pipeline",
        "description": "ChEMBL data import, processing, and integration tasks",
        "patterns": ["chembl", "import", "pipeline", "data", "integration"]
    },
    "pubchem": {
        "title": "PubChem Integration and Molecule Management",
        "description": "PubChem data import, molecule handling, and property management",
        "patterns": ["pubchem", "molecule", "property", "cid", "chemical"]
    },
    "api": {
        "title": "API Development and Standardization",
        "description": "API endpoints, standardization, documentation, and testing",
        "patterns": ["api", "endpoint", "route", "resource", "rest", "flask"]
    },
    "frontend": {
        "title": "Frontend Implementation and User Experience",
        "description": "Frontend features, UI components, and user experience improvements",
        "patterns": ["frontend", "ui", "ux", "react", "component", "page", "view"]
    },
    "authentication": {
        "title": "Authentication and Security Implementation",
        "description": "User authentication, authorization, RLS policies, and security",
        "patterns": ["auth", "jwt", "security", "rls", "policy", "role", "permission"]
    },
    "rdkit": {
        "title": "RDKit Integration and Chemical Functionality",
        "description": "RDKit features, chemical property calculations, and molecular analysis",
        "patterns": ["rdkit", "chemical", "property", "calculation", "molecular"]
    },
    "infrastructure": {
        "title": "Infrastructure and Deployment",
        "description": "Containerization, deployment, environment setup, and DevOps",
        "patterns": ["container", "docker", "podman", "fedora", "deployment", "environment"]
    },
    "testing": {
        "title": "Testing and Quality Assurance",
        "description": "Test frameworks, validation procedures, and quality assurance processes",
        "patterns": ["test", "validation", "quality", "pytest", "coverage"]
    }
}

# Map pattern words to category keys for matching
PATTERN_TO_CATEGORY = {}
for cat, info in EPIC_CATEGORIES.items():
    for pattern in info["patterns"]:
        PATTERN_TO_CATEGORY[pattern] = cat

# Issues to absolutely keep (super high priority)
ESSENTIAL_ISSUE_KEYWORDS = [
    "roadmap", 
    "project plan", 
    "milestone", 
    "phase 5",
    "critical bug",
    "urgent security",
    "major release"
]

def run_gh_command(cmd: List[str], capture_output: bool = True) -> str:
    """Run a GitHub CLI command and return the output."""
    try:
        result = subprocess.run(
            ["gh"] + cmd,
            capture_output=capture_output,
            text=True,
            check=True
        )
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error running GitHub CLI command: {e}")
        print(f"Error output: {e.stderr}")
        sys.exit(1)

def get_all_issues() -> List[Dict]:
    """Get all open issues."""
    try:
        return json.loads(run_gh_command([
            "issue", 
            "list", 
            "--json", "number,title,body,createdAt,updatedAt,author,labels,comments",
            "--state", "open",
            "--limit", "500"
        ]))
    except json.JSONDecodeError:
        print("Error parsing JSON from GitHub CLI output.")
        sys.exit(1)

def create_epic_issue(title: str, body: str, labels: List[str]) -> int:
    """Create a new epic issue and return its number."""
    # Create separate label commands for each label
    label_args = []
    for label in labels:
        label_args.extend(["--label", label])
    
    issue_json = run_gh_command([
        "issue", "create",
        "--title", title,
        "--body", body,
        *label_args
    ])
    try:
        match = re.search(r'#(\d+)', issue_json)
        return int(match.group(1)) if match else 0
    except:
        print(f"Could not parse issue number from output: {issue_json}")
        return 0

def close_issue(number: int, comment: str) -> None:
    """Close an issue with a comment."""
    run_gh_command([
        "issue", "close", str(number),
        "--comment", comment
    ], capture_output=False)

def is_essential_issue(issue: Dict) -> bool:
    """Determine if an issue is absolutely essential and must be kept."""
    content = (issue.get("title", "") + " " + issue.get("body", "")).lower()
    
    # Check for essential keywords
    for keyword in ESSENTIAL_ISSUE_KEYWORDS:
        if keyword.lower() in content:
            return True
    
    # Check for very recent activity (last 2 days)
    updated_at = issue.get("updatedAt", "")
    if updated_at:
        try:
            updated_time = datetime.fromisoformat(updated_at.replace("Z", "+00:00"))
            now = datetime.now().astimezone()
            days_ago = (now - updated_time).days
            if days_ago <= 2:  # Very recent activity
                return True
        except:
            pass
    
    # Check for "phase" in the title
    title = issue.get("title", "").lower()
    if re.search(r'phase\s+[0-9]', title):
        return True
    
    # Check for high-priority labels
    high_priority = False
    for label in issue.get("labels", []):
        label_name = label.get("name", "").lower()
        if "priority:high" in label_name:
            high_priority = True
            
            # Only keep high priority issues created in last 30 days
            created_at = issue.get("createdAt", "")
            if created_at:
                try:
                    created_time = datetime.fromisoformat(created_at.replace("Z", "+00:00"))
                    now = datetime.now().astimezone()
                    days_ago = (now - created_time).days
                    if days_ago <= 30:  # Created in last month
                        return True
                except:
                    pass
    
    return False

def categorize_issue(issue: Dict) -> List[str]:
    """Categorize an issue based on its content and return matching categories."""
    categories = set()
    content = (issue.get("title", "") + " " + issue.get("body", "")).lower()
    
    # Check for explicit area labels
    for label in issue.get("labels", []):
        label_name = label.get("name", "").lower()
        if label_name.startswith("area:"):
            area = label_name[5:]
            for cat in EPIC_CATEGORIES:
                if area in EPIC_CATEGORIES[cat]["patterns"]:
                    categories.add(cat)
    
    # Check content against patterns
    for pattern in PATTERN_TO_CATEGORY:
        if pattern in content:
            categories.add(PATTERN_TO_CATEGORY[pattern])
    
    return list(categories) if categories else ["uncategorized"]

def extract_key_info(issue: Dict) -> str:
    """Extract the most important information from an issue."""
    summary = []
    title = issue.get("title", "")
    body = issue.get("body", "")
    
    # Add issue info
    summary.append(f"### From #{issue['number']}: {title}")
    
    # Extract key sections from body
    sections = re.split(r'#{1,3} ', body)
    for section in sections:
        if not section.strip():
            continue
        
        # Try to identify important sections
        section_title = section.split("\n")[0].strip().lower()
        section_content = "\n".join(section.split("\n")[1:]).strip()
        
        if any(key in section_title for key in [
            "description", "use case", "problem", "solution", 
            "approach", "implementation", "acceptance"
        ]):
            # Include first 250 chars of important sections
            if section_content:
                summary.append(f"**{section_title.title()}**: {section_content[:250]}...")
    
    # If no structured content was found, include first 200 chars of body
    if len(summary) < 2 and body:
        summary.append(body[:200] + "..." if len(body) > 200 else body)
    
    return "\n\n".join(summary)

def main():
    parser = argparse.ArgumentParser(description="Mass consolidate GitHub issues")
    parser.add_argument("--execute", action="store_true", help="Execute the consolidation (without this, runs in dry-run mode)")
    parser.add_argument("--quiet", action="store_true", help="Run with minimal output")
    args = parser.parse_args()
    
    print("Fetching all open issues...")
    all_issues = get_all_issues()
    print(f"Found {len(all_issues)} open issues.")
    
    if not args.quiet:
        print("\nAnalyzing issues...")
    
    # Step 1: Categorize all issues
    categorized_issues = defaultdict(list)
    essential_issues = []
    
    for issue in all_issues:
        # Check if it's an essential issue first
        if is_essential_issue(issue):
            essential_issues.append(issue)
            continue
        
        # Categorize remaining issues
        categories = categorize_issue(issue)
        categorized_issues[categories[0]].append(issue)
    
    # Print summary
    print("\n=== CATEGORIZATION SUMMARY ===")
    print(f"Essential issues to keep: {len(essential_issues)}")
    print("\nIssues by category:")
    for cat in sorted(categorized_issues.keys()):
        print(f"  {cat}: {len(categorized_issues[cat])}")
    
    # Calculate how many issues would remain after consolidation
    total_epics = len(EPIC_CATEGORIES)
    total_essential = len(essential_issues)
    consolidated_count = total_epics + total_essential
    
    print(f"\nAfter consolidation: ~{consolidated_count} issues")
    print(f"  {total_epics} epic issues")
    print(f"  {total_essential} essential issues")
    
    if not args.execute:
        print("\nDRY RUN: No changes have been made.")
        print("Run with --execute to perform the consolidation.")
        return
    
    print("\n=== EXECUTING CONSOLIDATION ===")
    
    # Step 2: Create epic issues
    epic_numbers = {}
    for cat, info in EPIC_CATEGORIES.items():
        print(f"Creating epic issue for {cat}...")
        body = f"# {info['title']}\n\n{info['description']}\n\n## Consolidated Issues\n\n"
        
        # Add related issues to the body
        if cat in categorized_issues:
            for issue in categorized_issues[cat]:
                extracted_info = extract_key_info(issue)
                body += f"{extracted_info}\n\n---\n\n"
        
        labels = ["type:epic", "status:planning"]
        if cat in ["database", "chembl_integration", "pubchem"]:
            labels.append("area:database")
        elif cat in ["api"]:
            labels.append("area:api")
        elif cat in ["frontend"]:
            labels.append("area:ui")
        elif cat in ["authentication"]:
            labels.append("area:auth")
        elif cat in ["testing"]:
            labels.append("area:testing")
        elif cat in ["rdkit"]:
            labels.append("area:chembl")
        elif cat in ["infrastructure"]:
            labels.append("area:devops")
        
        # Add priority
        labels.append("priority:high")
        
        epic_num = create_epic_issue(info["title"], body, labels)
        if epic_num:
            epic_numbers[cat] = epic_num
            print(f"Created epic issue #{epic_num}")
            
            # Avoid rate limiting
            time.sleep(2)
    
    # Step 3: Close issues that have been consolidated (almost all issues)
    closed_count = 0
    
    for cat, issues in categorized_issues.items():
        if cat not in epic_numbers:
            continue
            
        epic_num = epic_numbers[cat]
        for issue in issues:
            # Skip essential issues (should be empty intersection anyway)
            if any(e["number"] == issue["number"] for e in essential_issues):
                continue
                
            issue_num = issue["number"]
            print(f"Closing issue #{issue_num} (consolidated into #{epic_num})...")
            close_issue(
                issue_num, 
                f"This issue has been consolidated into epic #{epic_num} as part of repository cleanup and organization.\n\n"
                f"Please refer to the epic for tracking this work going forward."
            )
            closed_count += 1
            time.sleep(1)
    
    # Final report
    open_issues = len(all_issues) - closed_count
    print(f"\n=== CONSOLIDATION COMPLETE ===")
    print(f"Created {len(epic_numbers)} epic issues")
    print(f"Closed {closed_count} issues")
    print(f"Remaining open issues: {open_issues}")
    
    print("\nRecommended next steps:")
    print("1. Review the created epic issues and ensure all key information was captured")
    print("2. Review any remaining open issues to ensure they're still relevant")
    print("3. Use the project board to organize the remaining issues")

if __name__ == "__main__":
    main()