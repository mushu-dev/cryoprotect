#!/usr/bin/env python
"""
GitHub Issue Epic Organizer

This script helps organize GitHub issues into the 9 epic categories.
It analyzes issue content to suggest appropriate epic categories,
creates epic issues if needed, and groups related issues together.

Usage:
    python organize_by_epic.py [--dry-run] [--create-epics] [--interactive]

Options:
    --dry-run       Run in dry run mode (no changes applied)
    --create-epics  Create any missing epic issues
    --interactive   Prompt for confirmation before making changes
"""

import argparse
import json
import os
import re
import subprocess
import sys
from datetime import datetime
from typing import Dict, List, Tuple, Optional, Any

# Epic categories
EPIC_CATEGORIES = [
    "Database Implementation and Optimization",
    "ChEMBL Integration and Data Pipeline",
    "PubChem Integration and Molecule Management",
    "API Development and Standardization",
    "Frontend Implementation and User Experience",
    "Authentication and Security Implementation",
    "RDKit Integration and Chemical Functionality",
    "Infrastructure and Deployment",
    "Testing and Quality Assurance"
]

# Keywords for each epic category
CATEGORY_KEYWORDS = {
    "Database Implementation and Optimization": 
        ["database", "schema", "sql", "query", "supabase", "postgres", "migration", "integrity", 
         "db", "pool", "connection pool", "field", "index", "foreign key", "constraint"],
    
    "ChEMBL Integration and Data Pipeline": 
        ["chembl", "chem", "integration", "import", "pipeline", "etl", "drug", "compound", 
         "data import", "chembl_", "data quality", "reconcile", "identifier"],
    
    "PubChem Integration and Molecule Management": 
        ["pubchem", "molecule", "cid", "pubchem_", "property", "consolidate", "duplicate", 
         "molecule management", "standardize", "formula", "structure", "smiles", "cryoprotectant"],
    
    "API Development and Standardization": 
        ["api", "endpoint", "flask", "route", "restful", "standardize api", "openapi", "swagger", 
         "response", "request", "pagination", "api standard", "resource", "middleware", "rate limit"],
    
    "Frontend Implementation and User Experience": 
        ["frontend", "ui", "ux", "react", "next", "vercel", "component", "page", "route", "interface", 
         "navigation", "dashboard", "visualization", "style", "css", "layout", "usability", "accessibility"],
    
    "Authentication and Security Implementation": 
        ["auth", "jwt", "token", "security", "permission", "role", "rbac", "authentication", "authorization", 
         "login", "logout", "session", "credential", "access control", "service role", "rls", "policy"],
    
    "RDKit Integration and Chemical Functionality": 
        ["rdkit", "chemical", "molecule", "property calculation", "cryoprotection", "prediction", 
         "toxicity", "container", "docker", "rdkit_", "calculation", "molecular property", "formula"],
    
    "Infrastructure and Deployment": 
        ["deploy", "infrastructure", "container", "docker", "podman", "heroku", "vercel", "ci/cd", 
         "environment", "config", "production", "staging", "selinux", "systemd", "fedora", "backup"],
    
    "Testing and Quality Assurance": 
        ["test", "qa", "quality", "coverage", "unit test", "integration test", "validation", "verification", 
         "benchmark", "performance", "load test", "assertion", "mock", "fixture", "reliability"]
}

class EpicOrganizer:
    def __init__(self, dry_run=True, create_epics=False, interactive=False):
        self.dry_run = dry_run
        self.create_epics = create_epics
        self.interactive = interactive
        self.issues = []
        self.epics = {}  # Epic title -> issue number
        self.issue_epic_map = {}  # issue number -> suggested epic
        self.epic_issues = {}  # epic category -> epic issue number
        
        # Output directory
        os.makedirs("output", exist_ok=True)
        
        # Report data
        self.report = {
            "timestamp": datetime.now().isoformat(),
            "total_issues": 0,
            "epic_issues": {},
            "issues_by_epic": {},
            "unclassified_issues": []
        }
    
    def run_gh_command(self, args, check=True, get_output=True):
        """Run a GitHub CLI command"""
        cmd = ["gh"] + args
        
        if get_output:
            result = subprocess.run(cmd, check=check, capture_output=True, text=True)
            return result.stdout.strip() if result.stdout else ""
        else:
            subprocess.run(cmd, check=check)
            return None
    
    def fetch_issues(self):
        """Fetch all open issues from GitHub"""
        print("Fetching issues from GitHub...")
        issues_json = self.run_gh_command(["issue", "list", "--state", "open", "--json", 
                                          "number,title,body,labels,assignees,url"])
        self.issues = json.loads(issues_json)
        self.report["total_issues"] = len(self.issues)
        print(f"Found {len(self.issues)} open issues")
        
        # Identify existing epic issues
        for issue in self.issues:
            labels = [label["name"] for label in issue.get("labels", [])]
            if "type:epic" in labels:
                epic_name = issue["title"]
                if ":" in epic_name:
                    epic_name = epic_name.split(":", 1)[1].strip()
                self.epics[epic_name] = issue["number"]
                
                # Record in report
                if epic_name in EPIC_CATEGORIES:
                    self.epic_issues[epic_name] = issue["number"]
                    self.report["epic_issues"][epic_name] = {
                        "number": issue["number"],
                        "url": issue["url"]
                    }
        
        print(f"Found {len(self.epics)} existing epic issues")
        
        # Check which epic categories don't have epic issues
        missing_epics = []
        for category in EPIC_CATEGORIES:
            if not any(category in epic_name for epic_name in self.epics.keys()):
                missing_epics.append(category)
        
        if missing_epics:
            print(f"Missing epic issues for {len(missing_epics)} categories:")
            for category in missing_epics:
                print(f"  - {category}")
            
            if self.create_epics:
                self.create_missing_epics(missing_epics)
        
        return self.issues
    
    def create_missing_epics(self, missing_epics):
        """Create epic issues for missing categories"""
        print("\nCreating missing epic issues...")
        
        for category in missing_epics:
            if self.dry_run:
                print(f"  [DRY RUN] Would create epic issue for: {category}")
                continue
            
            if self.interactive:
                choice = input(f"Create epic issue for '{category}'? [y/N]: ").lower()
                if choice != 'y':
                    print(f"  Skipping creation of epic for {category}")
                    continue
            
            print(f"  Creating epic issue for: {category}")
            
            # Prepare epic issue content
            title = f"Epic: {category}"
            body = f"# {category} Epic\n\n"
            body += "This epic tracks all issues related to this category.\n\n"
            body += "## Overview\n\n"
            
            # Add category-specific content
            if "Database" in category:
                body += "This epic covers database schema design, query optimization, connection pooling, and data integrity.\n"
            elif "ChEMBL" in category:
                body += "This epic covers integration with ChEMBL database, data import processes, and ETL pipelines.\n"
            elif "PubChem" in category:
                body += "This epic covers integration with PubChem, molecule properties, and data consolidation.\n"
            elif "API" in category:
                body += "This epic covers API endpoints, standardization, documentation, and gateway services.\n"
            elif "Frontend" in category:
                body += "This epic covers UI components, responsive design, user workflows, and visual styling.\n"
            elif "Authentication" in category:
                body += "This epic covers JWT auth, RLS policies, RBAC, and security enhancements.\n"
            elif "RDKit" in category:
                body += "This epic covers RDKit services, chemical calculations, and molecular analysis.\n"
            elif "Infrastructure" in category:
                body += "This epic covers deployment processes, container orchestration, CI/CD, and hosting configuration.\n"
            elif "Testing" in category:
                body += "This epic covers test frameworks, CI validation, and quality assurance procedures.\n"
            
            body += "\n## Related Keywords\n\n"
            if category in CATEGORY_KEYWORDS:
                body += ", ".join(CATEGORY_KEYWORDS[category])
            
            # Create the issue
            try:
                result = self.run_gh_command([
                    "issue", "create",
                    "--title", title,
                    "--body", body,
                    "--label", "type:epic",
                    "--label", "priority:high"
                ])
                
                # Extract the issue number from the result
                match = re.search(r"#(\d+)", result)
                if match:
                    issue_number = match.group(1)
                    self.epics[category] = issue_number
                    self.epic_issues[category] = issue_number
                    
                    # Add to report
                    self.report["epic_issues"][category] = {
                        "number": issue_number,
                        "url": f"https://github.com/owner/repo/issues/{issue_number}"  # Placeholder URL
                    }
                    
                    print(f"  Created epic issue #{issue_number} for {category}")
                    
                    # Also add appropriate area label
                    if "Database" in category:
                        self.run_gh_command(["issue", "edit", issue_number, "--add-label", "area:database"], get_output=False)
                    elif "API" in category:
                        self.run_gh_command(["issue", "edit", issue_number, "--add-label", "area:api"], get_output=False)
                    elif "Frontend" in category:
                        self.run_gh_command(["issue", "edit", issue_number, "--add-label", "area:frontend"], get_output=False)
                    elif "Infrastructure" in category or "Deployment" in category:
                        self.run_gh_command(["issue", "edit", issue_number, "--add-label", "area:deployment"], get_output=False)
                    elif "Authentication" in category or "Security" in category:
                        self.run_gh_command(["issue", "edit", issue_number, "--add-label", "area:security"], get_output=False)
                    elif "Testing" in category:
                        self.run_gh_command(["issue", "edit", issue_number, "--add-label", "area:testing"], get_output=False)
                    elif "RDKit" in category:
                        self.run_gh_command(["issue", "edit", issue_number, "--add-label", "area:rdkit"], get_output=False)
                    
            except subprocess.CalledProcessError as e:
                print(f"  Error creating epic issue for {category}: {e}")
    
    def suggest_epic_for_issue(self, issue):
        """Suggest which epic category an issue belongs to"""
        title = issue["title"].lower()
        body = issue.get("body", "").lower() if issue.get("body") else ""
        combined_text = title + " " + body
        
        # Check if already has an epic label
        labels = [label["name"] for label in issue.get("labels", [])]
        for label in labels:
            if label.startswith("epic:"):
                epic_name = label[5:]  # Remove "epic:" prefix
                # Convert from slug to title case
                epic_name = " ".join(word.capitalize() for word in epic_name.split("-"))
                for category in EPIC_CATEGORIES:
                    if epic_name in category or category in epic_name:
                        return category
        
        # Score each category based on keyword matches
        scores = {category: 0 for category in EPIC_CATEGORIES}
        
        for category, keywords in CATEGORY_KEYWORDS.items():
            for keyword in keywords:
                if keyword in title:
                    scores[category] += 3  # Title matches are weighted higher
                if keyword in body:
                    scores[category] += 1
        
        # Get the category with the highest score
        if max(scores.values()) > 0:
            return max(scores.items(), key=lambda x: x[1])[0]
        
        return None
    
    def analyze_issues(self):
        """Analyze issues and suggest epic categories"""
        print("\nAnalyzing issues and suggesting epic categories...")
        
        for issue in self.issues:
            # Skip if it's an epic issue itself
            labels = [label["name"] for label in issue.get("labels", [])]
            if "type:epic" in labels:
                continue
            
            issue_num = issue["number"]
            suggested_epic = self.suggest_epic_for_issue(issue)
            
            if suggested_epic:
                self.issue_epic_map[issue_num] = suggested_epic
                
                # Initialize epic category in report if not exists
                if suggested_epic not in self.report["issues_by_epic"]:
                    self.report["issues_by_epic"][suggested_epic] = []
                
                # Add issue to report
                self.report["issues_by_epic"][suggested_epic].append({
                    "number": issue_num,
                    "title": issue["title"],
                    "url": issue["url"]
                })
                
                print(f"  Issue #{issue_num}: {issue['title']} -> {suggested_epic}")
            else:
                print(f"  Issue #{issue_num}: {issue['title']} -> [No suggestion]")
                
                # Add to unclassified issues
                self.report["unclassified_issues"].append({
                    "number": issue_num,
                    "title": issue["title"],
                    "url": issue["url"]
                })
        
        # Count issues by epic
        total_classified = sum(len(issues) for issues in self.report["issues_by_epic"].values())
        total_unclassified = len(self.report["unclassified_issues"])
        
        print(f"\nClassification summary:")
        print(f"  Total issues analyzed: {len(self.issues) - len(self.epics)}")
        print(f"  Issues classified: {total_classified}")
        print(f"  Issues unclassified: {total_unclassified}")
        
        for epic, issues in self.report["issues_by_epic"].items():
            print(f"  {epic}: {len(issues)} issues")
    
    def organize_issues(self):
        """Organize issues under their epic categories"""
        if self.dry_run:
            print("\n[DRY RUN] Would organize issues under epics:")
            for issue_num, epic in self.issue_epic_map.items():
                issue = next((i for i in self.issues if i["number"] == issue_num), None)
                if not issue:
                    continue
                
                print(f"  Would link issue #{issue_num} ({issue['title']}) to epic {epic}")
            return
        
        print("\nOrganizing issues under epics...")
        
        for issue_num, epic in self.issue_epic_map.items():
            issue = next((i for i in self.issues if i["number"] == issue_num), None)
            if not issue:
                continue
            
            if self.interactive:
                print(f"\nIssue #{issue_num}: {issue['title']}")
                print(f"Suggested epic: {epic}")
                choice = input("Link this issue to the suggested epic? [y/N]: ").lower()
                if choice != 'y':
                    print("  Skipping this issue")
                    continue
            
            # Check if we have an epic issue for this category
            if epic in self.epic_issues:
                epic_issue_num = self.epic_issues[epic]
                
                # Add comment linking to epic
                try:
                    comment_body = f"This issue is part of the [{epic}](https://github.com/owner/repo/issues/{epic_issue_num}) epic."
                    self.run_gh_command([
                        "issue", "comment", str(issue_num),
                        "--body", comment_body
                    ], get_output=False)
                    
                    # Add epic label
                    epic_label = f"epic:{epic.lower().replace(' ', '-')}"
                    self.run_gh_command([
                        "issue", "edit", str(issue_num),
                        "--add-label", epic_label
                    ], get_output=False)
                    
                    print(f"  Linked issue #{issue_num} to epic #{epic_issue_num}")
                    
                except subprocess.CalledProcessError as e:
                    print(f"  Error linking issue #{issue_num} to epic: {e}")
            else:
                print(f"  No epic issue found for {epic}, skipping issue #{issue_num}")
    
    def save_report(self):
        """Save the analysis report"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"output/epic_organization_report_{timestamp}.json"
        
        with open(filename, 'w') as f:
            json.dump(self.report, f, indent=2)
        
        print(f"\nReport saved to {filename}")
        
        # Also create a markdown version
        md_filename = f"output/epic_organization_report_{timestamp}.md"
        with open(md_filename, 'w') as f:
            f.write("# Epic Organization Report\n\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write("## Summary\n\n")
            f.write(f"- Total issues: {self.report['total_issues']}\n")
            f.write(f"- Epic issues: {len(self.report['epic_issues'])}\n")
            f.write(f"- Classified issues: {sum(len(issues) for issues in self.report['issues_by_epic'].values())}\n")
            f.write(f"- Unclassified issues: {len(self.report['unclassified_issues'])}\n\n")
            
            f.write("## Epic Issues\n\n")
            for epic, data in self.report["epic_issues"].items():
                f.write(f"- #{data['number']}: [{epic}]({data['url']})\n")
            
            f.write("\n## Issues by Epic\n\n")
            for epic, issues in self.report["issues_by_epic"].items():
                f.write(f"### {epic}\n\n")
                for issue in issues:
                    f.write(f"- #{issue['number']}: [{issue['title']}]({issue['url']})\n")
                f.write("\n")
            
            if self.report["unclassified_issues"]:
                f.write("## Unclassified Issues\n\n")
                for issue in self.report["unclassified_issues"]:
                    f.write(f"- #{issue['number']}: [{issue['title']}]({issue['url']})\n")
        
        print(f"Markdown report saved to {md_filename}")


def main():
    parser = argparse.ArgumentParser(description="GitHub Issue Epic Organizer")
    parser.add_argument("--dry-run", action="store_true", help="Run in dry run mode (no changes applied)")
    parser.add_argument("--create-epics", action="store_true", help="Create any missing epic issues")
    parser.add_argument("--interactive", action="store_true", help="Prompt for confirmation before making changes")
    
    args = parser.parse_args()
    
    # Check if GitHub CLI is installed and authenticated
    try:
        subprocess.run(["gh", "--version"], check=True, capture_output=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("Error: GitHub CLI (gh) is not installed or not in PATH.")
        print("Please install it from https://cli.github.com/")
        sys.exit(1)
    
    try:
        subprocess.run(["gh", "auth", "status"], check=True, capture_output=True)
    except subprocess.CalledProcessError:
        print("Error: GitHub CLI is not authenticated.")
        print("Please run 'gh auth login' to authenticate.")
        sys.exit(1)
    
    # Initialize organizer
    organizer = EpicOrganizer(
        dry_run=args.dry_run,
        create_epics=args.create_epics,
        interactive=args.interactive
    )
    
    try:
        # Fetch issues
        organizer.fetch_issues()
        
        # Analyze issues
        organizer.analyze_issues()
        
        # Organize issues
        organizer.organize_issues()
        
        # Save report
        organizer.save_report()
        
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()