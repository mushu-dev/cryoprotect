#!/usr/bin/env python
"""
GitHub Issue Update Audit Script

This script identifies issues that need updating to conform to the current
project structure, organization standards, and epic categories. It complements
the existing audit tools by focusing specifically on:

1. Ensuring issues are properly categorized into the 9 epic categories
2. Updating stale issues and closing those that are no longer relevant
3. Ensuring all issues have consistent labels according to current standards
4. Cleaning up old issues from before the project reorganization

Usage:
    python update_issue_audit.py [--dry-run] [--apply]

Options:
    --dry-run    Run in dry run mode (shows changes but doesn't apply them)
    --apply      Apply recommended changes after confirmation
"""

import argparse
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from typing import Dict, List, Tuple, Optional, Any

# Current epic categories
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

# Current label schema
AREA_LABELS = [
    "area:database", "area:api", "area:frontend", "area:deployment", 
    "area:security", "area:testing", "area:documentation", 
    "area:rdkit", "area:infrastructure"
]

STATUS_LABELS = [
    "status:planning", "status:ready", "status:in-progress", 
    "status:review", "status:blocked", "status:completed", "status:stale"
]

PRIORITY_LABELS = [
    "priority:critical", "priority:high", "priority:medium", "priority:low"
]

RESOLUTION_LABELS = [
    "resolution:completed", "resolution:wontfix", "resolution:duplicate"
]

# Phase labels (old system)
OLD_PHASE_LABELS = [
    "phase:1", "phase:2", "phase:2.1", "phase:2.2", "phase:2.3",
    "phase:3", "phase:3.1", "phase:3.2", "phase:3.3", "phase:4", "phase:5"
]

class IssueAuditor:
    def __init__(self, dry_run: bool = True):
        self.dry_run = dry_run
        self.issues = []
        self.epics = []
        self.updated_count = 0
        self.closed_count = 0
        self.modifications = []
        self.report = {
            "timestamp": datetime.now().isoformat(),
            "total_issues": 0,
            "issues_requiring_updates": 0,
            "issues_closed": 0,
            "detailed_changes": []
        }
    
    def fetch_issues(self):
        """Fetch all issues from GitHub"""
        print("Fetching issues from GitHub...")
        result = subprocess.run(
            ["gh", "issue", "list", "--state", "all", "--json", 
             "number,title,state,labels,milestone,assignees,author,createdAt,updatedAt,closedAt,body,url"],
            capture_output=True, text=True, check=True
        )
        self.issues = json.loads(result.stdout)
        self.report["total_issues"] = len(self.issues)
        print(f"Found {len(self.issues)} issues")
        
        # Identify epic issues
        self.epics = [issue for issue in self.issues 
                       if any(label["name"] == "type:epic" for label in issue.get("labels", []))]
        print(f"Found {len(self.epics)} epic issues")
        
        return self.issues
    
    def identify_epic_category(self, issue: Dict) -> Optional[str]:
        """Identify which epic category an issue belongs to based on title and content"""
        title = issue["title"].lower()
        body = issue.get("body", "").lower()
        
        # Map of keywords to epic categories
        category_keywords = {
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
        
        # Score each category based on keyword matches
        scores = {category: 0 for category in EPIC_CATEGORIES}
        
        for category, keywords in category_keywords.items():
            for keyword in keywords:
                if keyword in title:
                    scores[category] += 3  # Title matches are weighted higher
                if keyword in body:
                    scores[category] += 1
        
        # Get the category with the highest score
        if max(scores.values()) > 0:
            return max(scores.items(), key=lambda x: x[1])[0]
        
        return None
    
    def analyze_issue(self, issue: Dict) -> Dict:
        """Analyze an issue and suggest updates"""
        issue_number = issue["number"]
        title = issue["title"]
        state = issue["state"]
        labels = [label["name"] for label in issue.get("labels", [])]
        
        updates = {
            "issue_number": issue_number,
            "title": title,
            "url": issue["url"],
            "current_state": state,
            "current_labels": labels,
            "needed_updates": [],
            "recommended_labels": [],
            "recommended_actions": [],
            "suggested_epic": None,
            "is_stale": False
        }
        
        # Check if issue belongs to old phase labeling system
        has_old_phase_label = any(label in OLD_PHASE_LABELS for label in labels)
        
        # Check if issue has required labels
        has_area_label = any(label in AREA_LABELS for label in labels)
        has_status_label = any(label in STATUS_LABELS for label in labels)
        has_priority_label = any(label in PRIORITY_LABELS for label in labels)
        
        # Check if issue is stale
        if state == "OPEN":
            updated_at = datetime.fromisoformat(issue["updatedAt"].replace("Z", "+00:00"))
            now = datetime.now(timezone.utc)
            days_since_update = (now - updated_at).days
            
            if days_since_update > 30:
                updates["is_stale"] = True
                updates["needed_updates"].append(f"Stale issue (not updated in {days_since_update} days)")
                if "status:stale" not in labels:
                    updates["recommended_labels"].append("status:stale")
        
        # Identify which epic category this issue belongs to
        suggested_epic = self.identify_epic_category(issue)
        if suggested_epic:
            updates["suggested_epic"] = suggested_epic
        
        # Check for label updates needed
        if not has_area_label:
            updates["needed_updates"].append("Missing area label")
            # Suggest area label based on content analysis
            if suggested_epic:
                if "Database" in suggested_epic:
                    updates["recommended_labels"].append("area:database")
                elif "API" in suggested_epic:
                    updates["recommended_labels"].append("area:api")
                elif "Frontend" in suggested_epic:
                    updates["recommended_labels"].append("area:frontend")
                elif "Infrastructure" in suggested_epic or "Deployment" in suggested_epic:
                    updates["recommended_labels"].append("area:deployment")
                elif "Authentication" in suggested_epic or "Security" in suggested_epic:
                    updates["recommended_labels"].append("area:security")
                elif "Testing" in suggested_epic:
                    updates["recommended_labels"].append("area:testing")
                elif "RDKit" in suggested_epic:
                    updates["recommended_labels"].append("area:rdkit")
        
        if not has_status_label:
            updates["needed_updates"].append("Missing status label")
            if state == "OPEN":
                updates["recommended_labels"].append("status:planning")
            elif state == "CLOSED":
                updates["recommended_labels"].append("status:completed")
        
        if not has_priority_label:
            updates["needed_updates"].append("Missing priority label")
            updates["recommended_labels"].append("priority:medium")
        
        # Check for closed issues missing resolution labels
        if state == "CLOSED" and not any(label in RESOLUTION_LABELS for label in labels):
            updates["needed_updates"].append("Closed issue missing resolution label")
            updates["recommended_labels"].append("resolution:completed")
        
        # Check for old phase labels that need to be updated
        if has_old_phase_label:
            updates["needed_updates"].append("Uses old phase labeling system")
            updates["recommended_actions"].append("Remove old phase labels and categorize into epic")
        
        # If issue is very old and has no activity, suggest closing
        created_at = datetime.fromisoformat(issue["createdAt"].replace("Z", "+00:00"))
        now = datetime.now(timezone.utc)
        days_since_creation = (now - created_at).days
        
        if state == "OPEN" and days_since_creation > 90 and updates["is_stale"]:
            updates["recommended_actions"].append("Close issue due to inactivity (created over 90 days ago and stale)")
        
        return updates
    
    def audit_issues(self):
        """Audit all issues and generate update recommendations"""
        print("Auditing issues...")
        updates_needed = []
        
        for issue in self.issues:
            analysis = self.analyze_issue(issue)
            if analysis["needed_updates"] or analysis["recommended_actions"]:
                updates_needed.append(analysis)
        
        self.report["issues_requiring_updates"] = len(updates_needed)
        print(f"Found {len(updates_needed)} issues requiring updates")
        
        return updates_needed
    
    def apply_updates(self, updates_needed: List[Dict], confirm: bool = True):
        """Apply the recommended updates to issues"""
        if self.dry_run:
            print("\nDRY RUN MODE - No changes will be applied")
            for update in updates_needed:
                print(f"\nIssue #{update['issue_number']}: {update['title']}")
                print(f"  URL: {update['url']}")
                print(f"  Current state: {update['current_state']}")
                print(f"  Current labels: {', '.join(update['current_labels'])}")
                print(f"  Needed updates: {', '.join(update['needed_updates'])}")
                if update['recommended_labels']:
                    print(f"  Recommended labels: {', '.join(update['recommended_labels'])}")
                if update['recommended_actions']:
                    print(f"  Recommended actions: {', '.join(update['recommended_actions'])}")
                if update['suggested_epic']:
                    print(f"  Suggested epic category: {update['suggested_epic']}")
            return
        
        if confirm:
            print(f"\nReady to apply updates to {len(updates_needed)} issues.")
            choice = input("Proceed? [y/N]: ").lower()
            if choice != 'y':
                print("Update canceled.")
                return
        
        for update in updates_needed:
            issue_num = update['issue_number']
            print(f"\nUpdating issue #{issue_num}: {update['title']}")
            
            # Add recommended labels
            if update['recommended_labels']:
                labels_to_add = ','.join(update['recommended_labels'])
                print(f"  Adding labels: {labels_to_add}")
                try:
                    subprocess.run(
                        ["gh", "issue", "edit", str(issue_num), "--add-label", labels_to_add],
                        check=True, capture_output=True
                    )
                    self.updated_count += 1
                    self.modifications.append(f"Added labels to issue #{issue_num}: {labels_to_add}")
                except subprocess.CalledProcessError as e:
                    print(f"  Error adding labels: {e.stderr.decode()}")
            
            # Handle recommended actions
            for action in update['recommended_actions']:
                if "Close issue" in action and update['current_state'] == "OPEN":
                    print(f"  Closing issue as stale")
                    try:
                        comment = "Closing this issue due to inactivity. If this issue is still relevant, please reopen it and add updated information."
                        subprocess.run(
                            ["gh", "issue", "comment", str(issue_num), "--body", comment],
                            check=True, capture_output=True
                        )
                        subprocess.run(
                            ["gh", "issue", "close", str(issue_num), "--reason", "not_planned"],
                            check=True, capture_output=True
                        )
                        self.closed_count += 1
                        self.modifications.append(f"Closed issue #{issue_num} as stale")
                    except subprocess.CalledProcessError as e:
                        print(f"  Error closing issue: {e.stderr.decode()}")
                
                elif "Remove old phase labels" in action:
                    # Find old phase labels to remove
                    old_labels = [label for label in update['current_labels'] if label in OLD_PHASE_LABELS]
                    if old_labels:
                        labels_to_remove = ','.join(old_labels)
                        print(f"  Removing old phase labels: {labels_to_remove}")
                        try:
                            subprocess.run(
                                ["gh", "issue", "edit", str(issue_num), "--remove-label", labels_to_remove],
                                check=True, capture_output=True
                            )
                            self.modifications.append(f"Removed old labels from issue #{issue_num}: {labels_to_remove}")
                        except subprocess.CalledProcessError as e:
                            print(f"  Error removing labels: {e.stderr.decode()}")
                    
                    # Add epic category label if suggested
                    if update['suggested_epic']:
                        epic_label = f"epic:{update['suggested_epic'].lower().replace(' ', '-')}"
                        print(f"  Adding epic category: {epic_label}")
                        try:
                            subprocess.run(
                                ["gh", "issue", "edit", str(issue_num), "--add-label", epic_label],
                                check=True, capture_output=True
                            )
                            self.modifications.append(f"Added epic category to issue #{issue_num}: {epic_label}")
                        except subprocess.CalledProcessError as e:
                            print(f"  Error adding epic category: {e.stderr.decode()}")
            
            # Record the changes in the report
            self.report["detailed_changes"].append({
                "issue_number": issue_num,
                "title": update['title'],
                "url": update['url'],
                "changes_made": [
                    f"Added labels: {', '.join(update['recommended_labels'])}" if update['recommended_labels'] else None,
                    f"Categorized into epic: {update['suggested_epic']}" if update['suggested_epic'] else None,
                    "Closed as stale" if any("Close issue" in action for action in update['recommended_actions']) else None,
                    f"Removed old phase labels: {', '.join([label for label in update['current_labels'] if label in OLD_PHASE_LABELS])}" 
                    if any("Remove old phase labels" in action for action in update['recommended_actions']) else None
                ]
            })
        
        self.report["issues_closed"] = self.closed_count
        print(f"\nCompleted updates: {self.updated_count} issues updated, {self.closed_count} issues closed")
    
    def save_report(self):
        """Save the audit report to a file"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"issue_update_report_{timestamp}.json"
        
        # Add summary of modifications to report
        self.report["modifications"] = self.modifications
        
        with open(filename, 'w') as f:
            json.dump(self.report, f, indent=2)
        
        print(f"\nAudit report saved to {filename}")
        
        # Also save a markdown summary
        md_filename = f"issue_update_report_{timestamp}.md"
        with open(md_filename, 'w') as f:
            f.write(f"# GitHub Issue Update Report\n\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            f.write(f"## Summary\n\n")
            f.write(f"- Total issues: {self.report['total_issues']}\n")
            f.write(f"- Issues requiring updates: {self.report['issues_requiring_updates']}\n")
            f.write(f"- Issues closed: {self.report['issues_closed']}\n\n")
            
            f.write(f"## Changes Made\n\n")
            for mod in self.modifications:
                f.write(f"- {mod}\n")
            
            f.write(f"\n## Detailed Changes\n\n")
            for change in self.report["detailed_changes"]:
                f.write(f"### Issue #{change['issue_number']}: {change['title']}\n\n")
                f.write(f"- URL: {change['url']}\n")
                f.write("- Changes:\n")
                for item in change['changes_made']:
                    if item:
                        f.write(f"  - {item}\n")
                f.write("\n")
        
        print(f"Markdown summary saved to {md_filename}")


def main():
    parser = argparse.ArgumentParser(description="GitHub Issue Update Audit Script")
    parser.add_argument("--dry-run", action="store_true", help="Run in dry run mode (no changes applied)")
    parser.add_argument("--apply", action="store_true", help="Apply recommended changes after confirmation")
    
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
    
    auditor = IssueAuditor(dry_run=not args.apply)
    
    try:
        # Fetch and audit issues
        auditor.fetch_issues()
        updates_needed = auditor.audit_issues()
        
        # Apply updates if requested
        if updates_needed:
            auditor.apply_updates(updates_needed, confirm=True)
        else:
            print("No issues require updates.")
        
        # Save report
        auditor.save_report()
        
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()