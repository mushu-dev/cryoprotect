#!/usr/bin/env python
"""
GitHub Issue Audit Script

This script audits GitHub issues in the repository to ensure they are properly tagged,
organized, and consistent with the current project structure. It helps maintain a clean
project board and avoid confusion with older issues.

Usage:
    python audit_issues.py [--dry-run] [--fix]

Options:
    --dry-run    Run in dry run mode (no changes applied)
    --fix        Automatically fix issues where possible
"""

import argparse
import os
import sys
import json
import subprocess
from datetime import datetime

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

# Area labels
AREA_LABELS = [
    "area:database",
    "area:api",
    "area:frontend",
    "area:deployment",
    "area:security",
    "area:testing",
    "area:documentation",
    "area:rdkit",
    "area:infrastructure"
]

# Status labels
STATUS_LABELS = [
    "status:planning",
    "status:ready",
    "status:in-progress",
    "status:review",
    "status:blocked",
    "status:completed"
]

# Priority labels
PRIORITY_LABELS = [
    "priority:critical",
    "priority:high",
    "priority:medium",
    "priority:low"
]

def run_gh_command(cmd):
    """Run a GitHub CLI command and return the output as JSON"""
    try:
        result = subprocess.run(
            ["gh"] + cmd, 
            check=True, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE,
            text=True
        )
        return json.loads(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {e}")
        print(f"stderr: {e.stderr}")
        return None
    except json.JSONDecodeError:
        print(f"Error parsing JSON output from command: {' '.join(cmd)}")
        return None

def get_all_issues():
    """Get all issues from the repository"""
    return run_gh_command(["issue", "list", "--state", "all", "--json", 
                          "number,title,state,labels,milestone,assignees,author,createdAt,updatedAt,closedAt,url"])

def get_project_items():
    """Get all items in the project board"""
    # Note: This requires the project ID or URL to be configured
    # Adjust this command based on your project structure
    return run_gh_command(["project", "item-list", "--format", "json"])

def analyze_issue(issue):
    """Analyze an issue and identify any problems"""
    problems = []
    
    # Extract labels
    labels = [label["name"] for label in issue.get("labels", [])]
    
    # Check if issue has at least one area label
    has_area_label = any(label in AREA_LABELS for label in labels)
    if not has_area_label:
        problems.append("Missing area label")
    
    # Check if issue has a status label
    has_status_label = any(label in STATUS_LABELS for label in labels)
    if not has_status_label:
        problems.append("Missing status label")
    
    # Check if issue has a priority label
    has_priority_label = any(label in PRIORITY_LABELS for label in labels)
    if not has_priority_label:
        problems.append("Missing priority label")
    
    # Check for stale issues (open but not updated in 30 days)
    if issue["state"] == "OPEN":
        updated_at = datetime.fromisoformat(issue["updatedAt"].replace("Z", "+00:00"))
        now = datetime.now().astimezone()
        days_since_update = (now - updated_at).days
        if days_since_update > 30:
            problems.append(f"Stale issue (not updated in {days_since_update} days)")
    
    # Check if closed issue has a resolution label or comment
    if issue["state"] == "CLOSED" and issue["closedAt"]:
        has_resolution = "resolution:completed" in labels or "resolution:wontfix" in labels or "resolution:duplicate" in labels
        if not has_resolution:
            problems.append("Closed issue missing resolution label")
    
    return problems

def suggest_fixes(issue, problems):
    """Suggest fixes for issue problems"""
    fixes = []
    issue_number = issue["number"]
    
    if "Missing area label" in problems:
        # Try to determine area from title or content
        title = issue["title"].lower()
        if any(kw in title for kw in ["database", "db", "supabase", "sql", "query"]):
            fixes.append(f"gh issue edit {issue_number} --add-label 'area:database'")
        elif any(kw in title for kw in ["api", "endpoint", "route", "flask"]):
            fixes.append(f"gh issue edit {issue_number} --add-label 'area:api'")
        elif any(kw in title for kw in ["frontend", "ui", "ux", "react", "next", "vercel"]):
            fixes.append(f"gh issue edit {issue_number} --add-label 'area:frontend'")
        elif any(kw in title for kw in ["deploy", "production", "heroku", "cloud"]):
            fixes.append(f"gh issue edit {issue_number} --add-label 'area:deployment'")
        elif any(kw in title for kw in ["security", "auth", "jwt", "token", "permission"]):
            fixes.append(f"gh issue edit {issue_number} --add-label 'area:security'")
        elif any(kw in title for kw in ["test", "ci", "quality", "coverage", "spec"]):
            fixes.append(f"gh issue edit {issue_number} --add-label 'area:testing'")
        elif any(kw in title for kw in ["rdkit", "molecule", "chemical", "property"]):
            fixes.append(f"gh issue edit {issue_number} --add-label 'area:rdkit'")
        else:
            fixes.append(f"gh issue edit {issue_number} --add-label 'needs:triage'")
    
    if "Missing status label" in problems:
        if issue["state"] == "OPEN":
            fixes.append(f"gh issue edit {issue_number} --add-label 'status:planning'")
    
    if "Missing priority label" in problems:
        fixes.append(f"gh issue edit {issue_number} --add-label 'priority:medium'")
    
    if "Stale issue" in problems:
        fixes.append(f"gh issue comment {issue_number} --body 'This issue appears to be stale. Please update or consider closing it.'")
        fixes.append(f"gh issue edit {issue_number} --add-label 'status:stale'")
    
    if "Closed issue missing resolution label" in problems:
        fixes.append(f"gh issue edit {issue_number} --add-label 'resolution:completed'")
    
    return fixes

def apply_fixes(fixes, dry_run=True):
    """Apply fixes to issues"""
    if dry_run:
        print("\nSuggested fixes (dry run):")
        for fix in fixes:
            print(f"  {fix}")
    else:
        print("\nApplying fixes:")
        for fix in fixes:
            print(f"  {fix}")
            try:
                subprocess.run(fix, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(f"    Error applying fix: {e}")

def main():
    parser = argparse.ArgumentParser(description="Audit GitHub issues")
    parser.add_argument("--dry-run", action="store_true", help="Run in dry run mode (no changes applied)")
    parser.add_argument("--fix", action="store_true", help="Automatically fix issues where possible")
    
    args = parser.parse_args()
    
    print("Fetching issues...")
    issues = get_all_issues()
    
    if not issues:
        print("No issues found or error fetching issues")
        return
    
    print(f"Found {len(issues)} issues")
    
    issues_with_problems = []
    all_fixes = []
    
    for issue in issues:
        problems = analyze_issue(issue)
        if problems:
            issues_with_problems.append((issue, problems))
            fixes = suggest_fixes(issue, problems)
            all_fixes.extend(fixes)
    
    # Generate audit report
    report = {
        "total_issues": len(issues),
        "issues_with_problems": len(issues_with_problems),
        "timestamp": datetime.now().isoformat(),
        "problems_by_type": {},
        "issue_details": []
    }
    
    # Count problem types
    for issue, problems in issues_with_problems:
        for problem in problems:
            if problem not in report["problems_by_type"]:
                report["problems_by_type"][problem] = 0
            report["problems_by_type"][problem] += 1
    
    # Add issue details
    for issue, problems in issues_with_problems:
        report["issue_details"].append({
            "number": issue["number"],
            "title": issue["title"],
            "url": issue["url"],
            "state": issue["state"],
            "problems": problems
        })
    
    # Print summary
    print(f"\nIssues with problems: {len(issues_with_problems)} / {len(issues)}")
    print("\nProblem types:")
    for problem, count in report["problems_by_type"].items():
        print(f"  {problem}: {count}")
    
    # Print detailed report
    print("\nDetailed issue problems:")
    for issue, problems in issues_with_problems:
        print(f"  #{issue['number']} - {issue['title']}")
        for problem in problems:
            print(f"    - {problem}")
    
    # Apply fixes if requested
    if args.fix:
        apply_fixes(all_fixes, dry_run=False)
    else:
        apply_fixes(all_fixes, dry_run=True)
    
    # Save report to file
    report_file = f"issue_audit_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    with open(report_file, "w") as f:
        json.dump(report, f, indent=2)
    
    print(f"\nAudit report saved to {report_file}")
    
    # Recommend cleanup actions
    print("\nRecommended actions:")
    print("1. Review and verify all applied label changes")
    print("2. Close or update stale issues")
    print("3. Add issues to the project board if missing")
    print("4. Group related issues under epics")
    print("5. Update milestone assignments")

if __name__ == "__main__":
    main()