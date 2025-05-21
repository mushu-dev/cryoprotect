#!/usr/bin/env python3
"""
AI-Generated Issue Detection Script

This script helps identify potentially AI-generated issues by analyzing
patterns, repetition, and other indicators of mass-generated content.

Usage:
    python detect-ai-spam.py [--threshold=0.7] [--close-confirmed]

Options:
    --threshold     Confidence threshold for flagging issues (default: 0.7)
    --close-confirmed   Automatically close issues confirmed as spam (default: False)
"""

import argparse
import json
import re
import subprocess
import sys
from datetime import datetime, timedelta
from typing import Dict, List, Tuple

# Define patterns that indicate AI-generated content
AI_PATTERNS = [
    r"test\s+issue",
    r"lorem\s+ipsum",
    r"this\s+is\s+a\s+test",
    r"sample\s+issue",
    r"demo\s+issue",
    r"placeholder",
    r"verifying\s+functionality"
]

# Define heuristics for detecting AI-generated issues
HEURISTICS = {
    "generic_title": 0.3,  # Generic one-word or two-word titles
    "pattern_match": 0.4,  # Matches known AI patterns
    "repetitive_content": 0.3,  # Issue description repeats title or has repetitive phrases
    "creation_time_cluster": 0.2,  # Created in a cluster with many other issues
    "lacks_context": 0.3,  # No references to actual code or specific details
    "template_unmodified": 0.5,  # Issue template appears unmodified
}

def run_gh_command(cmd: List[str]) -> Dict:
    """Run a GitHub CLI command and return the JSON result."""
    try:
        result = subprocess.run(
            ["gh"] + cmd,
            capture_output=True,
            text=True,
            check=True
        )
        return json.loads(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"Error running GitHub CLI command: {e}")
        print(f"Error output: {e.stderr}")
        sys.exit(1)
    except json.JSONDecodeError:
        print(f"Error parsing JSON from GitHub CLI output.")
        sys.exit(1)

def get_all_issues(days_back: int = 30) -> List[Dict]:
    """Get all issues from the last X days."""
    date_threshold = (datetime.now() - timedelta(days=days_back)).strftime("%Y-%m-%d")
    return run_gh_command([
        "issue", 
        "list", 
        "--json", "number,title,body,createdAt,author,labels", 
        "--search", f"created:>={date_threshold}"
    ])

def is_generic_title(title: str) -> bool:
    """Check if the title is too generic."""
    word_count = len(title.split())
    return word_count <= 2 and all(len(word) < 10 for word in title.split())

def matches_ai_patterns(text: str) -> bool:
    """Check if the text matches known AI patterns."""
    text_lower = text.lower()
    for pattern in AI_PATTERNS:
        if re.search(pattern, text_lower):
            return True
    return False

def has_repetitive_content(title: str, body: str) -> bool:
    """Check if the issue has repetitive content."""
    # Check if title is repeated in body
    if title.lower() in body.lower():
        return True
    
    # Check for repeated phrases
    words = body.lower().split()
    if len(words) < 10:
        return False
    
    # Check for repeated sentences
    sentences = re.split(r'[.!?]', body)
    sentences = [s.strip() for s in sentences if s.strip()]
    
    for i in range(len(sentences)):
        for j in range(i+1, len(sentences)):
            if sentences[i] and sentences[j] and sentences[i] == sentences[j]:
                return True
    
    return False

def part_of_creation_cluster(created_at: str, all_issues: List[Dict]) -> bool:
    """Check if the issue was created in a cluster with many other issues."""
    created_time = datetime.fromisoformat(created_at.replace("Z", "+00:00"))
    
    # Count issues created within 10 minutes
    cluster_count = 0
    for issue in all_issues:
        other_time = datetime.fromisoformat(issue["createdAt"].replace("Z", "+00:00"))
        if abs((created_time - other_time).total_seconds()) < 600:  # 10 minutes
            cluster_count += 1
    
    return cluster_count >= 5  # If 5 or more issues were created within 10 minutes

def lacks_context(body: str) -> bool:
    """Check if the issue lacks specific context about the project."""
    # Look for file paths, code snippets, etc.
    has_code_snippet = "```" in body
    has_file_reference = re.search(r'[a-zA-Z0-9_-]+\.[a-zA-Z]{1,5}', body) is not None
    has_specific_errors = re.search(r'error|exception|bug|crash|fail', body.lower()) is not None
    
    return not (has_code_snippet or has_file_reference or has_specific_errors)

def has_unmodified_template(body: str) -> bool:
    """Check if the issue appears to have an unmodified template."""
    template_indicators = [
        "## Description",
        "## Steps to Reproduce",
        "## Expected Behavior",
        "1. Go to '...'",
        "Add any other context",
        "A clear and concise description of what you expected to happen",
        "## Actual Behavior"
    ]
    
    matches = 0
    for indicator in template_indicators:
        if indicator in body:
            matches += 1
    
    return matches >= 3  # If 3 or more template phrases are present

def analyze_issue(issue: Dict, all_issues: List[Dict]) -> Tuple[float, Dict[str, bool]]:
    """Analyze an issue and return the AI spam confidence score and reasons."""
    title = issue.get("title", "")
    body = issue.get("body", "")
    created_at = issue.get("createdAt", "")
    
    # Apply heuristics
    score = 0.0
    reasons = {}
    
    if is_generic_title(title):
        score += HEURISTICS["generic_title"]
        reasons["generic_title"] = True
    
    if matches_ai_patterns(title) or matches_ai_patterns(body):
        score += HEURISTICS["pattern_match"]
        reasons["pattern_match"] = True
    
    if has_repetitive_content(title, body):
        score += HEURISTICS["repetitive_content"]
        reasons["repetitive_content"] = True
    
    if part_of_creation_cluster(created_at, all_issues):
        score += HEURISTICS["creation_time_cluster"]
        reasons["creation_time_cluster"] = True
    
    if lacks_context(body):
        score += HEURISTICS["lacks_context"]
        reasons["lacks_context"] = True
    
    if has_unmodified_template(body):
        score += HEURISTICS["template_unmodified"]
        reasons["template_unmodified"] = True
    
    return score, reasons

def main():
    parser = argparse.ArgumentParser(description="Detect AI-generated GitHub issues")
    parser.add_argument("--threshold", type=float, default=0.7, 
                        help="Confidence threshold for flagging issues (default: 0.7)")
    parser.add_argument("--close-confirmed", action="store_true",
                        help="Automatically close issues confirmed as spam")
    args = parser.parse_args()
    
    print("Fetching recent issues...")
    all_issues = get_all_issues(30)  # Get issues from the last 30 days
    
    print(f"Analyzing {len(all_issues)} issues for AI-generated content...")
    print("")
    
    # Analyze each issue
    potential_spam = []
    
    for issue in all_issues:
        score, reasons = analyze_issue(issue, all_issues)
        if score >= args.threshold:
            issue["spam_score"] = score
            issue["spam_reasons"] = reasons
            potential_spam.append(issue)
    
    # Print results
    if not potential_spam:
        print("No potential AI-generated issues found.")
        return
    
    print(f"Found {len(potential_spam)} potential AI-generated issues:")
    print("")
    
    for issue in sorted(potential_spam, key=lambda x: x["spam_score"], reverse=True):
        print(f"#{issue['number']} - {issue['title']}")
        print(f"  Created by: {issue['author']['login']}")
        print(f"  Created at: {issue['createdAt']}")
        print(f"  Spam confidence: {issue['spam_score']:.2f}")
        print(f"  Reasons:")
        for reason, value in issue['spam_reasons'].items():
            if value:
                print(f"    - {reason.replace('_', ' ').title()}")
        print("")
    
    # Handle closing confirmed spam issues
    if args.close_confirmed and potential_spam:
        high_confidence_spam = [i for i in potential_spam if i["spam_score"] > 0.85]
        
        if not high_confidence_spam:
            print("No high-confidence spam issues to close.")
            return
        
        print(f"Found {len(high_confidence_spam)} high-confidence spam issues.")
        close_input = input("Do you want to close these issues? (y/n): ")
        
        if close_input.lower() == 'y':
            for issue in high_confidence_spam:
                issue_number = issue["number"]
                print(f"Closing issue #{issue_number}...")
                
                try:
                    subprocess.run([
                        "gh", "issue", "close", str(issue_number),
                        "--comment", "Closing as this appears to be an AI-generated test issue."
                    ], check=True)
                    
                    # Add a label
                    subprocess.run([
                        "gh", "issue", "edit", str(issue_number),
                        "--add-label", "ai-generated,spam"
                    ], check=True)
                    
                    print(f"Successfully closed issue #{issue_number}")
                except subprocess.CalledProcessError as e:
                    print(f"Error closing issue #{issue_number}: {e}")
            
            print("Finished closing spam issues.")

if __name__ == "__main__":
    main()