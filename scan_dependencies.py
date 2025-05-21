import subprocess
import json
import os
from datetime import datetime

def run_command(cmd, output_file=None):
    """Run a shell command and optionally write output to a file."""
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=False)
        if output_file:
            with open(output_file, "w", encoding="utf-8") as f:
                f.write(result.stdout)
        return result.stdout
    except Exception as e:
        return str(e)

def run_bandit():
    output_file = "bandit-scan.json"
    cmd = "bandit -r ./api ./scripts ./ --exclude tests -f json"
    print(f"Running Bandit: {cmd}")
    run_command(cmd, output_file)
    return output_file

def run_safety():
    output_file = "safety-scan.json"
    req_file = "requirements.txt"
    if not os.path.exists(req_file):
        print("requirements.txt not found, skipping Safety scan.")
        return None
    cmd = f"safety check -r {req_file} --json"
    print(f"Running Safety: {cmd}")
    run_command(cmd, output_file)
    return output_file

def run_dependency_check():
    output_file = "dependency-check-scan.json"
    # Assumes dependency-check is installed and available in PATH
    cmd = (
        "dependency-check --project 'CryoProtect v2' --scan . --format JSON "
        f"--out {output_file}"
    )
    print(f"Running OWASP Dependency-Check: {cmd}")
    run_command(cmd)
    return os.path.join(output_file, "dependency-check-report.json") if os.path.isdir(output_file) else output_file

def parse_critical_vulns_bandit(bandit_file):
    try:
        with open(bandit_file, "r", encoding="utf-8") as f:
            data = json.load(f)
        critical = [issue for issue in data.get("results", []) if issue.get("issue_severity") == "HIGH"]
        return critical
    except Exception:
        return []

def parse_critical_vulns_safety(safety_file):
    try:
        with open(safety_file, "r", encoding="utf-8") as f:
            data = json.load(f)
        critical = [v for v in data if v.get("severity", "").lower() == "high" or v.get("vuln_type", "").lower() == "vulnerability"]
        return critical
    except Exception:
        return []

def parse_critical_vulns_dependency_check(dep_file):
    try:
        with open(dep_file, "r", encoding="utf-8") as f:
            data = json.load(f)
        critical = []
        for dep in data.get("dependencies", []):
            for vuln in dep.get("vulnerabilities", []):
                if vuln.get("severity", "").upper() in ("CRITICAL", "HIGH"):
                    critical.append({
                        "dependency": dep.get("fileName"),
                        "vulnerability": vuln
                    })
        return critical
    except Exception:
        return []

def main():
    print("=== CryoProtect v2 Runtime Dependency Scan ===")
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_file = f"scan_report_{timestamp}.log"
    alert_file = f"scan_alert_{timestamp}.json"

    bandit_file = run_bandit()
    safety_file = run_safety()
    dep_file = run_dependency_check()

    critical_bandit = parse_critical_vulns_bandit(bandit_file)
    critical_safety = parse_critical_vulns_safety(safety_file)
    critical_dep = parse_critical_vulns_dependency_check(dep_file)

    summary = {
        "timestamp": timestamp,
        "bandit_critical": critical_bandit,
        "safety_critical": critical_safety,
        "dependency_check_critical": critical_dep,
    }

    # Write summary log
    with open(log_file, "w", encoding="utf-8") as f:
        f.write("=== CryoProtect v2 Vulnerability Scan Report ===\n")
        f.write(f"Timestamp: {timestamp}\n\n")
        f.write(f"Bandit critical issues: {len(critical_bandit)}\n")
        f.write(f"Safety critical issues: {len(critical_safety)}\n")
        f.write(f"Dependency-Check critical issues: {len(critical_dep)}\n\n")
        if critical_bandit:
            f.write("Bandit critical findings:\n")
            for issue in critical_bandit:
                f.write(json.dumps(issue, indent=2) + "\n")
        if critical_safety:
            f.write("Safety critical findings:\n")
            for issue in critical_safety:
                f.write(json.dumps(issue, indent=2) + "\n")
        if critical_dep:
            f.write("Dependency-Check critical findings:\n")
            for issue in critical_dep:
                f.write(json.dumps(issue, indent=2) + "\n")

    # Write alert file if any critical issues found
    if critical_bandit or critical_safety or critical_dep:
        with open(alert_file, "w", encoding="utf-8") as f:
            json.dump(summary, f, indent=2)
        print(f"ALERT: Critical vulnerabilities found! See {alert_file}")
    else:
        print("No critical vulnerabilities found.")

    print(f"Scan complete. See {log_file} for details.")

if __name__ == "__main__":
    main()