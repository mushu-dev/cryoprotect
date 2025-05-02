import argparse
import json
import os
import sys
import logging
import datetime
from glob import glob

# Setup logging
LOG_FILE = "logs/validation_report.log"
os.makedirs(os.path.dirname(LOG_FILE), exist_ok=True)
logging.basicConfig(
    filename=LOG_FILE,
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
logger = logging.getLogger(__name__)

def find_latest_report(report_dir, prefix):
    """Find the latest report file in a directory with a given prefix."""
    pattern = os.path.join(report_dir, f"{prefix}_*.json")
    files = glob(pattern)
    if not files:
        return None
    latest_file = max(files, key=os.path.getmtime)
    return latest_file

def load_json_report(report_path):
    """Load a JSON report from file."""
    try:
        with open(report_path, "r") as f:
            return json.load(f)
    except Exception as e:
        logger.error(f"Failed to load report {report_path}: {e}")
        return None

def summarize_reports(verification, restore):
    """Create a summary dictionary from verification and restore reports."""
    summary = {
        "timestamp": datetime.datetime.now().isoformat(),
        "verification": verification,
        "restore": restore
    }
    return summary

def format_markdown(summary):
    """Format the summary as Markdown."""
    lines = [
        f"# Backup Validation Report",
        f"**Generated:** {summary['timestamp']}",
        "",
        "## Backup Verification",
    ]
    verification = summary.get("verification")
    if verification:
        lines += [
            f"- **Total Backups:** {verification.get('total_backups', 'N/A')}",
            f"- **Passed:** {verification.get('passed', 'N/A')}",
            f"- **Failed:** {verification.get('failed', 'N/A')}",
            "",
            "### Details:",
        ]
        for b in verification.get("verified_backups", []):
            lines.append(
                f"- `{b.get('backup_name', 'N/A')}`: **{b.get('status', 'N/A').upper()}**"
            )
    else:
        lines.append("_No verification report found._")

    lines += [
        "",
        "## Restore Test",
    ]
    restore = summary.get("restore")
    if restore:
        lines += [
            f"- **Backup File:** {restore.get('backup_file', 'N/A')}",
            f"- **Restore Success:** {restore.get('success', 'N/A')}",
            f"- **Error:** {restore.get('error', 'None')}",
        ]
    else:
        lines.append("_No restore report found._")
    return "\n".join(lines)

# Notification stubs
def notify_email(summary):
    logger.info("Stub: Sending email notification (not implemented).")

def notify_slack(summary):
    logger.info("Stub: Sending Slack notification (not implemented).")

def notify_webhook(summary):
    logger.info("Stub: Sending webhook notification (not implemented).")

def main():
    parser = argparse.ArgumentParser(
        description="Generate and report on backup validation status."
    )
    parser.add_argument(
        "--output", choices=["json", "md", "both"], default="both",
        help="Output format for the report (default: both)"
    )
    parser.add_argument(
        "--notify", choices=["none", "email", "slack", "webhook", "all"], default="none",
        help="Send notifications via specified channel(s) (default: none)"
    )
    parser.add_argument(
        "--verification-report", type=str, default=None,
        help="Path to a specific backup verification report (optional)"
    )
    parser.add_argument(
        "--restore-report", type=str, default=None,
        help="Path to a specific restore report (optional)"
    )
    args = parser.parse_args()

    logger.info("Starting validation report generation.")

    # Find latest reports if not specified
    verification_report_path = (
        args.verification_report or
        find_latest_report("reports/backup_verification", "verification")
    )
    restore_report_path = (
        args.restore_report or
        find_latest_report("reports/database_restore", "restore")
    )

    verification = load_json_report(verification_report_path) if verification_report_path else None
    restore = load_json_report(restore_report_path) if restore_report_path else None

    summary = summarize_reports(verification, restore)

    # Output report(s)
    output_dir = "reports/validation"
    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    output_files = []

    if args.output in ("json", "both"):
        json_path = os.path.join(output_dir, f"validation_report_{timestamp}.json")
        try:
            with open(json_path, "w") as f:
                json.dump(summary, f, indent=2)
            logger.info(f"Validation report written (JSON): {json_path}")
            output_files.append(json_path)
        except Exception as e:
            logger.error(f"Failed to write JSON report: {e}")

    if args.output in ("md", "both"):
        md_path = os.path.join(output_dir, f"validation_report_{timestamp}.md")
        try:
            with open(md_path, "w") as f:
                f.write(format_markdown(summary))
            logger.info(f"Validation report written (Markdown): {md_path}")
            output_files.append(md_path)
        except Exception as e:
            logger.error(f"Failed to write Markdown report: {e}")

    # Notifications
    if args.notify != "none":
        if args.notify in ("email", "all"):
            notify_email(summary)
        if args.notify in ("slack", "all"):
            notify_slack(summary)
        if args.notify in ("webhook", "all"):
            notify_webhook(summary)
        logger.info(f"Notifications sent via: {args.notify}")

    logger.info("Validation report generation completed.")
    print("Validation report(s) generated:")
    for f in output_files:
        print(f" - {f}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logger.error(f"Fatal error in validation_report.py: {e}")
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)