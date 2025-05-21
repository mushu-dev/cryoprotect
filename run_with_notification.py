import os
import sys
import subprocess
from notify import notify_failure

def get_log_excerpt(log_file, num_lines=40):
    """Return the last num_lines from the log file, if it exists."""
    if not log_file or not os.path.exists(log_file):
        return None
    try:
        with open(log_file, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()
        return "".join(lines[-num_lines:])
    except Exception:
        return None

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Run a command and notify maintainers on failure (nonzero exit)."
    )
    parser.add_argument("job_name", help="Name of the job (for notification)")
    parser.add_argument("command", nargs=argparse.REMAINDER, help="Command to run (e.g., python update_pubchem_data.py)")
    parser.add_argument("--log-file", help="Path to log file to include excerpt in notification", default=os.getenv("LOG_FILE"))
    parser.add_argument("--lines", type=int, default=40, help="Number of log lines to include in notification")
    args = parser.parse_args()

    if not args.command:
        print("No command specified to run.")
        sys.exit(2)

    # Run the command
    print(f"Running: {' '.join(args.command)}")
    result = subprocess.run(args.command)
    exit_code = result.returncode

    if exit_code != 0:
        log_excerpt = get_log_excerpt(args.log_file, args.lines) if args.log_file else None
        notify_failure(
            job_name=args.job_name,
            exit_code=exit_code,
            log_excerpt=log_excerpt,
            extra_message=f"Command: {' '.join(args.command)}"
        )
        print(f"Job failed (exit code {exit_code}). Notification sent.")
    else:
        print("Job succeeded.")

    sys.exit(exit_code)

if __name__ == "__main__":
    main()