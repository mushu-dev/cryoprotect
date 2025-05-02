import os
import smtplib
import json
import requests
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart

def send_email(subject, body, config=None):
    """
    Send an email notification using SMTP.
    Config can be provided as a dict or read from environment variables.
    """
    if config is None:
        config = {
            "smtp_server": os.getenv("NOTIFY_SMTP_SERVER"),
            "smtp_port": int(os.getenv("NOTIFY_SMTP_PORT", "587")),
            "smtp_user": os.getenv("NOTIFY_SMTP_USER"),
            "smtp_password": os.getenv("NOTIFY_SMTP_PASSWORD"),
            "from_addr": os.getenv("NOTIFY_FROM_ADDR"),
            "to_addrs": os.getenv("NOTIFY_TO_ADDRS"),
            "use_tls": os.getenv("NOTIFY_SMTP_USE_TLS", "true").lower() == "true",
        }
    to_addrs = [a.strip() for a in config["to_addrs"].split(",") if a.strip()]
    msg = MIMEMultipart()
    msg["From"] = config["from_addr"]
    msg["To"] = ", ".join(to_addrs)
    msg["Subject"] = subject
    msg.attach(MIMEText(body, "plain"))

    try:
        server = smtplib.SMTP(config["smtp_server"], config["smtp_port"], timeout=10)
        if config["use_tls"]:
            server.starttls()
        if config["smtp_user"]:
            server.login(config["smtp_user"], config["smtp_password"])
        server.sendmail(config["from_addr"], to_addrs, msg.as_string())
        server.quit()
        return True
    except Exception as e:
        print(f"Failed to send email notification: {e}")
        return False

def send_webhook(payload, config=None):
    """
    Send a webhook notification (POST request with JSON payload).
    Config can be provided as a dict or read from environment variables.
    """
    if config is None:
        config = {
            "webhook_url": os.getenv("NOTIFY_WEBHOOK_URL"),
            "webhook_headers": os.getenv("NOTIFY_WEBHOOK_HEADERS", ""),
        }
    url = config["webhook_url"]
    headers = {}
    if config["webhook_headers"]:
        try:
            headers = json.loads(config["webhook_headers"])
        except Exception:
            pass
    try:
        resp = requests.post(url, json=payload, headers=headers, timeout=10)
        return resp.status_code in (200, 201, 204)
    except Exception as e:
        print(f"Failed to send webhook notification: {e}")
        return False

def notify_failure(job_name, exit_code, log_excerpt=None, extra_message=None):
    """
    Notify maintainers of a job failure via email and/or webhook.
    Reads config from environment variables.
    """
    subject = f"[CryoProtect] Scheduled Job Failed: {job_name} (exit code {exit_code})"
    body = f"Scheduled job '{job_name}' failed with exit code {exit_code}."
    if extra_message:
        body += f"\n\n{extra_message}"
    if log_excerpt:
        body += f"\n\nLog excerpt:\n{log_excerpt}"

    # Email notification
    if os.getenv("NOTIFY_TO_ADDRS"):
        send_email(subject, body)

    # Webhook notification
    if os.getenv("NOTIFY_WEBHOOK_URL"):
        payload = {
            "job": job_name,
            "status": "failure",
            "exit_code": exit_code,
            "message": extra_message,
            "log_excerpt": log_excerpt,
        }
        send_webhook(payload)