#!/usr/bin/env python3
"""
Monitoring Integrations Module for CryoProtect

This module provides integrations between the unified monitoring system and external services
for alerting and notification purposes, including:
1. Email notifications
2. Slack webhook integration
3. PagerDuty incident creation
4. SMS notifications (via Twilio)
5. Custom webhook support

Usage:
    from monitoring_integrations import (
        EmailNotifier, 
        SlackNotifier,
        PagerDutyNotifier,
        TwilioSMSNotifier,
        WebhookNotifier
    )
    
    # Initialize a notifier
    slack = SlackNotifier(webhook_url="https://hooks.slack.com/services/...")
    
    # Register with monitoring service
    monitor.add_alert_handler(slack.handle_alert)
"""

import os
import json
import logging
import requests
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from typing import Dict, Any, Optional, List, Callable
from datetime import datetime
from urllib.parse import urljoin

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('monitoring_integrations')

class BaseNotifier:
    """Base class for notification services."""
    
    def __init__(self, 
                 name: str,
                 min_severity: str = "warning",
                 throttle_seconds: int = 300,
                 include_context: bool = True,
                 custom_formatter: Optional[Callable] = None):
        """
        Initialize the base notifier.
        
        Args:
            name: Name of the notifier for identification
            min_severity: Minimum severity level to send notifications for
            throttle_seconds: Minimum time between identical alerts in seconds
            include_context: Whether to include context data in notifications
            custom_formatter: Optional custom formatter function
        """
        self.name = name
        self.min_severity = min_severity.lower()
        self.throttle_seconds = throttle_seconds
        self.include_context = include_context
        self.custom_formatter = custom_formatter
        
        # Track sent alerts to prevent duplicates
        self.last_sent_alerts = {}
        
        logger.info(f"Initialized {name} notifier")
    
    def handle_alert(self, alert: Dict[str, Any]) -> bool:
        """
        Handle an alert from the monitoring system.
        
        Args:
            alert: Alert data dictionary
            
        Returns:
            bool: True if notification was sent, False otherwise
        """
        # Check severity
        severity = alert.get('severity', 'error').lower()
        if not self._should_send_by_severity(severity):
            return False
        
        # Check throttling
        alert_key = f"{alert.get('alert_type')}:{alert.get('source')}"
        if self._is_throttled(alert_key):
            return False
        
        # Format the alert message
        message = self._format_alert(alert)
        
        # Send the notification
        success = self._send_notification(message, alert)
        
        # Record successful send for throttling
        if success:
            self.last_sent_alerts[alert_key] = datetime.now().timestamp()
            logger.info(f"{self.name} sent notification for alert: {alert_key}")
        else:
            logger.error(f"{self.name} failed to send notification for alert: {alert_key}")
        
        return success
    
    def _should_send_by_severity(self, severity: str) -> bool:
        """
        Check if alert should be sent based on severity.
        
        Args:
            severity: Alert severity
            
        Returns:
            bool: True if alert should be sent, False otherwise
        """
        severity_levels = {
            "debug": 0,
            "info": 1,
            "warning": 2,
            "error": 3,
            "critical": 4
        }
        
        alert_level = severity_levels.get(severity.lower(), 2)
        min_level = severity_levels.get(self.min_severity, 2)
        
        return alert_level >= min_level
    
    def _is_throttled(self, alert_key: str) -> bool:
        """
        Check if alert should be throttled.
        
        Args:
            alert_key: Alert identification key
            
        Returns:
            bool: True if alert should be throttled, False otherwise
        """
        if alert_key in self.last_sent_alerts:
            last_sent = self.last_sent_alerts[alert_key]
            now = datetime.now().timestamp()
            
            if now - last_sent < self.throttle_seconds:
                logger.debug(f"Alert {alert_key} throttled (sent {now - last_sent:.1f}s ago)")
                return True
        
        return False
    
    def _format_alert(self, alert: Dict[str, Any]) -> str:
        """
        Format alert data into a message.
        
        Args:
            alert: Alert data dictionary
            
        Returns:
            str: Formatted message
        """
        # Use custom formatter if provided
        if self.custom_formatter:
            return self.custom_formatter(alert)
        
        # Default formatting
        timestamp = alert.get('formatted_time', datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        severity = alert.get('severity', 'unknown').upper()
        source = alert.get('source', 'unknown')
        message = alert.get('message', 'No message provided')
        
        formatted = f"[{timestamp}] {severity} alert from {source}: {message}"
        
        # Add context if enabled
        if self.include_context and 'context' in alert and alert['context']:
            context_str = json.dumps(alert['context'], indent=2)
            formatted += f"\nContext: {context_str}"
        
        return formatted
    
    def _send_notification(self, message: str, alert: Dict[str, Any]) -> bool:
        """
        Send a notification.
        
        This method should be implemented by subclasses.
        
        Args:
            message: Formatted message
            alert: Original alert data
            
        Returns:
            bool: True if notification was sent successfully, False otherwise
        """
        raise NotImplementedError("Subclasses must implement _send_notification")


class EmailNotifier(BaseNotifier):
    """Email notification service."""
    
    def __init__(self, 
                 smtp_server: str,
                 smtp_port: int,
                 sender_email: str,
                 recipient_emails: List[str],
                 username: Optional[str] = None,
                 password: Optional[str] = None,
                 use_tls: bool = True,
                 subject_prefix: str = "[ALERT] CryoProtect",
                 **kwargs):
        """
        Initialize email notifier.
        
        Args:
            smtp_server: SMTP server address
            smtp_port: SMTP server port
            sender_email: Email address to send from
            recipient_emails: List of email addresses to send to
            username: SMTP authentication username (if different from sender_email)
            password: SMTP authentication password
            use_tls: Whether to use TLS for SMTP connection
            subject_prefix: Prefix for email subjects
            **kwargs: Additional arguments for BaseNotifier
        """
        super().__init__(name="Email", **kwargs)
        self.smtp_server = smtp_server
        self.smtp_port = smtp_port
        self.sender_email = sender_email
        self.recipient_emails = recipient_emails
        self.username = username or sender_email
        self.password = password
        self.use_tls = use_tls
        self.subject_prefix = subject_prefix
        
        logger.info(f"Configured Email notifier with server {smtp_server}:{smtp_port}")
    
    def _send_notification(self, message: str, alert: Dict[str, Any]) -> bool:
        """
        Send an email notification.
        
        Args:
            message: Formatted message
            alert: Original alert data
            
        Returns:
            bool: True if notification was sent successfully, False otherwise
        """
        try:
            # Create message
            email = MIMEMultipart("alternative")
            
            # Create subject with alert info
            alert_type = alert.get('alert_type', 'unknown')
            source = alert.get('source', 'unknown')
            severity = alert.get('severity', 'unknown').upper()
            
            subject = f"{self.subject_prefix} {severity} - {source} - {alert_type}"
            email["Subject"] = subject
            email["From"] = self.sender_email
            email["To"] = ", ".join(self.recipient_emails)
            
            # Add plain text and HTML parts
            text_part = MIMEText(message, "plain")
            email.attach(text_part)
            
            # Simple HTML formatting
            html_message = message.replace("\n", "<br>")
            html_part = MIMEText(f"<html><body><pre>{html_message}</pre></body></html>", "html")
            email.attach(html_part)
            
            # Connect to server and send
            with smtplib.SMTP(self.smtp_server, self.smtp_port) as server:
                if self.use_tls:
                    server.starttls()
                
                if self.username and self.password:
                    server.login(self.username, self.password)
                
                server.sendmail(self.sender_email, self.recipient_emails, email.as_string())
            
            return True
        except Exception as e:
            logger.error(f"Failed to send email notification: {str(e)}")
            return False


class SlackNotifier(BaseNotifier):
    """Slack webhook notification service."""
    
    def __init__(self, 
                 webhook_url: str,
                 channel: Optional[str] = None,
                 username: str = "CryoProtect Monitoring",
                 icon_emoji: str = ":warning:",
                 **kwargs):
        """
        Initialize Slack notifier.
        
        Args:
            webhook_url: Slack webhook URL
            channel: Optional channel override
            username: Bot username to display
            icon_emoji: Emoji icon for the bot
            **kwargs: Additional arguments for BaseNotifier
        """
        super().__init__(name="Slack", **kwargs)
        self.webhook_url = webhook_url
        self.channel = channel
        self.username = username
        self.icon_emoji = icon_emoji
        
        logger.info(f"Configured Slack notifier with webhook URL")
    
    def _send_notification(self, message: str, alert: Dict[str, Any]) -> bool:
        """
        Send a Slack notification.
        
        Args:
            message: Formatted message
            alert: Original alert data
            
        Returns:
            bool: True if notification was sent successfully, False otherwise
        """
        try:
            # Determine color based on severity
            severity = alert.get('severity', 'error').lower()
            color = self._get_color_for_severity(severity)
            
            # Create payload
            payload = {
                "username": self.username,
                "icon_emoji": self.icon_emoji,
                "attachments": [
                    {
                        "color": color,
                        "title": f"Alert: {alert.get('alert_type', 'Unknown Alert')}",
                        "text": message,
                        "fields": []
                    }
                ]
            }
            
            # Add channel if specified
            if self.channel:
                payload["channel"] = self.channel
            
            # Add fields for important alert data
            source = alert.get('source')
            if source:
                payload["attachments"][0]["fields"].append({
                    "title": "Source",
                    "value": source,
                    "short": True
                })
            
            severity_display = alert.get('severity', 'unknown').upper()
            payload["attachments"][0]["fields"].append({
                "title": "Severity",
                "value": severity_display,
                "short": True
            })
            
            # Add timestamp
            payload["attachments"][0]["ts"] = alert.get('timestamp', datetime.now().timestamp())
            
            # Send to webhook
            response = requests.post(
                self.webhook_url,
                json=payload,
                headers={"Content-Type": "application/json"},
                timeout=5
            )
            
            if response.status_code != 200:
                logger.error(f"Slack API error: {response.status_code} - {response.text}")
                return False
            
            return True
        except Exception as e:
            logger.error(f"Failed to send Slack notification: {str(e)}")
            return False
    
    def _get_color_for_severity(self, severity: str) -> str:
        """
        Get Slack attachment color for severity level.
        
        Args:
            severity: Severity level
            
        Returns:
            str: Hex color code
        """
        colors = {
            "debug": "#999999",    # Gray
            "info": "#3498db",     # Blue
            "warning": "#f39c12",  # Orange
            "error": "#e74c3c",    # Red
            "critical": "#7b241c"  # Dark Red
        }
        
        return colors.get(severity.lower(), "#f39c12")


class PagerDutyNotifier(BaseNotifier):
    """PagerDuty notification service."""
    
    def __init__(self, 
                 integration_key: str,
                 severity_mapping: Optional[Dict[str, str]] = None,
                 component: str = "CryoProtect",
                 group: str = "Monitoring",
                 class_type: str = "application",
                 **kwargs):
        """
        Initialize PagerDuty notifier.
        
        Args:
            integration_key: PagerDuty integration/routing key
            severity_mapping: Optional mapping from internal to PagerDuty severities
            component: Component name
            group: Group name
            class_type: Class value for the event
            **kwargs: Additional arguments for BaseNotifier
        """
        super().__init__(name="PagerDuty", **kwargs)
        self.integration_key = integration_key
        self.component = component
        self.group = group
        self.class_type = class_type
        
        # PagerDuty severities: info, warning, error, critical
        self.severity_mapping = severity_mapping or {
            "debug": "info",
            "info": "info",
            "warning": "warning",
            "error": "error",
            "critical": "critical"
        }
        
        logger.info(f"Configured PagerDuty notifier with integration key")
    
    def _send_notification(self, message: str, alert: Dict[str, Any]) -> bool:
        """
        Send a PagerDuty notification.
        
        Args:
            message: Formatted message
            alert: Original alert data
            
        Returns:
            bool: True if notification was sent successfully, False otherwise
        """
        try:
            # Map internal severity to PagerDuty severity
            internal_severity = alert.get('severity', 'error').lower()
            pd_severity = self.severity_mapping.get(internal_severity, "error")
            
            # Generate a deduplication key based on alert data
            alert_type = alert.get('alert_type', 'unknown')
            source = alert.get('source', 'unknown')
            dedup_key = f"{self.component}:{source}:{alert_type}"
            
            # Create event payload
            payload = {
                "routing_key": self.integration_key,
                "event_action": "trigger",
                "dedup_key": dedup_key,
                "payload": {
                    "summary": alert.get('message', 'No message provided'),
                    "source": source,
                    "severity": pd_severity,
                    "component": self.component,
                    "group": self.group,
                    "class": self.class_type,
                    "custom_details": {
                        "alert_type": alert_type,
                        "message": message
                    }
                }
            }
            
            # Add context if available
            if self.include_context and 'context' in alert and alert['context']:
                payload["payload"]["custom_details"]["context"] = alert["context"]
            
            # Send to PagerDuty Events API
            response = requests.post(
                "https://events.pagerduty.com/v2/enqueue",
                json=payload,
                headers={"Content-Type": "application/json"},
                timeout=5
            )
            
            if response.status_code != 202:
                logger.error(f"PagerDuty API error: {response.status_code} - {response.text}")
                return False
            
            return True
        except Exception as e:
            logger.error(f"Failed to send PagerDuty notification: {str(e)}")
            return False


class TwilioSMSNotifier(BaseNotifier):
    """Twilio SMS notification service."""
    
    def __init__(self, 
                 account_sid: str,
                 auth_token: str,
                 from_number: str,
                 to_numbers: List[str],
                 **kwargs):
        """
        Initialize Twilio SMS notifier.
        
        Args:
            account_sid: Twilio account SID
            auth_token: Twilio auth token
            from_number: Twilio phone number to send from
            to_numbers: List of phone numbers to send to
            **kwargs: Additional arguments for BaseNotifier
        """
        super().__init__(name="Twilio SMS", **kwargs)
        self.account_sid = account_sid
        self.auth_token = auth_token
        self.from_number = from_number
        self.to_numbers = to_numbers
        
        # Import twilio here to make it optional
        try:
            from twilio.rest import Client
            self.client = Client(account_sid, auth_token)
            self.twilio_available = True
        except ImportError:
            logger.warning("Twilio package not available. SMS notifications will be simulated.")
            self.twilio_available = False
            self.client = None
        
        logger.info(f"Configured Twilio SMS notifier with account {account_sid}")
    
    def _format_alert(self, alert: Dict[str, Any]) -> str:
        """
        Format alert data into a message suitable for SMS.
        
        Args:
            alert: Alert data dictionary
            
        Returns:
            str: Formatted message
        """
        # Use custom formatter if provided
        if self.custom_formatter:
            return self.custom_formatter(alert)
        
        # SMS-specific formatting (shorter than default)
        severity = alert.get('severity', 'unknown').upper()
        source = alert.get('source', 'unknown')
        message = alert.get('message', 'No message provided')
        
        return f"CryoProtect {severity} ALERT - {source}: {message}"
    
    def _send_notification(self, message: str, alert: Dict[str, Any]) -> bool:
        """
        Send an SMS notification.
        
        Args:
            message: Formatted message
            alert: Original alert data
            
        Returns:
            bool: True if notification was sent successfully, False otherwise
        """
        try:
            if not self.twilio_available:
                logger.info(f"SIMULATED SMS to {self.to_numbers}: {message}")
                return True
            
            # Send to each recipient
            successful_sends = 0
            for to_number in self.to_numbers:
                sms = self.client.messages.create(
                    body=message,
                    from_=self.from_number,
                    to=to_number
                )
                
                if sms.sid:
                    successful_sends += 1
                    logger.debug(f"Sent SMS to {to_number}, SID: {sms.sid}")
            
            # Return success if at least one message was sent
            return successful_sends > 0
            
        except Exception as e:
            logger.error(f"Failed to send SMS notification: {str(e)}")
            return False


class WebhookNotifier(BaseNotifier):
    """Generic webhook notification service."""
    
    def __init__(self, 
                 webhook_url: str,
                 http_method: str = "POST",
                 headers: Optional[Dict[str, str]] = None,
                 payload_template: Optional[Dict[str, Any]] = None,
                 verify_ssl: bool = True,
                 timeout: int = 5,
                 **kwargs):
        """
        Initialize webhook notifier.
        
        Args:
            webhook_url: Webhook URL
            http_method: HTTP method to use (GET, POST, etc.)
            headers: HTTP headers to include
            payload_template: Optional template for JSON payload
            verify_ssl: Whether to verify SSL certificates
            timeout: Request timeout in seconds
            **kwargs: Additional arguments for BaseNotifier
        """
        super().__init__(name="Webhook", **kwargs)
        self.webhook_url = webhook_url
        self.http_method = http_method.upper()
        self.headers = headers or {"Content-Type": "application/json"}
        self.payload_template = payload_template
        self.verify_ssl = verify_ssl
        self.timeout = timeout
        
        logger.info(f"Configured Webhook notifier with URL {webhook_url}")
    
    def _send_notification(self, message: str, alert: Dict[str, Any]) -> bool:
        """
        Send a webhook notification.
        
        Args:
            message: Formatted message
            alert: Original alert data
            
        Returns:
            bool: True if notification was sent successfully, False otherwise
        """
        try:
            # Create payload
            if self.payload_template:
                # Use template and fill in alert data
                payload = self._fill_payload_template(message, alert)
            else:
                # Use alert data directly
                payload = {
                    "message": message,
                    "timestamp": alert.get('timestamp', datetime.now().timestamp()),
                    "formatted_time": alert.get('formatted_time', datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
                    "alert_type": alert.get('alert_type', 'unknown'),
                    "source": alert.get('source', 'unknown'),
                    "severity": alert.get('severity', 'unknown')
                }
                
                # Add context if enabled
                if self.include_context and 'context' in alert and alert['context']:
                    payload["context"] = alert["context"]
            
            # Send request
            if self.http_method == "GET":
                response = requests.get(
                    self.webhook_url,
                    params=payload,
                    headers=self.headers,
                    verify=self.verify_ssl,
                    timeout=self.timeout
                )
            else:  # POST, PUT, etc.
                request_func = getattr(requests, self.http_method.lower())
                response = request_func(
                    self.webhook_url,
                    json=payload,
                    headers=self.headers,
                    verify=self.verify_ssl,
                    timeout=self.timeout
                )
            
            if response.status_code >= 400:
                logger.error(f"Webhook API error: {response.status_code} - {response.text}")
                return False
            
            return True
        except Exception as e:
            logger.error(f"Failed to send webhook notification: {str(e)}")
            return False
    
    def _fill_payload_template(self, message: str, alert: Dict[str, Any]) -> Dict[str, Any]:
        """
        Fill payload template with alert data.
        
        Args:
            message: Formatted message
            alert: Alert data
            
        Returns:
            Dict[str, Any]: Filled payload
        """
        import copy
        import re
        
        # Deep copy template
        payload = copy.deepcopy(self.payload_template)
        
        # Add standard fields to alert data for template
        alert_data = dict(alert)
        alert_data["formatted_message"] = message
        
        # Define a recursive function to replace placeholders
        def process_dict(d):
            if isinstance(d, dict):
                for k, v in d.items():
                    if isinstance(v, (dict, list)):
                        process_dict(v)
                    elif isinstance(v, str):
                        d[k] = self._replace_placeholders(v, alert_data)
            elif isinstance(d, list):
                for i, item in enumerate(d):
                    if isinstance(item, (dict, list)):
                        process_dict(item)
                    elif isinstance(item, str):
                        d[i] = self._replace_placeholders(item, alert_data)
        
        # Process the template
        process_dict(payload)
        return payload
    
    def _replace_placeholders(self, text: str, data: Dict[str, Any]) -> str:
        """
        Replace placeholders in text with data values.
        
        Args:
            text: Text with placeholders
            data: Data dictionary
            
        Returns:
            str: Text with placeholders replaced
        """
        # Replace placeholders in the form {{key}}
        for key, value in data.items():
            if isinstance(value, (str, int, float, bool)):
                text = text.replace(f"{{{{{key}}}}}", str(value))
        
        return text


# Example usage and utility functions

def create_email_notifier(
    smtp_server: str,
    smtp_port: int,
    sender_email: str,
    recipient_emails: List[str],
    username: Optional[str] = None,
    password: Optional[str] = None,
    **kwargs
) -> EmailNotifier:
    """
    Create an email notifier from configuration.
    
    Args:
        smtp_server: SMTP server address
        smtp_port: SMTP server port
        sender_email: Email address to send from
        recipient_emails: List of email addresses to send to
        username: SMTP authentication username
        password: SMTP authentication password
        **kwargs: Additional arguments for EmailNotifier
        
    Returns:
        EmailNotifier: Configured email notifier
    """
    return EmailNotifier(
        smtp_server=smtp_server,
        smtp_port=smtp_port,
        sender_email=sender_email,
        recipient_emails=recipient_emails if isinstance(recipient_emails, list) else [recipient_emails],
        username=username,
        password=password,
        **kwargs
    )


def create_slack_notifier(webhook_url: str, **kwargs) -> SlackNotifier:
    """
    Create a Slack notifier from configuration.
    
    Args:
        webhook_url: Slack webhook URL
        **kwargs: Additional arguments for SlackNotifier
        
    Returns:
        SlackNotifier: Configured Slack notifier
    """
    return SlackNotifier(
        webhook_url=webhook_url,
        **kwargs
    )


def create_pagerduty_notifier(integration_key: str, **kwargs) -> PagerDutyNotifier:
    """
    Create a PagerDuty notifier from configuration.
    
    Args:
        integration_key: PagerDuty integration/routing key
        **kwargs: Additional arguments for PagerDutyNotifier
        
    Returns:
        PagerDutyNotifier: Configured PagerDuty notifier
    """
    return PagerDutyNotifier(
        integration_key=integration_key,
        **kwargs
    )


def create_sms_notifier(
    account_sid: str,
    auth_token: str,
    from_number: str,
    to_numbers: List[str],
    **kwargs
) -> TwilioSMSNotifier:
    """
    Create an SMS notifier from configuration.
    
    Args:
        account_sid: Twilio account SID
        auth_token: Twilio auth token
        from_number: Twilio phone number to send from
        to_numbers: List of phone numbers to send to
        **kwargs: Additional arguments for TwilioSMSNotifier
        
    Returns:
        TwilioSMSNotifier: Configured SMS notifier
    """
    return TwilioSMSNotifier(
        account_sid=account_sid,
        auth_token=auth_token,
        from_number=from_number,
        to_numbers=to_numbers if isinstance(to_numbers, list) else [to_numbers],
        **kwargs
    )


def create_webhook_notifier(webhook_url: str, **kwargs) -> WebhookNotifier:
    """
    Create a webhook notifier from configuration.
    
    Args:
        webhook_url: Webhook URL
        **kwargs: Additional arguments for WebhookNotifier
        
    Returns:
        WebhookNotifier: Configured webhook notifier
    """
    return WebhookNotifier(
        webhook_url=webhook_url,
        **kwargs
    )


def create_notifiers_from_config(config: Dict[str, Any]) -> List[BaseNotifier]:
    """
    Create notifiers from configuration dictionary.
    
    Args:
        config: Configuration dictionary
        
    Returns:
        List[BaseNotifier]: List of configured notifiers
    """
    notifiers = []
    
    # Create email notifier if configured
    if 'email' in config and config['email'].get('enabled', False):
        email_config = config['email']
        notifiers.append(create_email_notifier(**email_config))
    
    # Create Slack notifier if configured
    if 'slack' in config and config['slack'].get('enabled', False):
        slack_config = config['slack']
        notifiers.append(create_slack_notifier(**slack_config))
    
    # Create PagerDuty notifier if configured
    if 'pagerduty' in config and config['pagerduty'].get('enabled', False):
        pagerduty_config = config['pagerduty']
        notifiers.append(create_pagerduty_notifier(**pagerduty_config))
    
    # Create SMS notifier if configured
    if 'sms' in config and config['sms'].get('enabled', False):
        sms_config = config['sms']
        notifiers.append(create_sms_notifier(**sms_config))
    
    # Create webhook notifiers if configured
    if 'webhooks' in config:
        for webhook_config in config['webhooks']:
            if webhook_config.get('enabled', False):
                notifiers.append(create_webhook_notifier(**webhook_config))
    
    return notifiers


# Example configuration
EXAMPLE_CONFIG = {
    "email": {
        "enabled": True,
        "smtp_server": "smtp.example.com",
        "smtp_port": 587,
        "sender_email": "monitoring@example.com",
        "recipient_emails": ["admin@example.com"],
        "username": "monitoring@example.com",
        "password": "password123",
        "use_tls": True,
        "min_severity": "warning"
    },
    "slack": {
        "enabled": True,
        "webhook_url": "https://hooks.slack.com/services/XXX/YYY/ZZZ",
        "channel": "#monitoring",
        "min_severity": "warning"
    },
    "pagerduty": {
        "enabled": True,
        "integration_key": "abcdef123456",
        "min_severity": "error"
    },
    "sms": {
        "enabled": False,
        "account_sid": "ACXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
        "auth_token": "your_auth_token",
        "from_number": "+1234567890",
        "to_numbers": ["+1234567890"],
        "min_severity": "critical"
    },
    "webhooks": [
        {
            "enabled": True,
            "webhook_url": "https://example.com/webhook",
            "http_method": "POST",
            "min_severity": "warning"
        }
    ]
}


# If this module is executed directly, run an example
if __name__ == "__main__":
    import sys
    
    # Create notifiers from example config
    notifiers = create_notifiers_from_config(EXAMPLE_CONFIG)
    
    # Print configured notifiers
    print(f"Configured {len(notifiers)} notifiers:")
    for notifier in notifiers:
        print(f"- {notifier.name} (min severity: {notifier.min_severity})")
    
    # Create a test alert
    test_alert = {
        "alert_type": "example_alert",
        "message": "This is a test alert message",
        "source": "monitoring_integrations.py",
        "severity": "warning",
        "timestamp": datetime.now().timestamp(),
        "formatted_time": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "context": {
            "test": True,
            "value": 42
        }
    }
    
    # Send test alert to all notifiers
    print("\nSending test alert...")
    for notifier in notifiers:
        # Skip actual sending in test mode
        if len(sys.argv) > 1 and sys.argv[1] == "--test":
            print(f"Would send to {notifier.name} notifier")
            formatted = notifier._format_alert(test_alert)
            print(f"Formatted message: {formatted}")
        else:
            result = notifier.handle_alert(test_alert)
            print(f"Sent to {notifier.name} notifier: {'success' if result else 'failed'}")