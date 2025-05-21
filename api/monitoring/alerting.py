#!/usr/bin/env python3
"""
CryoProtect v2 - Alerting System

This module provides a centralized alerting system:
- Multi-channel alert notifications (console, email, Slack, etc.)
- Customizable alert thresholds and rules
- Alert aggregation and deduplication
- Alert status tracking and resolution
- On-call rotation and escalation policies
- Integration with error tracking and performance monitoring

Usage:
    from api.monitoring.alerting import (
        AlertManager,
        Alert,
        AlertSeverity,
        AlertChannel,
        send_alert,
        register_alert_handler
    )
"""

import os
import time
import uuid
import json
import enum
import smtplib
import logging
import threading
import traceback
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from typing import Dict, List, Any, Set, Optional, Union, Callable, Type, TypeVar, Tuple

# Import structured logging
from logging_enhanced import log_with_context, get_logger

# Type variables
T = TypeVar('T')  # Return type for functions

# Set up logger
logger = get_logger(__name__)


class AlertSeverity(enum.Enum):
    """Alert severity levels."""
    INFO = 0        # Informational, no action needed
    WARNING = 1     # Warning, might need attention
    ERROR = 2       # Error, needs attention
    CRITICAL = 3    # Critical, immediate attention required


class AlertChannel(enum.Enum):
    """Alert notification channels."""
    CONSOLE = "console"  # Log to console
    EMAIL = "email"      # Send email notification
    SLACK = "slack"      # Send Slack notification
    WEBHOOK = "webhook"  # Send webhook notification
    SMS = "sms"          # Send SMS notification
    CUSTOM = "custom"    # Custom notification channel


class AlertStatus(enum.Enum):
    """Alert status."""
    OPEN = "open"                # New alert, not acknowledged
    ACKNOWLEDGED = "acknowledged"  # Alert acknowledged
    RESOLVED = "resolved"        # Alert resolved
    MUTED = "muted"              # Alert muted (ignored)


class Alert:
    """
    Represents an alert notification.
    """
    
    def __init__(
        self,
        title: str,
        message: str,
        severity: AlertSeverity,
        source: str,
        context: Optional[Dict[str, Any]] = None,
        fingerprint: Optional[str] = None,
        alert_id: Optional[str] = None
    ):
        """
        Initialize an alert.
        
        Args:
            title: Alert title
            message: Alert message
            severity: Alert severity
            source: Alert source (e.g., "error_tracking", "performance")
            context: Additional context
            fingerprint: Alert fingerprint for deduplication
            alert_id: Alert ID (generated if None)
        """
        self.title = title
        self.message = message
        self.severity = severity
        self.source = source
        self.context = context or {}
        self.alert_id = alert_id or str(uuid.uuid4())
        self.created_at = time.time()
        self.updated_at = self.created_at
        self.status = AlertStatus.OPEN
        self.acknowledged_at: Optional[float] = None
        self.resolved_at: Optional[float] = None
        self.acknowledged_by: Optional[str] = None
        self.resolved_by: Optional[str] = None
        self.resolution_message: Optional[str] = None
        self.notification_history: List[Dict[str, Any]] = []
        
        # Generate fingerprint if not provided
        if fingerprint is None:
            # Create a fingerprint based on title and source
            self.fingerprint = f"{source}:{title}"
        else:
            self.fingerprint = fingerprint
    
    def acknowledge(self, user: Optional[str] = None) -> None:
        """
        Acknowledge the alert.
        
        Args:
            user: User who acknowledged the alert
        """
        if self.status == AlertStatus.OPEN:
            self.status = AlertStatus.ACKNOWLEDGED
            self.acknowledged_at = time.time()
            self.updated_at = self.acknowledged_at
            self.acknowledged_by = user
    
    def resolve(self, user: Optional[str] = None, message: Optional[str] = None) -> None:
        """
        Resolve the alert.
        
        Args:
            user: User who resolved the alert
            message: Resolution message
        """
        if self.status != AlertStatus.RESOLVED:
            self.status = AlertStatus.RESOLVED
            self.resolved_at = time.time()
            self.updated_at = self.resolved_at
            self.resolved_by = user
            self.resolution_message = message
    
    def mute(self) -> None:
        """Mute the alert."""
        if self.status != AlertStatus.RESOLVED:
            self.status = AlertStatus.MUTED
            self.updated_at = time.time()
    
    def add_notification(self, channel: AlertChannel, sent_at: float, success: bool, details: Optional[str] = None) -> None:
        """
        Add a notification record to the alert history.
        
        Args:
            channel: Notification channel
            sent_at: Notification timestamp
            success: Whether the notification was sent successfully
            details: Optional details about the notification
        """
        self.notification_history.append({
            'channel': channel.value,
            'sent_at': sent_at,
            'success': success,
            'details': details
        })
    
    def to_dict(self) -> Dict[str, Any]:
        """
        Convert the alert to a dictionary.
        
        Returns:
            Dictionary representation of the alert
        """
        return {
            'alert_id': self.alert_id,
            'fingerprint': self.fingerprint,
            'title': self.title,
            'message': self.message,
            'severity': self.severity.name,
            'source': self.source,
            'context': self.context,
            'created_at': self.created_at,
            'updated_at': self.updated_at,
            'status': self.status.value,
            'acknowledged_at': self.acknowledged_at,
            'resolved_at': self.resolved_at,
            'acknowledged_by': self.acknowledged_by,
            'resolved_by': self.resolved_by,
            'resolution_message': self.resolution_message,
            'notification_history': self.notification_history
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'Alert':
        """
        Create an Alert from a dictionary.
        
        Args:
            data: Dictionary representation of the alert
            
        Returns:
            Alert instance
        """
        alert = cls(
            title=data['title'],
            message=data['message'],
            severity=AlertSeverity[data['severity']],
            source=data['source'],
            context=data['context'],
            fingerprint=data['fingerprint'],
            alert_id=data['alert_id']
        )
        
        alert.created_at = data['created_at']
        alert.updated_at = data['updated_at']
        alert.status = AlertStatus(data['status'])
        alert.acknowledged_at = data['acknowledged_at']
        alert.resolved_at = data['resolved_at']
        alert.acknowledged_by = data['acknowledged_by']
        alert.resolved_by = data['resolved_by']
        alert.resolution_message = data['resolution_message']
        alert.notification_history = data['notification_history']
        
        return alert
    
    def to_json(self) -> str:
        """
        Convert the alert to a JSON string.
        
        Returns:
            JSON string representation of the alert
        """
        return json.dumps(self.to_dict(), default=str)
    
    @classmethod
    def from_json(cls, json_str: str) -> 'Alert':
        """
        Create an Alert from a JSON string.
        
        Args:
            json_str: JSON string representation of the alert
            
        Returns:
            Alert instance
        """
        data = json.loads(json_str)
        return cls.from_dict(data)


class NotificationChannel:
    """
    Base class for notification channels.
    """
    
    def send(self, alert: Alert) -> Tuple[bool, Optional[str]]:
        """
        Send an alert notification.
        
        Args:
            alert: Alert to send
            
        Returns:
            Tuple of (success, details)
        """
        raise NotImplementedError("Subclasses must implement send()")


class ConsoleChannel(NotificationChannel):
    """
    Console notification channel.
    """
    
    def __init__(self, level: int = logging.INFO):
        """
        Initialize console channel.
        
        Args:
            level: Logging level
        """
        self.level = level
        self.logger = get_logger(__name__)
    
    def send(self, alert: Alert) -> Tuple[bool, Optional[str]]:
        """
        Send an alert notification to the console.
        
        Args:
            alert: Alert to send
            
        Returns:
            Tuple of (success, details)
        """
        # Map severity to log level
        severity_to_level = {
            AlertSeverity.INFO: logging.INFO,
            AlertSeverity.WARNING: logging.WARNING,
            AlertSeverity.ERROR: logging.ERROR,
            AlertSeverity.CRITICAL: logging.CRITICAL
        }
        
        log_level = severity_to_level.get(alert.severity, self.level)
        
        # Create log message
        message = f"ALERT ({alert.severity.name}): {alert.title} - {alert.message}"
        
        # Log the alert
        log_with_context(
            self.logger,
            logging.getLevelName(log_level).lower(),
            message,
            context={
                'alerting': {
                    'alert_id': alert.alert_id,
                    'fingerprint': alert.fingerprint,
                    'severity': alert.severity.name,
                    'source': alert.source
                },
                **alert.context
            }
        )
        
        return True, None


class EmailChannel(NotificationChannel):
    """
    Email notification channel.
    """
    
    def __init__(
        self,
        smtp_server: str,
        smtp_port: int,
        smtp_username: str,
        smtp_password: str,
        sender: str,
        recipients: List[str],
        use_tls: bool = True
    ):
        """
        Initialize email channel.
        
        Args:
            smtp_server: SMTP server hostname
            smtp_port: SMTP server port
            smtp_username: SMTP username
            smtp_password: SMTP password
            sender: Sender email address
            recipients: List of recipient email addresses
            use_tls: Whether to use TLS
        """
        self.smtp_server = smtp_server
        self.smtp_port = smtp_port
        self.smtp_username = smtp_username
        self.smtp_password = smtp_password
        self.sender = sender
        self.recipients = recipients
        self.use_tls = use_tls
    
    def send(self, alert: Alert) -> Tuple[bool, Optional[str]]:
        """
        Send an alert notification via email.
        
        Args:
            alert: Alert to send
            
        Returns:
            Tuple of (success, details)
        """
        try:
            # Create message
            msg = MIMEMultipart()
            msg['From'] = self.sender
            msg['To'] = ', '.join(self.recipients)
            msg['Subject'] = f"[{alert.severity.name}] {alert.title}"
            
            # Create email body
            body = f"""
            Alert: {alert.title}
            Severity: {alert.severity.name}
            Source: {alert.source}
            Time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(alert.created_at))}
            
            Message:
            {alert.message}
            
            Context:
            {json.dumps(alert.context, indent=2, default=str)}
            
            Alert ID: {alert.alert_id}
            """
            
            msg.attach(MIMEText(body, 'plain'))
            
            # Connect to SMTP server
            server = smtplib.SMTP(self.smtp_server, self.smtp_port)
            
            if self.use_tls:
                server.starttls()
            
            server.login(self.smtp_username, self.smtp_password)
            
            # Send email
            server.send_message(msg)
            server.quit()
            
            return True, None
        except Exception as e:
            return False, str(e)


class SlackChannel(NotificationChannel):
    """
    Slack notification channel.
    """
    
    def __init__(self, webhook_url: str):
        """
        Initialize Slack channel.
        
        Args:
            webhook_url: Slack webhook URL
        """
        self.webhook_url = webhook_url
    
    def send(self, alert: Alert) -> Tuple[bool, Optional[str]]:
        """
        Send an alert notification to Slack.
        
        Args:
            alert: Alert to send
            
        Returns:
            Tuple of (success, details)
        """
        try:
            # Map severity to color
            severity_to_color = {
                AlertSeverity.INFO: "#2196F3",      # Blue
                AlertSeverity.WARNING: "#FFC107",   # Amber
                AlertSeverity.ERROR: "#F44336",     # Red
                AlertSeverity.CRITICAL: "#9C27B0"   # Purple
            }
            
            color = severity_to_color.get(alert.severity, "#2196F3")
            
            # Format created_at timestamp
            timestamp = time.strftime(
                '%Y-%m-%d %H:%M:%S',
                time.localtime(alert.created_at)
            )
            
            # Create Slack message payload
            payload = {
                "attachments": [
                    {
                        "fallback": f"{alert.severity.name}: {alert.title}",
                        "color": color,
                        "title": f"{alert.severity.name}: {alert.title}",
                        "text": alert.message,
                        "fields": [
                            {
                                "title": "Source",
                                "value": alert.source,
                                "short": True
                            },
                            {
                                "title": "Alert ID",
                                "value": alert.alert_id,
                                "short": True
                            },
                            {
                                "title": "Time",
                                "value": timestamp,
                                "short": True
                            }
                        ],
                        "footer": "CryoProtect Alerting System",
                        "ts": int(alert.created_at)
                    }
                ]
            }
            
            # Add context as fields if not too large
            if alert.context and len(json.dumps(alert.context)) < 1000:
                for key, value in alert.context.items():
                    if isinstance(value, (str, int, float, bool)):
                        payload["attachments"][0]["fields"].append({
                            "title": key,
                            "value": str(value),
                            "short": True
                        })
            
            # Send to Slack
            import requests
            response = requests.post(
                self.webhook_url,
                json=payload,
                headers={"Content-Type": "application/json"}
            )
            
            if response.status_code == 200:
                return True, None
            else:
                return False, f"Slack API error: {response.status_code} - {response.text}"
        except Exception as e:
            return False, str(e)


class WebhookChannel(NotificationChannel):
    """
    Webhook notification channel.
    """
    
    def __init__(self, webhook_url: str, headers: Optional[Dict[str, str]] = None):
        """
        Initialize webhook channel.
        
        Args:
            webhook_url: Webhook URL
            headers: HTTP headers
        """
        self.webhook_url = webhook_url
        self.headers = headers or {"Content-Type": "application/json"}
    
    def send(self, alert: Alert) -> Tuple[bool, Optional[str]]:
        """
        Send an alert notification to a webhook.
        
        Args:
            alert: Alert to send
            
        Returns:
            Tuple of (success, details)
        """
        try:
            # Create webhook payload
            payload = alert.to_dict()
            
            # Send to webhook
            import requests
            response = requests.post(
                self.webhook_url,
                json=payload,
                headers=self.headers
            )
            
            if response.status_code >= 200 and response.status_code < 300:
                return True, None
            else:
                return False, f"Webhook error: {response.status_code} - {response.text}"
        except Exception as e:
            return False, str(e)


class CustomChannel(NotificationChannel):
    """
    Custom notification channel.
    """
    
    def __init__(self, handler: Callable[[Alert], Tuple[bool, Optional[str]]]):
        """
        Initialize custom channel.
        
        Args:
            handler: Function to handle the notification
        """
        self.handler = handler
    
    def send(self, alert: Alert) -> Tuple[bool, Optional[str]]:
        """
        Send an alert notification using the custom handler.
        
        Args:
            alert: Alert to send
            
        Returns:
            Tuple of (success, details)
        """
        try:
            return self.handler(alert)
        except Exception as e:
            return False, str(e)


class AlertManager:
    """
    Manages alerts and notification channels.
    """
    
    _instance = None
    
    @classmethod
    def get_instance(cls) -> 'AlertManager':
        """
        Get or create the singleton instance of AlertManager.
        
        Returns:
            AlertManager instance
        """
        if cls._instance is None:
            cls._instance = AlertManager()
        return cls._instance
    
    def __init__(self):
        """Initialize alert manager."""
        # Ensure singleton pattern
        if AlertManager._instance is not None:
            raise RuntimeError("AlertManager is a singleton. Use get_instance() instead.")
        
        AlertManager._instance = self
        
        # Notification channels
        self.channels: Dict[AlertChannel, NotificationChannel] = {}
        
        # Default channels by severity
        self.default_channels: Dict[AlertSeverity, List[AlertChannel]] = {
            AlertSeverity.INFO: [AlertChannel.CONSOLE],
            AlertSeverity.WARNING: [AlertChannel.CONSOLE],
            AlertSeverity.ERROR: [AlertChannel.CONSOLE],
            AlertSeverity.CRITICAL: [AlertChannel.CONSOLE]
        }
        
        # Alerts storage
        self.alerts: Dict[str, Alert] = {}  # alert_id -> Alert
        self.active_alerts: Dict[str, Alert] = {}  # fingerprint -> Alert (most recent)
        
        # Alert handlers
        self.alert_handlers: List[Callable[[Alert], None]] = []
        
        # Rate limiting for alerts
        self.cooldown_periods: Dict[AlertSeverity, float] = {
            AlertSeverity.INFO: 3600,      # 1 hour
            AlertSeverity.WARNING: 1800,   # 30 minutes
            AlertSeverity.ERROR: 900,      # 15 minutes
            AlertSeverity.CRITICAL: 300    # 5 minutes
        }
        self.last_alert_times: Dict[str, float] = {}  # fingerprint -> timestamp
        
        # Thread safety
        self.lock = threading.RLock()
        
        # Create logger
        self.logger = get_logger(__name__)
        
        # Initialize default channels
        self._init_default_channels()
    
    def _init_default_channels(self) -> None:
        """Initialize default notification channels."""
        # Add console channel
        self.register_channel(AlertChannel.CONSOLE, ConsoleChannel())
    
    def register_channel(self, channel_type: AlertChannel, channel: NotificationChannel) -> None:
        """
        Register a notification channel.
        
        Args:
            channel_type: Channel type
            channel: Channel instance
        """
        with self.lock:
            self.channels[channel_type] = channel
    
    def set_default_channels(self, severity: AlertSeverity, channels: List[AlertChannel]) -> None:
        """
        Set default channels for a severity level.
        
        Args:
            severity: Alert severity
            channels: List of channel types
        """
        with self.lock:
            self.default_channels[severity] = channels
    
    def set_cooldown_period(self, severity: AlertSeverity, seconds: float) -> None:
        """
        Set the cooldown period for a severity level.
        
        Args:
            severity: Alert severity
            seconds: Cooldown period in seconds
        """
        with self.lock:
            self.cooldown_periods[severity] = seconds
    
    def send_alert(
        self,
        title: str,
        message: str,
        severity: AlertSeverity,
        source: str,
        context: Optional[Dict[str, Any]] = None,
        fingerprint: Optional[str] = None,
        channels: Optional[List[AlertChannel]] = None
    ) -> Optional[Alert]:
        """
        Send an alert notification.
        
        Args:
            title: Alert title
            message: Alert message
            severity: Alert severity
            source: Alert source
            context: Additional context
            fingerprint: Alert fingerprint for deduplication
            channels: Channels to send the alert to (uses default channels if None)
            
        Returns:
            Created alert or None if rate-limited
        """
        with self.lock:
            # Create alert
            alert = Alert(
                title=title,
                message=message,
                severity=severity,
                source=source,
                context=context,
                fingerprint=fingerprint
            )
            
            # Check for rate limiting
            if self._should_rate_limit(alert):
                self.logger.info(
                    f"Alert rate-limited: {alert.title} (fingerprint: {alert.fingerprint})"
                )
                return None
            
            # Store the alert
            self.alerts[alert.alert_id] = alert
            self.active_alerts[alert.fingerprint] = alert
            
            # Update last alert time
            self.last_alert_times[alert.fingerprint] = alert.created_at
            
            # Get channels to use
            if channels is None:
                channels = self.default_channels.get(severity, [AlertChannel.CONSOLE])
            
            # Send to channels
            for channel_type in channels:
                if channel_type in self.channels:
                    channel = self.channels[channel_type]
                    
                    try:
                        # Send the alert
                        success, details = channel.send(alert)
                        
                        # Record in alert history
                        alert.add_notification(
                            channel=channel_type,
                            sent_at=time.time(),
                            success=success,
                            details=details
                        )
                        
                        if not success:
                            self.logger.warning(
                                f"Failed to send alert to {channel_type.value}: {details}"
                            )
                    except Exception as e:
                        self.logger.error(
                            f"Error sending alert to {channel_type.value}: {str(e)}"
                        )
                        
                        # Record error in alert history
                        alert.add_notification(
                            channel=channel_type,
                            sent_at=time.time(),
                            success=False,
                            details=str(e)
                        )
            
            # Notify alert handlers
            for handler in self.alert_handlers:
                try:
                    handler(alert)
                except Exception as e:
                    self.logger.error(f"Error in alert handler: {str(e)}")
            
            return alert
    
    def _should_rate_limit(self, alert: Alert) -> bool:
        """
        Check if an alert should be rate-limited.
        
        Args:
            alert: Alert to check
            
        Returns:
            True if the alert should be rate-limited
        """
        # Check if we've sent this alert recently
        if alert.fingerprint in self.last_alert_times:
            last_time = self.last_alert_times[alert.fingerprint]
            cooldown = self.cooldown_periods.get(alert.severity, 3600)
            
            if time.time() - last_time < cooldown:
                return True
        
        return False
    
    def acknowledge_alert(self, alert_id: str, user: Optional[str] = None) -> bool:
        """
        Acknowledge an alert.
        
        Args:
            alert_id: Alert ID
            user: User who acknowledged the alert
            
        Returns:
            True if the alert was found and acknowledged
        """
        with self.lock:
            if alert_id in self.alerts:
                alert = self.alerts[alert_id]
                alert.acknowledge(user)
                return True
            return False
    
    def resolve_alert(
        self,
        alert_id: str,
        user: Optional[str] = None,
        message: Optional[str] = None
    ) -> bool:
        """
        Resolve an alert.
        
        Args:
            alert_id: Alert ID
            user: User who resolved the alert
            message: Resolution message
            
        Returns:
            True if the alert was found and resolved
        """
        with self.lock:
            if alert_id in self.alerts:
                alert = self.alerts[alert_id]
                alert.resolve(user, message)
                
                # Remove from active alerts if it's the most recent for this fingerprint
                if (
                    alert.fingerprint in self.active_alerts and
                    self.active_alerts[alert.fingerprint].alert_id == alert_id
                ):
                    del self.active_alerts[alert.fingerprint]
                
                return True
            return False
    
    def resolve_by_fingerprint(
        self,
        fingerprint: str,
        user: Optional[str] = None,
        message: Optional[str] = None
    ) -> List[str]:
        """
        Resolve all alerts with a given fingerprint.
        
        Args:
            fingerprint: Alert fingerprint
            user: User who resolved the alerts
            message: Resolution message
            
        Returns:
            List of resolved alert IDs
        """
        with self.lock:
            resolved_ids = []
            
            # Find all alerts with the fingerprint
            for alert_id, alert in self.alerts.items():
                if alert.fingerprint == fingerprint and alert.status != AlertStatus.RESOLVED:
                    alert.resolve(user, message)
                    resolved_ids.append(alert_id)
            
            # Remove from active alerts
            if fingerprint in self.active_alerts:
                del self.active_alerts[fingerprint]
            
            return resolved_ids
    
    def mute_alert(self, alert_id: str) -> bool:
        """
        Mute an alert.
        
        Args:
            alert_id: Alert ID
            
        Returns:
            True if the alert was found and muted
        """
        with self.lock:
            if alert_id in self.alerts:
                alert = self.alerts[alert_id]
                alert.mute()
                return True
            return False
    
    def get_alert(self, alert_id: str) -> Optional[Alert]:
        """
        Get an alert by ID.
        
        Args:
            alert_id: Alert ID
            
        Returns:
            Alert or None if not found
        """
        with self.lock:
            return self.alerts.get(alert_id)
    
    def get_alerts(
        self,
        status: Optional[AlertStatus] = None,
        severity: Optional[AlertSeverity] = None,
        source: Optional[str] = None,
        start_time: Optional[float] = None,
        end_time: Optional[float] = None,
        limit: Optional[int] = None
    ) -> List[Alert]:
        """
        Get alerts filtered by various criteria.
        
        Args:
            status: Filter by alert status
            severity: Filter by alert severity
            source: Filter by alert source
            start_time: Filter by created_at >= start_time
            end_time: Filter by created_at <= end_time
            limit: Maximum number of alerts to return
            
        Returns:
            List of alerts matching the criteria
        """
        with self.lock:
            # Filter alerts
            filtered_alerts = []
            
            for alert in self.alerts.values():
                # Apply filters
                if status is not None and alert.status != status:
                    continue
                
                if severity is not None and alert.severity != severity:
                    continue
                
                if source is not None and alert.source != source:
                    continue
                
                if start_time is not None and alert.created_at < start_time:
                    continue
                
                if end_time is not None and alert.created_at > end_time:
                    continue
                
                filtered_alerts.append(alert)
            
            # Sort by created_at (newest first)
            filtered_alerts.sort(key=lambda a: a.created_at, reverse=True)
            
            # Apply limit
            if limit is not None:
                filtered_alerts = filtered_alerts[:limit]
            
            return filtered_alerts
    
    def get_active_alerts(
        self,
        severity: Optional[AlertSeverity] = None,
        source: Optional[str] = None
    ) -> List[Alert]:
        """
        Get active alerts (not resolved or muted).
        
        Args:
            severity: Filter by alert severity
            source: Filter by alert source
            
        Returns:
            List of active alerts
        """
        with self.lock:
            # Filter by status
            return self.get_alerts(
                status=None,  # Any status
                severity=severity,
                source=source
            )
    
    def register_alert_handler(self, handler: Callable[[Alert], None]) -> None:
        """
        Register a handler function for alerts.
        
        Args:
            handler: Function to call when an alert is created
        """
        with self.lock:
            self.alert_handlers.append(handler)
    
    def clear_alerts(self, older_than: Optional[float] = None) -> int:
        """
        Clear alerts from memory.
        
        Args:
            older_than: Clear alerts older than this timestamp
            
        Returns:
            Number of alerts cleared
        """
        with self.lock:
            # Determine which alerts to clear
            if older_than is None:
                # Clear all alerts
                count = len(self.alerts)
                self.alerts.clear()
                self.active_alerts.clear()
                self.last_alert_times.clear()
                return count
            
            # Clear alerts older than the specified time
            to_remove = []
            
            for alert_id, alert in self.alerts.items():
                if alert.created_at < older_than:
                    to_remove.append(alert_id)
            
            # Remove alerts
            for alert_id in to_remove:
                alert = self.alerts[alert_id]
                del self.alerts[alert_id]
                
                # Remove from active alerts if it's the most recent for this fingerprint
                if (
                    alert.fingerprint in self.active_alerts and
                    self.active_alerts[alert.fingerprint].alert_id == alert_id
                ):
                    del self.active_alerts[alert.fingerprint]
                
                # Remove from last alert times if it's the most recent for this fingerprint
                # (this will allow new alerts with this fingerprint)
                if (
                    alert.fingerprint in self.last_alert_times and
                    abs(self.last_alert_times[alert.fingerprint] - alert.created_at) < 0.1
                ):
                    del self.last_alert_times[alert.fingerprint]
            
            return len(to_remove)
    
    def get_stats(self) -> Dict[str, Any]:
        """
        Get alerting statistics.
        
        Returns:
            Dictionary with alerting statistics
        """
        with self.lock:
            # Count alerts by status
            status_counts = {status.value: 0 for status in AlertStatus}
            for alert in self.alerts.values():
                status_counts[alert.status.value] += 1
            
            # Count alerts by severity
            severity_counts = {severity.name: 0 for severity in AlertSeverity}
            for alert in self.alerts.values():
                severity_counts[alert.severity.name] += 1
            
            # Count alerts by source
            source_counts = {}
            for alert in self.alerts.values():
                if alert.source not in source_counts:
                    source_counts[alert.source] = 0
                source_counts[alert.source] += 1
            
            # Calculate alert age statistics
            now = time.time()
            ages = [now - alert.created_at for alert in self.alerts.values()]
            
            if ages:
                avg_age = sum(ages) / len(ages)
                max_age = max(ages)
                min_age = min(ages)
            else:
                avg_age = 0
                max_age = 0
                min_age = 0
            
            return {
                'total_alerts': len(self.alerts),
                'active_alerts': len(self.active_alerts),
                'status_counts': status_counts,
                'severity_counts': severity_counts,
                'source_counts': source_counts,
                'avg_age': avg_age,
                'max_age': max_age,
                'min_age': min_age
            }


# Global alert manager instance
alert_manager = AlertManager.get_instance()


def send_alert(
    title: str,
    message: str,
    severity: AlertSeverity,
    source: str,
    context: Optional[Dict[str, Any]] = None,
    fingerprint: Optional[str] = None,
    channels: Optional[List[AlertChannel]] = None
) -> Optional[Alert]:
    """
    Send an alert using the global alert manager.
    
    Args:
        title: Alert title
        message: Alert message
        severity: Alert severity
        source: Alert source
        context: Additional context
        fingerprint: Alert fingerprint for deduplication
        channels: Channels to send the alert to
        
    Returns:
        Created alert or None if rate-limited
    """
    return alert_manager.send_alert(
        title=title,
        message=message,
        severity=severity,
        source=source,
        context=context,
        fingerprint=fingerprint,
        channels=channels
    )


def register_alert_handler(handler: Callable[[Alert], None]) -> None:
    """
    Register a handler function for alerts with the global alert manager.
    
    Args:
        handler: Function to call when an alert is created
    """
    alert_manager.register_alert_handler(handler)


def register_email_channel(
    smtp_server: str,
    smtp_port: int,
    smtp_username: str,
    smtp_password: str,
    sender: str,
    recipients: List[str],
    use_tls: bool = True
) -> None:
    """
    Register an email notification channel with the global alert manager.
    
    Args:
        smtp_server: SMTP server hostname
        smtp_port: SMTP server port
        smtp_username: SMTP username
        smtp_password: SMTP password
        sender: Sender email address
        recipients: List of recipient email addresses
        use_tls: Whether to use TLS
    """
    channel = EmailChannel(
        smtp_server=smtp_server,
        smtp_port=smtp_port,
        smtp_username=smtp_username,
        smtp_password=smtp_password,
        sender=sender,
        recipients=recipients,
        use_tls=use_tls
    )
    
    alert_manager.register_channel(AlertChannel.EMAIL, channel)


def register_slack_channel(webhook_url: str) -> None:
    """
    Register a Slack notification channel with the global alert manager.
    
    Args:
        webhook_url: Slack webhook URL
    """
    channel = SlackChannel(webhook_url=webhook_url)
    alert_manager.register_channel(AlertChannel.SLACK, channel)


def register_webhook_channel(webhook_url: str, headers: Optional[Dict[str, str]] = None) -> None:
    """
    Register a webhook notification channel with the global alert manager.
    
    Args:
        webhook_url: Webhook URL
        headers: HTTP headers
    """
    channel = WebhookChannel(webhook_url=webhook_url, headers=headers)
    alert_manager.register_channel(AlertChannel.WEBHOOK, channel)


def acknowledge_alert(alert_id: str, user: Optional[str] = None) -> bool:
    """
    Acknowledge an alert using the global alert manager.
    
    Args:
        alert_id: Alert ID
        user: User who acknowledged the alert
        
    Returns:
        True if the alert was found and acknowledged
    """
    return alert_manager.acknowledge_alert(alert_id, user)


def resolve_alert(
    alert_id: str,
    user: Optional[str] = None,
    message: Optional[str] = None
) -> bool:
    """
    Resolve an alert using the global alert manager.
    
    Args:
        alert_id: Alert ID
        user: User who resolved the alert
        message: Resolution message
        
    Returns:
        True if the alert was found and resolved
    """
    return alert_manager.resolve_alert(alert_id, user, message)