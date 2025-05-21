"""
Monitoring and observability for the CryoProtect Unified Importer.

This module integrates progress tracking, observability, and alerting
capabilities for the unified importer.
"""

import os
import time
import json
import logging
import datetime
import threading
import smtplib
import tempfile
import webbrowser
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from typing import Dict, Any, Optional, List, Union, Callable, Set, Tuple
from dataclasses import dataclass, field

from .progress_tracking import ProgressTracker, ProgressStats, ReportFormat

class Alert:
    """Base class for alerting capabilities."""
    
    def __init__(self, name: str, description: str, severity: str = "info"):
        """Initialize an alert.
        
        Args:
            name: Name of the alert
            description: Description of the alert
            severity: Alert severity (info, warning, error, critical)
        """
        self.name = name
        self.description = description
        self.severity = severity
        self.timestamp = time.time()
        self.resolved = False
        self.resolved_timestamp = None
    
    def resolve(self) -> None:
        """Mark the alert as resolved."""
        self.resolved = True
        self.resolved_timestamp = time.time()
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert the alert to a dictionary.
        
        Returns:
            Dictionary representation of the alert
        """
        return {
            "name": self.name,
            "description": self.description,
            "severity": self.severity,
            "timestamp": self.timestamp,
            "timestamp_formatted": datetime.datetime.fromtimestamp(
                self.timestamp
            ).strftime("%Y-%m-%d %H:%M:%S"),
            "resolved": self.resolved,
            "resolved_timestamp": self.resolved_timestamp,
            "resolved_timestamp_formatted": (
                datetime.datetime.fromtimestamp(
                    self.resolved_timestamp
                ).strftime("%Y-%m-%d %H:%M:%S") if self.resolved_timestamp else None
            )
        }
    
    def __str__(self) -> str:
        """String representation of the alert.
        
        Returns:
            Formatted alert string
        """
        resolved_str = " (Resolved)" if self.resolved else ""
        timestamp = datetime.datetime.fromtimestamp(self.timestamp).strftime(
            "%Y-%m-%d %H:%M:%S"
        )
        return f"[{self.severity.upper()}] {self.name}{resolved_str} - {timestamp}: {self.description}"

class ThresholdAlert(Alert):
    """Alert based on a threshold condition."""
    
    def __init__(
        self,
        name: str,
        description: str,
        metric: str,
        threshold: float,
        comparison: str,
        current_value: float,
        severity: str = "warning"
    ):
        """Initialize a threshold alert.
        
        Args:
            name: Name of the alert
            description: Description of the alert
            metric: Name of the metric being monitored
            threshold: Threshold value
            comparison: Comparison operator (>, <, >=, <=, ==, !=)
            current_value: Current value of the metric
            severity: Alert severity (info, warning, error, critical)
        """
        super().__init__(name, description, severity)
        self.metric = metric
        self.threshold = threshold
        self.comparison = comparison
        self.current_value = current_value
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert the threshold alert to a dictionary.
        
        Returns:
            Dictionary representation of the alert
        """
        data = super().to_dict()
        data.update({
            "metric": self.metric,
            "threshold": self.threshold,
            "comparison": self.comparison,
            "current_value": self.current_value
        })
        return data
    
    def __str__(self) -> str:
        """String representation of the threshold alert.
        
        Returns:
            Formatted alert string
        """
        base_str = super().__str__()
        return f"{base_str} ({self.metric} {self.comparison} {self.threshold}, current: {self.current_value})"

class ErrorAlert(Alert):
    """Alert for error conditions."""
    
    def __init__(
        self,
        name: str,
        description: str,
        error: Exception,
        context: Dict[str, Any],
        severity: str = "error"
    ):
        """Initialize an error alert.
        
        Args:
            name: Name of the alert
            description: Description of the alert
            error: Exception that triggered the alert
            context: Additional context about the error
            severity: Alert severity (info, warning, error, critical)
        """
        super().__init__(name, description, severity)
        self.error = error
        self.error_type = type(error).__name__
        self.error_message = str(error)
        self.context = context
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert the error alert to a dictionary.
        
        Returns:
            Dictionary representation of the alert
        """
        data = super().to_dict()
        data.update({
            "error_type": self.error_type,
            "error_message": self.error_message,
            "context": self.context
        })
        return data
    
    def __str__(self) -> str:
        """String representation of the error alert.
        
        Returns:
            Formatted alert string
        """
        base_str = super().__str__()
        return f"{base_str} - {self.error_type}: {self.error_message}"

class AlertManager:
    """Manages alerts and notifications."""
    
    def __init__(
        self,
        logger: Optional[logging.Logger] = None,
        alert_log_file: Optional[str] = None,
        email_config: Optional[Dict[str, Any]] = None,
        webhook_urls: Optional[List[str]] = None,
        severity_thresholds: Optional[Dict[str, str]] = None
    ):
        """Initialize the alert manager.
        
        Args:
            logger: Optional logger instance
            alert_log_file: Optional file to log alerts
            email_config: Optional email configuration for notifications
            webhook_urls: Optional webhook URLs for notifications
            severity_thresholds: Optional severity thresholds for alerts
        """
        self.logger = logger or logging.getLogger(__name__)
        self.alert_log_file = alert_log_file
        self.email_config = email_config or {}
        self.webhook_urls = webhook_urls or []
        self.severity_thresholds = severity_thresholds or {
            "email": "error",
            "webhook": "warning",
            "log": "info"
        }
        
        self.alerts: List[Alert] = []
        self.active_thresholds: Dict[str, Dict[str, Any]] = {}
    
    def add_alert(self, alert: Alert) -> None:
        """Add an alert to the manager and trigger notifications.
        
        Args:
            alert: Alert to add
        """
        self.alerts.append(alert)
        
        # Log the alert
        if alert.severity == "critical":
            self.logger.critical(str(alert))
        elif alert.severity == "error":
            self.logger.error(str(alert))
        elif alert.severity == "warning":
            self.logger.warning(str(alert))
        else:
            self.logger.info(str(alert))
        
        # Write to alert log file if configured
        if self.alert_log_file:
            try:
                with open(self.alert_log_file, 'a') as f:
                    f.write(f"{str(alert)}\n")
            except Exception as e:
                self.logger.error(f"Failed to write to alert log file: {e}")
        
        # Send notifications based on severity thresholds
        self._send_notifications(alert)
    
    def _send_notifications(self, alert: Alert) -> None:
        """Send notifications based on alert severity.
        
        Args:
            alert: Alert to send notifications for
        """
        # Map severity string to numeric level for comparison
        severity_levels = {
            "info": 0,
            "warning": 1,
            "error": 2,
            "critical": 3
        }
        
        alert_level = severity_levels.get(alert.severity, 0)
        
        # Check email threshold
        email_threshold = severity_levels.get(
            self.severity_thresholds.get("email", "error"),
            2  # Default to error
        )
        
        if alert_level >= email_threshold and self.email_config:
            self._send_email_notification(alert)
        
        # Check webhook threshold
        webhook_threshold = severity_levels.get(
            self.severity_thresholds.get("webhook", "warning"),
            1  # Default to warning
        )
        
        if alert_level >= webhook_threshold and self.webhook_urls:
            self._send_webhook_notification(alert)
    
    def _send_email_notification(self, alert: Alert) -> None:
        """Send email notification for an alert.
        
        Args:
            alert: Alert to send notification for
        """
        if not self.email_config:
            return
        
        try:
            # Extract email config
            smtp_server = self.email_config.get("smtp_server", "localhost")
            smtp_port = self.email_config.get("smtp_port", 25)
            smtp_user = self.email_config.get("smtp_user")
            smtp_password = self.email_config.get("smtp_password")
            sender = self.email_config.get("sender", "unified_importer@example.com")
            recipients = self.email_config.get("recipients", [])
            
            if not recipients:
                self.logger.warning("No email recipients configured, skipping notification")
                return
            
            # Create message
            msg = MIMEMultipart()
            msg["From"] = sender
            msg["To"] = ", ".join(recipients)
            msg["Subject"] = f"[{alert.severity.upper()}] CryoProtect Alert: {alert.name}"
            
            # Create message body
            body = f"""
            <html>
            <body>
                <h2>CryoProtect Unified Importer Alert</h2>
                <p><strong>Name:</strong> {alert.name}</p>
                <p><strong>Severity:</strong> {alert.severity.upper()}</p>
                <p><strong>Time:</strong> {datetime.datetime.fromtimestamp(alert.timestamp).strftime('%Y-%m-%d %H:%M:%S')}</p>
                <p><strong>Description:</strong> {alert.description}</p>
            """
            
            # Add additional details based on alert type
            if isinstance(alert, ThresholdAlert):
                body += f"""
                <p><strong>Metric:</strong> {alert.metric}</p>
                <p><strong>Threshold:</strong> {alert.comparison} {alert.threshold}</p>
                <p><strong>Current Value:</strong> {alert.current_value}</p>
                """
            elif isinstance(alert, ErrorAlert):
                body += f"""
                <p><strong>Error Type:</strong> {alert.error_type}</p>
                <p><strong>Error Message:</strong> {alert.error_message}</p>
                <p><strong>Context:</strong> {json.dumps(alert.context, indent=2)}</p>
                """
            
            body += """
            </body>
            </html>
            """
            
            msg.attach(MIMEText(body, "html"))
            
            # Send email
            with smtplib.SMTP(smtp_server, smtp_port) as server:
                if smtp_user and smtp_password:
                    server.login(smtp_user, smtp_password)
                server.send_message(msg)
            
            self.logger.info(f"Sent email notification for alert: {alert.name}")
        except Exception as e:
            self.logger.error(f"Failed to send email notification: {e}")
    
    def _send_webhook_notification(self, alert: Alert) -> None:
        """Send webhook notification for an alert.
        
        Args:
            alert: Alert to send notification for
        """
        if not self.webhook_urls:
            return
        
        try:
            import requests
            
            # Create payload
            payload = {
                "alert": alert.to_dict()
            }
            
            # Send to all webhook URLs
            for url in self.webhook_urls:
                try:
                    response = requests.post(
                        url,
                        json=payload,
                        headers={"Content-Type": "application/json"}
                    )
                    
                    if response.status_code >= 400:
                        self.logger.warning(
                            f"Webhook notification failed with status {response.status_code}: {response.text}"
                        )
                    else:
                        self.logger.info(f"Sent webhook notification to {url}")
                except Exception as e:
                    self.logger.error(f"Failed to send webhook notification to {url}: {e}")
        except ImportError:
            self.logger.error("Failed to import requests module for webhook notifications")
    
    def add_threshold_alert(
        self,
        name: str,
        metric: str,
        threshold: float,
        comparison: str,
        current_value: float,
        description: str = None,
        severity: str = "warning"
    ) -> Alert:
        """Add a threshold alert.
        
        Args:
            name: Name of the alert
            metric: Name of the metric being monitored
            threshold: Threshold value
            comparison: Comparison operator (>, <, >=, <=, ==, !=)
            current_value: Current value of the metric
            description: Optional description of the alert
            severity: Alert severity (info, warning, error, critical)
            
        Returns:
            Created alert
        """
        # Generate description if not provided
        if description is None:
            if comparison == ">":
                description = f"{metric} exceeds threshold of {threshold} (current: {current_value})"
            elif comparison == "<":
                description = f"{metric} below threshold of {threshold} (current: {current_value})"
            elif comparison == ">=":
                description = f"{metric} equals or exceeds threshold of {threshold} (current: {current_value})"
            elif comparison == "<=":
                description = f"{metric} equals or is below threshold of {threshold} (current: {current_value})"
            elif comparison == "==":
                description = f"{metric} equals threshold of {threshold}"
            elif comparison == "!=":
                description = f"{metric} does not equal expected value of {threshold} (current: {current_value})"
            else:
                description = f"{metric} {comparison} {threshold} (current: {current_value})"
        
        # Create alert
        alert = ThresholdAlert(
            name=name,
            description=description,
            metric=metric,
            threshold=threshold,
            comparison=comparison,
            current_value=current_value,
            severity=severity
        )
        
        # Add to manager
        self.add_alert(alert)
        
        return alert
    
    def add_error_alert(
        self,
        name: str,
        error: Exception,
        context: Dict[str, Any] = None,
        description: str = None,
        severity: str = "error"
    ) -> Alert:
        """Add an error alert.
        
        Args:
            name: Name of the alert
            error: Exception that triggered the alert
            context: Additional context about the error
            description: Optional description of the alert
            severity: Alert severity (info, warning, error, critical)
            
        Returns:
            Created alert
        """
        # Generate description if not provided
        if description is None:
            description = f"Error: {type(error).__name__}: {str(error)}"
        
        # Create alert
        alert = ErrorAlert(
            name=name,
            description=description,
            error=error,
            context=context or {},
            severity=severity
        )
        
        # Add to manager
        self.add_alert(alert)
        
        return alert
    
    def check_threshold(
        self,
        name: str,
        metric: str,
        threshold: float,
        comparison: str,
        current_value: float,
        description: str = None,
        severity: str = "warning",
        auto_resolve: bool = True
    ) -> Optional[Alert]:
        """Check a threshold and create an alert if exceeded.
        
        Args:
            name: Name of the alert
            metric: Name of the metric being monitored
            threshold: Threshold value
            comparison: Comparison operator (>, <, >=, <=, ==, !=)
            current_value: Current value of the metric
            description: Optional description of the alert
            severity: Alert severity (info, warning, error, critical)
            auto_resolve: Whether to automatically resolve alerts when threshold not exceeded
            
        Returns:
            Created alert if threshold is exceeded, otherwise None
        """
        # Check if threshold is exceeded
        threshold_exceeded = False
        
        if comparison == ">":
            threshold_exceeded = current_value > threshold
        elif comparison == "<":
            threshold_exceeded = current_value < threshold
        elif comparison == ">=":
            threshold_exceeded = current_value >= threshold
        elif comparison == "<=":
            threshold_exceeded = current_value <= threshold
        elif comparison == "==":
            threshold_exceeded = current_value == threshold
        elif comparison == "!=":
            threshold_exceeded = current_value != threshold
        else:
            self.logger.warning(f"Unknown comparison operator: {comparison}")
            return None
        
        # Check if alert exists for this threshold
        threshold_key = f"{name}:{metric}:{comparison}:{threshold}"
        existing_alert = self.active_thresholds.get(threshold_key)
        
        if threshold_exceeded:
            if existing_alert and not existing_alert["alert"].resolved:
                # Threshold still exceeded, update the alert
                existing_alert["last_checked"] = time.time()
                existing_alert["current_value"] = current_value
                return existing_alert["alert"]
            else:
                # New threshold exceeded, create alert
                alert = self.add_threshold_alert(
                    name=name,
                    metric=metric,
                    threshold=threshold,
                    comparison=comparison,
                    current_value=current_value,
                    description=description,
                    severity=severity
                )
                
                # Add to active thresholds
                self.active_thresholds[threshold_key] = {
                    "alert": alert,
                    "last_checked": time.time(),
                    "current_value": current_value
                }
                
                return alert
        elif existing_alert and not existing_alert["alert"].resolved and auto_resolve:
            # Threshold no longer exceeded, resolve the alert
            existing_alert["alert"].resolve()
            
            # Log resolution
            self.logger.info(
                f"Resolved alert {existing_alert['alert'].name}: "
                f"{metric} no longer {comparison} {threshold} (current: {current_value})"
            )
            
            # Remove from active thresholds
            self.active_thresholds.pop(threshold_key, None)
        
        return None
    
    def get_active_alerts(self) -> List[Alert]:
        """Get all active (unresolved) alerts.
        
        Returns:
            List of active alerts
        """
        return [alert for alert in self.alerts if not alert.resolved]
    
    def get_alerts_by_severity(self, severity: str) -> List[Alert]:
        """Get alerts by severity.
        
        Args:
            severity: Severity level to filter by
            
        Returns:
            List of alerts with the specified severity
        """
        return [alert for alert in self.alerts if alert.severity == severity]
    
    def get_alert_summary(self) -> Dict[str, Any]:
        """Get a summary of alerts.
        
        Returns:
            Dictionary with alert summary statistics
        """
        # Get counts by severity
        severity_counts = {
            "info": 0,
            "warning": 0,
            "error": 0,
            "critical": 0
        }
        
        active_severity_counts = severity_counts.copy()
        
        for alert in self.alerts:
            severity_counts[alert.severity] = severity_counts.get(alert.severity, 0) + 1
            
            if not alert.resolved:
                active_severity_counts[alert.severity] = active_severity_counts.get(alert.severity, 0) + 1
        
        return {
            "total_alerts": len(self.alerts),
            "active_alerts": len(self.get_active_alerts()),
            "severity_counts": severity_counts,
            "active_severity_counts": active_severity_counts,
            "oldest_active_alert": min([a.timestamp for a in self.get_active_alerts()]) if self.get_active_alerts() else None,
            "newest_active_alert": max([a.timestamp for a in self.get_active_alerts()]) if self.get_active_alerts() else None
        }

class DashboardGenerator:
    """Generates HTML dashboards for monitoring."""
    
    def __init__(
        self,
        title: str = "CryoProtect Unified Importer Dashboard",
        refresh_interval: int = 30,
        logger: Optional[logging.Logger] = None
    ):
        """Initialize the dashboard generator.
        
        Args:
            title: Dashboard title
            refresh_interval: Page refresh interval in seconds
            logger: Optional logger instance
        """
        self.title = title
        self.refresh_interval = refresh_interval
        self.logger = logger or logging.getLogger(__name__)
    
    def generate_dashboard(
        self,
        progress_trackers: List[ProgressTracker],
        alert_manager: Optional[AlertManager] = None,
        output_file: Optional[str] = None,
        open_browser: bool = False
    ) -> str:
        """Generate an HTML dashboard.
        
        Args:
            progress_trackers: List of progress trackers to include
            alert_manager: Optional alert manager for alerts section
            output_file: Optional file to save the dashboard
            open_browser: Whether to open the dashboard in a browser
            
        Returns:
            HTML dashboard content
        """
        # Generate HTML content
        html = self._generate_html(progress_trackers, alert_manager)
        
        # Save to file if specified
        if output_file:
            try:
                with open(output_file, 'w') as f:
                    f.write(html)
                self.logger.info(f"Saved dashboard to {output_file}")
            except Exception as e:
                self.logger.error(f"Failed to save dashboard: {e}")
        
        # Open in browser if requested
        if open_browser:
            self._open_in_browser(html, output_file)
        
        return html
    
    def _generate_html(
        self,
        progress_trackers: List[ProgressTracker],
        alert_manager: Optional[AlertManager] = None
    ) -> str:
        """Generate HTML content for the dashboard.
        
        Args:
            progress_trackers: List of progress trackers to include
            alert_manager: Optional alert manager for alerts section
            
        Returns:
            HTML content
        """
        # Start HTML content
        html = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <meta charset="UTF-8">
            <meta http-equiv="refresh" content="{self.refresh_interval}">
            <title>{self.title}</title>
            <style>
                body {{
                    font-family: Arial, sans-serif;
                    margin: 20px;
                    line-height: 1.6;
                }}
                h1, h2, h3 {{
                    color: #333;
                }}
                .dashboard-container {{
                    max-width: 1200px;
                    margin: 0 auto;
                }}
                .progress-container {{
                    margin-bottom: 30px;
                    border: 1px solid #ddd;
                    border-radius: 5px;
                    padding: 15px;
                }}
                .progress-header {{
                    display: flex;
                    justify-content: space-between;
                    align-items: center;
                    margin-bottom: 10px;
                }}
                .progress-bar-container {{
                    width: 100%;
                    height: 20px;
                    background-color: #f3f3f3;
                    border-radius: 4px;
                    overflow: hidden;
                    margin-bottom: 10px;
                }}
                .progress-bar {{
                    height: 100%;
                    background-color: #4CAF50;
                    text-align: center;
                    color: white;
                    font-weight: bold;
                }}
                .progress-stats {{
                    display: flex;
                    flex-wrap: wrap;
                    gap: 20px;
                    margin-bottom: 15px;
                }}
                .stat-item {{
                    flex: 1;
                    min-width: 150px;
                    background-color: #f9f9f9;
                    padding: 10px;
                    border-radius: 4px;
                }}
                .stat-value {{
                    font-size: 1.2em;
                    font-weight: bold;
                }}
                .errors-warning {{
                    background-color: #FFF9C4;
                }}
                .checkpoint-container {{
                    margin-top: 15px;
                }}
                .alert-container {{
                    margin-bottom: 30px;
                    border: 1px solid #ddd;
                    border-radius: 5px;
                    padding: 15px;
                }}
                .alert {{
                    padding: 10px;
                    border-radius: 4px;
                    margin-bottom: 10px;
                }}
                .alert-info {{
                    background-color: #E3F2FD;
                    border-left: 4px solid #2196F3;
                }}
                .alert-warning {{
                    background-color: #FFF9C4;
                    border-left: 4px solid #FFC107;
                }}
                .alert-error {{
                    background-color: #FFEBEE;
                    border-left: 4px solid #F44336;
                }}
                .alert-critical {{
                    background-color: #7E0000;
                    color: white;
                    border-left: 4px solid #B71C1C;
                }}
                .timestamp {{
                    font-size: 0.8em;
                    color: #666;
                }}
                .footer {{
                    text-align: center;
                    margin-top: 30px;
                    font-size: 0.8em;
                    color: #666;
                }}
                details {{
                    margin-top: 10px;
                }}
                summary {{
                    cursor: pointer;
                    font-weight: bold;
                }}
                table {{
                    width: 100%;
                    border-collapse: collapse;
                }}
                th, td {{
                    text-align: left;
                    padding: 8px;
                    border-bottom: 1px solid #ddd;
                }}
                th {{
                    background-color: #f3f3f3;
                }}
                .resolved {{
                    text-decoration: line-through;
                    opacity: 0.7;
                }}
                .dashboard-header {{
                    display: flex;
                    justify-content: space-between;
                    align-items: center;
                    margin-bottom: 20px;
                }}
                .dashboard-summary {{
                    display: flex;
                    gap: 15px;
                }}
                .summary-item {{
                    padding: 8px 15px;
                    border-radius: 4px;
                    font-weight: bold;
                }}
            </style>
        </head>
        <body>
            <div class="dashboard-container">
                <div class="dashboard-header">
                    <h1>{self.title}</h1>
                    <div class="dashboard-summary">
        """
        
        # Add alert summary if available
        if alert_manager:
            alert_summary = alert_manager.get_alert_summary()
            active_alerts = alert_summary["active_alerts"]
            
            if active_alerts > 0:
                alert_class = "alert-critical" if alert_summary["active_severity_counts"].get("critical", 0) > 0 else \
                             "alert-error" if alert_summary["active_severity_counts"].get("error", 0) > 0 else \
                             "alert-warning" if alert_summary["active_severity_counts"].get("warning", 0) > 0 else \
                             "alert-info"
                html += f"""
                        <div class="summary-item {alert_class}">
                            Active Alerts: {active_alerts}
                        </div>
                """
        
        # Add import progress summary
        total_items = sum(tracker.stats.total_items for tracker in progress_trackers)
        processed_items = sum(tracker.stats.processed_items for tracker in progress_trackers)
        percent_complete = (processed_items / total_items * 100) if total_items > 0 else 0
        
        html += f"""
                        <div class="summary-item">
                            Overall Progress: {percent_complete:.1f}%
                        </div>
                    </div>
                </div>
        """
        
        # Add alerts section if available
        if alert_manager:
            active_alerts = alert_manager.get_active_alerts()
            
            html += f"""
                <div class="alert-container">
                    <h2>Alerts ({len(active_alerts)} active)</h2>
            """
            
            if active_alerts:
                for alert in active_alerts:
                    alert_class = f"alert-{alert.severity}"
                    timestamp = datetime.datetime.fromtimestamp(alert.timestamp).strftime("%Y-%m-%d %H:%M:%S")
                    
                    html += f"""
                        <div class="alert {alert_class}">
                            <div><strong>{alert.name}</strong></div>
                            <div>{alert.description}</div>
                            <div class="timestamp">{timestamp}</div>
                    """
                    
                    # Add additional details based on alert type
                    if isinstance(alert, ThresholdAlert):
                        html += f"""
                            <div>
                                Metric: {alert.metric}, Threshold: {alert.comparison} {alert.threshold}, 
                                Current: {alert.current_value}
                            </div>
                        """
                    elif isinstance(alert, ErrorAlert):
                        html += f"""
                            <div>
                                Error: {alert.error_type}: {alert.error_message}
                            </div>
                        """
                    
                    html += """
                        </div>
                    """
            else:
                html += "<p>No active alerts</p>"
            
            html += """
                </div>
            """
        
        # Add progress trackers
        for i, tracker in enumerate(progress_trackers):
            stats = tracker.stats
            
            # Calculate percent complete
            percent_complete = (stats.processed_items / stats.total_items * 100) if stats.total_items > 0 else 0
            
            html += f"""
                <div class="progress-container">
                    <div class="progress-header">
                        <h2>{tracker.description or f'Import {i+1}'}</h2>
                        <div>
                            {stats.processed_items}/{stats.total_items} ({percent_complete:.1f}%)
                        </div>
                    </div>
                    
                    <div class="progress-bar-container">
                        <div class="progress-bar" style="width: {percent_complete}%">
                            {percent_complete:.1f}%
                        </div>
                    </div>
                    
                    <div class="progress-stats">
                        <div class="stat-item">
                            <div>Elapsed Time</div>
                            <div class="stat-value">{self._format_duration(stats.elapsed_time)}</div>
                        </div>
            """
            
            # Add ETA if available
            if stats.estimated_completion_time:
                eta_str = self._format_duration(stats.estimated_time_remaining)
                completion_time = datetime.datetime.fromtimestamp(stats.estimated_completion_time).strftime("%H:%M:%S")
                
                html += f"""
                        <div class="stat-item">
                            <div>Estimated Time Remaining</div>
                            <div class="stat-value">{eta_str}</div>
                            <div>Completion at {completion_time}</div>
                        </div>
                """
            
            # Add success/error counts
            errors_class = " errors-warning" if stats.error_count > 0 else ""
            html += f"""
                        <div class="stat-item{errors_class}">
                            <div>Results</div>
                            <div class="stat-value">
                                {stats.success_count} successful, {stats.error_count} errors
                            </div>
                        </div>
            """
            
            # Add performance metrics
            html += f"""
                        <div class="stat-item">
                            <div>Performance</div>
                            <div class="stat-value">
                                {stats.items_per_second:.2f} items/sec
                            </div>
                        </div>
                    </div>
            """
            
            # Add checkpoints if available
            if stats.checkpoint_times:
                html += """
                    <details class="checkpoint-container">
                        <summary>Checkpoints</summary>
                        <table>
                            <tr>
                                <th>Checkpoint</th>
                                <th>Time</th>
                                <th>Elapsed</th>
                                <th>Items Processed</th>
                            </tr>
                """
                
                for cp in stats.checkpoint_times:
                    cp_time = datetime.datetime.fromtimestamp(cp["time"]).strftime("%H:%M:%S")
                    cp_elapsed = self._format_duration(cp["elapsed"])
                    
                    html += f"""
                            <tr>
                                <td>{cp["name"]}</td>
                                <td>{cp_time}</td>
                                <td>{cp_elapsed}</td>
                                <td>{cp["processed_items"]}/{stats.total_items}</td>
                            </tr>
                    """
                
                html += """
                        </table>
                    </details>
                """
            
            html += """
                </div>
            """
        
        # Add footer and close HTML
        update_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        html += f"""
                <div class="footer">
                    <p>Last updated: {update_time} · Auto-refresh every {self.refresh_interval} seconds</p>
                </div>
            </div>
        </body>
        </html>
        """
        
        return html
    
    def _format_duration(self, seconds: float) -> str:
        """Format duration in a human-readable way.
        
        Args:
            seconds: Duration in seconds
            
        Returns:
            Formatted duration string
        """
        if seconds < 1:
            return f"{seconds*1000:.0f}ms"
        elif seconds < 60:
            return f"{seconds:.1f}s"
        elif seconds < 3600:
            minutes = int(seconds / 60)
            secs = seconds % 60
            return f"{minutes}m {secs:.0f}s"
        elif seconds < 86400:
            hours = int(seconds / 3600)
            minutes = int((seconds % 3600) / 60)
            return f"{hours}h {minutes}m"
        else:
            days = int(seconds / 86400)
            hours = int((seconds % 86400) / 3600)
            return f"{days}d {hours}h"
    
    def _open_in_browser(self, html: str, output_file: Optional[str] = None) -> None:
        """Open the dashboard in a web browser.
        
        Args:
            html: HTML content to display
            output_file: Optional file path to open directly
        """
        try:
            if output_file:
                # Open the saved file
                webbrowser.open(f"file://{os.path.abspath(output_file)}")
            else:
                # Create a temporary file
                with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as temp_file:
                    temp_file_path = temp_file.name
                    temp_file.write(html.encode('utf-8'))
                
                # Open the temporary file
                webbrowser.open(f"file://{temp_file_path}")
                
                # Schedule deletion after a delay
                def delete_temp_file():
                    time.sleep(60)  # Wait for 60 seconds
                    try:
                        os.remove(temp_file_path)
                    except Exception:
                        pass
                
                threading.Thread(target=delete_temp_file, daemon=True).start()
        except Exception as e:
            self.logger.error(f"Failed to open dashboard in browser: {e}")

class MonitoringManager:
    """Manages monitoring and observability for the unified importer."""
    
    def __init__(
        self,
        logger: Optional[logging.Logger] = None,
        checkpoint_dir: Optional[str] = None,
        dashboard_dir: Optional[str] = None,
        alert_log_file: Optional[str] = None,
        email_config: Optional[Dict[str, Any]] = None,
        webhook_urls: Optional[List[str]] = None
    ):
        """Initialize the monitoring manager.
        
        Args:
            logger: Optional logger instance
            checkpoint_dir: Optional directory for progress checkpoints
            dashboard_dir: Optional directory for dashboards
            alert_log_file: Optional file to log alerts
            email_config: Optional email configuration for notifications
            webhook_urls: Optional webhook URLs for notifications
        """
        self.logger = logger or logging.getLogger(__name__)
        self.checkpoint_dir = checkpoint_dir
        self.dashboard_dir = dashboard_dir
        
        # Create directories if they don't exist
        if checkpoint_dir and not os.path.exists(checkpoint_dir):
            os.makedirs(checkpoint_dir)
        
        if dashboard_dir and not os.path.exists(dashboard_dir):
            os.makedirs(dashboard_dir)
        
        # Initialize components
        self.alert_manager = AlertManager(
            logger=self.logger,
            alert_log_file=alert_log_file,
            email_config=email_config,
            webhook_urls=webhook_urls
        )
        
        self.dashboard_generator = DashboardGenerator(
            title="CryoProtect Unified Importer Dashboard",
            logger=self.logger
        )
        
        self.progress_trackers: Dict[str, ProgressTracker] = {}
        self.monitoring_thread = None
        self.stop_monitoring = False
    
    def create_progress_tracker(
        self,
        name: str,
        total_items: int,
        description: str = None,
        update_interval: float = 1.0
    ) -> ProgressTracker:
        """Create a new progress tracker.
        
        Args:
            name: Unique name for the tracker
            total_items: Total number of items to process
            description: Optional description
            update_interval: Update interval in seconds
            
        Returns:
            Created progress tracker
        """
        checkpoint_file = None
        if self.checkpoint_dir:
            checkpoint_file = os.path.join(self.checkpoint_dir, f"{name}_checkpoint.json")
        
        tracker = ProgressTracker(
            total_items=total_items,
            description=description or name,
            logger=self.logger,
            update_interval=update_interval,
            checkpoint_file=checkpoint_file,
            report_format=ReportFormat.CONSOLE,
            on_update=lambda stats: self._on_progress_update(name, stats)
        )
        
        self.progress_trackers[name] = tracker
        return tracker
    
    def _on_progress_update(self, name: str, stats: ProgressStats) -> None:
        """Handle progress updates and check for threshold alerts.
        
        Args:
            name: Name of the tracker
            stats: Updated progress stats
        """
        # Check for slow progress
        if (
            stats.processed_items > 0 and 
            stats.elapsed_time > 60 and  # Only check after 1 minute
            stats.items_per_second < 0.1  # Less than 6 items per minute
        ):
            self.alert_manager.check_threshold(
                name=f"Slow Progress - {name}",
                metric="items_per_second",
                threshold=0.1,
                comparison="<",
                current_value=stats.items_per_second,
                description=f"Import {name} is processing items very slowly ({stats.items_per_second:.2f} items/sec)",
                severity="warning"
            )
        
        # Check for high error rate
        if stats.processed_items > 10:  # Only check after processing some items
            error_rate = stats.error_count / stats.processed_items
            if error_rate > 0.2:  # More than 20% errors
                self.alert_manager.check_threshold(
                    name=f"High Error Rate - {name}",
                    metric="error_rate",
                    threshold=0.2,
                    comparison=">",
                    current_value=error_rate,
                    description=f"Import {name} has a high error rate ({error_rate:.1%})",
                    severity="error" if error_rate > 0.5 else "warning"
                )
    
    def generate_dashboard(
        self,
        output_file: Optional[str] = None,
        open_browser: bool = False
    ) -> str:
        """Generate a monitoring dashboard.
        
        Args:
            output_file: Optional file to save the dashboard
            open_browser: Whether to open the dashboard in a browser
            
        Returns:
            HTML dashboard content
        """
        # Generate default output file if not specified
        if output_file is None and self.dashboard_dir:
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            output_file = os.path.join(self.dashboard_dir, f"dashboard_{timestamp}.html")
        
        # Generate dashboard
        return self.dashboard_generator.generate_dashboard(
            progress_trackers=list(self.progress_trackers.values()),
            alert_manager=self.alert_manager,
            output_file=output_file,
            open_browser=open_browser
        )
    
    def start_monitoring(
        self,
        dashboard_interval: int = 60,
        output_file: Optional[str] = None,
        open_browser: bool = False
    ) -> None:
        """Start continuous monitoring in a background thread.
        
        Args:
            dashboard_interval: Interval for dashboard updates in seconds
            output_file: Optional file to save the dashboard
            open_browser: Whether to open the dashboard in a browser
        """
        if self.monitoring_thread and self.monitoring_thread.is_alive():
            self.logger.warning("Monitoring thread is already running")
            return
        
        self.stop_monitoring = False
        
        # Define monitoring function
        def monitoring_thread_func():
            self.logger.info("Started monitoring thread")
            
            last_dashboard_time = 0
            
            while not self.stop_monitoring:
                current_time = time.time()
                
                # Generate dashboard at specified interval
                if current_time - last_dashboard_time >= dashboard_interval:
                    try:
                        self.generate_dashboard(
                            output_file=output_file,
                            open_browser=open_browser and last_dashboard_time == 0  # Only open browser first time
                        )
                        last_dashboard_time = current_time
                    except Exception as e:
                        self.logger.error(f"Error generating dashboard: {e}")
                
                # Sleep briefly to avoid high CPU usage
                time.sleep(1)
            
            self.logger.info("Stopped monitoring thread")
        
        # Start monitoring thread
        self.monitoring_thread = threading.Thread(
            target=monitoring_thread_func,
            daemon=True
        )
        self.monitoring_thread.start()
    
    def stop_monitoring_thread(self) -> None:
        """Stop the monitoring thread."""
        if self.monitoring_thread and self.monitoring_thread.is_alive():
            self.stop_monitoring = True
            self.monitoring_thread.join(timeout=5.0)
            
            if self.monitoring_thread.is_alive():
                self.logger.warning("Monitoring thread did not stop cleanly")
        
        self.monitoring_thread = None
    
    def add_alert(self, alert: Alert) -> None:
        """Add an alert.
        
        Args:
            alert: Alert to add
        """
        self.alert_manager.add_alert(alert)
    
    def add_threshold_alert(
        self,
        name: str,
        metric: str,
        threshold: float,
        comparison: str,
        current_value: float,
        description: str = None,
        severity: str = "warning"
    ) -> Alert:
        """Add a threshold alert.
        
        Args:
            name: Name of the alert
            metric: Name of the metric being monitored
            threshold: Threshold value
            comparison: Comparison operator (>, <, >=, <=, ==, !=)
            current_value: Current value of the metric
            description: Optional description of the alert
            severity: Alert severity (info, warning, error, critical)
            
        Returns:
            Created alert
        """
        return self.alert_manager.add_threshold_alert(
            name=name,
            metric=metric,
            threshold=threshold,
            comparison=comparison,
            current_value=current_value,
            description=description,
            severity=severity
        )
    
    def add_error_alert(
        self,
        name: str,
        error: Exception,
        context: Dict[str, Any] = None,
        description: str = None,
        severity: str = "error"
    ) -> Alert:
        """Add an error alert.
        
        Args:
            name: Name of the alert
            error: Exception that triggered the alert
            context: Additional context about the error
            description: Optional description of the alert
            severity: Alert severity (info, warning, error, critical)
            
        Returns:
            Created alert
        """
        return self.alert_manager.add_error_alert(
            name=name,
            error=error,
            context=context,
            description=description,
            severity=severity
        )
    
    def check_threshold(
        self,
        name: str,
        metric: str,
        threshold: float,
        comparison: str,
        current_value: float,
        description: str = None,
        severity: str = "warning",
        auto_resolve: bool = True
    ) -> Optional[Alert]:
        """Check a threshold and create an alert if exceeded.
        
        Args:
            name: Name of the alert
            metric: Name of the metric being monitored
            threshold: Threshold value
            comparison: Comparison operator (>, <, >=, <=, ==, !=)
            current_value: Current value of the metric
            description: Optional description of the alert
            severity: Alert severity (info, warning, error, critical)
            auto_resolve: Whether to automatically resolve alerts when threshold not exceeded
            
        Returns:
            Created alert if threshold is exceeded, otherwise None
        """
        return self.alert_manager.check_threshold(
            name=name,
            metric=metric,
            threshold=threshold,
            comparison=comparison,
            current_value=current_value,
            description=description,
            severity=severity,
            auto_resolve=auto_resolve
        )
    
    def get_active_alerts(self) -> List[Alert]:
        """Get all active (unresolved) alerts.
        
        Returns:
            List of active alerts
        """
        return self.alert_manager.get_active_alerts()
    
    def get_alert_summary(self) -> Dict[str, Any]:
        """Get a summary of alerts.
        
        Returns:
            Dictionary with alert summary statistics
        """
        return self.alert_manager.get_alert_summary()