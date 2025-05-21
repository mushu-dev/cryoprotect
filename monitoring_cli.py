#!/usr/bin/env python3
"""
Monitoring CLI Tool for CryoProtect

This command-line tool provides access to the unified monitoring system,
allowing users to:
1. Check system health status
2. View performance metrics
3. Track ongoing operations
4. Manage alerts and notifications
5. Generate reports

Usage:
    python monitoring_cli.py [command] [options]

Commands:
    status              Show current system health status
    metrics             Show performance metrics
    progress            Show progress of ongoing operations
    alerts              Show recent alerts
    errors              Show recent errors
    silence             Silence alerts by type
    unsilence           Unsilence alerts by type
    report              Generate monitoring report
    dashboard           Open the monitoring dashboard in a browser
"""

import os
import sys
import json
import time
import argparse
import socket
import webbrowser
from datetime import datetime
from typing import Dict, Any, List, Optional, Union

# Check if monitoring module exists
try:
    from unified_monitoring import MonitoringService
    monitoring_available = True
except ImportError:
    monitoring_available = False

# Terminal colors for better readability
class Colors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def colored(text, color):
    """Add color to text if supported by terminal."""
    if os.environ.get('TERM') == 'dumb' or not sys.stdout.isatty():
        return text
    return f"{color}{text}{Colors.ENDC}"

def get_status_color(status):
    """Get color for status."""
    status = status.lower()
    if status in ('ok', 'healthy'):
        return Colors.GREEN
    elif status in ('warning'):
        return Colors.YELLOW
    elif status in ('error', 'unhealthy', 'critical'):
        return Colors.RED
    else:
        return Colors.CYAN

def format_timestamp(timestamp):
    """Format timestamp as human-readable."""
    try:
        return datetime.fromtimestamp(timestamp).strftime('%Y-%m-%d %H:%M:%S')
    except (TypeError, ValueError):
        return str(timestamp)

def format_duration(seconds):
    """Format duration in seconds as human-readable."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        minutes = seconds / 60
        return f"{int(minutes)}m {int(seconds % 60)}s"
    else:
        hours = seconds / 3600
        minutes = (seconds % 3600) / 60
        return f"{int(hours)}h {int(minutes)}m"

def format_bytes(size):
    """Format bytes as human-readable."""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size < 1024 or unit == 'GB':
            return f"{size:.2f} {unit}"
        size /= 1024

def get_monitoring_data():
    """Get monitoring data from service or file."""
    if monitoring_available:
        # Get data from running monitoring service
        monitor = MonitoringService.get_instance()
        return monitor.get_monitoring_data()
    else:
        # Try to load data from file
        try:
            monitoring_dir = "monitoring"
            if not os.path.exists(monitoring_dir):
                return None
            
            # Try to load health status
            health_file = os.path.join(monitoring_dir, "health_status.json")
            if os.path.exists(health_file):
                with open(health_file, 'r') as f:
                    health_status = json.load(f)
            else:
                health_status = None
            
            # Try to load database metrics
            db_file = os.path.join(monitoring_dir, "database_metrics.json")
            if os.path.exists(db_file):
                with open(db_file, 'r') as f:
                    db_metrics = json.load(f)
            else:
                db_metrics = None
            
            # Try to load API metrics
            api_file = os.path.join(monitoring_dir, "api_metrics.json")
            if os.path.exists(api_file):
                with open(api_file, 'r') as f:
                    api_metrics = json.load(f)
            else:
                api_metrics = None
            
            # Try to load system metrics
            system_file = os.path.join(monitoring_dir, "system_metrics.json")
            if os.path.exists(system_file):
                with open(system_file, 'r') as f:
                    system_metrics = json.load(f)
            else:
                system_metrics = None
            
            # Try to load resource metrics
            resource_file = os.path.join(monitoring_dir, "resource_metrics.json")
            if os.path.exists(resource_file):
                with open(resource_file, 'r') as f:
                    resource_metrics = json.load(f)
            else:
                resource_metrics = None
            
            # Try to load errors
            errors_file = os.path.join(monitoring_dir, "errors.json")
            if os.path.exists(errors_file):
                with open(errors_file, 'r') as f:
                    errors_data = json.load(f)
                    errors = errors_data.get("errors", [])
            else:
                errors = []
            
            # Try to load alerts
            alerts_file = os.path.join(monitoring_dir, "alerts.json")
            if os.path.exists(alerts_file):
                with open(alerts_file, 'r') as f:
                    alerts_data = json.load(f)
                    alerts = alerts_data.get("alerts", [])
            else:
                alerts = []
            
            # Combine data
            return {
                "health": health_status or {},
                "performance": {
                    "database": db_metrics or {},
                    "api": api_metrics or {},
                    "system": system_metrics or {}
                },
                "resources": resource_metrics or {},
                "errors": errors,
                "alerts": alerts,
                "system_info": {
                    "hostname": socket.gethostname(),
                    "monitoring_service": False
                },
                "last_updated": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            }
        except Exception as e:
            print(f"Error loading monitoring data: {str(e)}")
            return None

def show_status(args):
    """Show current system health status."""
    data = get_monitoring_data()
    if not data:
        print("No monitoring data available. Is the monitoring service running?")
        return
    
    print(colored("\n=== CryoProtect System Health Status ===", Colors.HEADER + Colors.BOLD))
    print(f"Last updated: {data.get('last_updated', 'Unknown')}")
    print()
    
    # Show overall system health
    health = data.get("health", {})
    
    # Database status
    db_status = health.get("database", {}).get("status", "Unknown")
    db_status_color = get_status_color(db_status)
    print(f"Database:    {colored(db_status, db_status_color)}")
    
    # API status
    api_status = health.get("api", {}).get("status", "Unknown")
    api_status_color = get_status_color(api_status)
    print(f"API:         {colored(api_status, api_status_color)}")
    
    # System status
    system_status = health.get("system", {}).get("status", "Unknown")
    system_status_color = get_status_color(system_status)
    print(f"System:      {colored(system_status, system_status_color)}")
    
    print("\n--- Resource Usage ---")
    system_health = health.get("system", {})
    cpu_usage = system_health.get("cpu_usage", 0)
    memory_usage = system_health.get("memory_usage", 0)
    disk_usage = system_health.get("disk_usage", 0)
    
    # Colorize usage values
    cpu_color = Colors.GREEN if cpu_usage < 70 else Colors.YELLOW if cpu_usage < 90 else Colors.RED
    memory_color = Colors.GREEN if memory_usage < 70 else Colors.YELLOW if memory_usage < 90 else Colors.RED
    disk_color = Colors.GREEN if disk_usage < 70 else Colors.YELLOW if disk_usage < 90 else Colors.RED
    
    print(f"CPU:         {colored(f'{cpu_usage:.1f}%', cpu_color)}")
    print(f"Memory:      {colored(f'{memory_usage:.1f}%', memory_color)}")
    print(f"Disk:        {colored(f'{disk_usage:.1f}%', disk_color)}")
    
    # Show database health details
    if args.verbose:
        print("\n--- Database Details ---")
        db_health = health.get("database", {})
        active_conn_type = db_health.get("active_connection_type", "Unknown")
        is_healthy = db_health.get("is_healthy", False)
        last_check_time = db_health.get("last_check_time", 0)
        
        print(f"Connection:  {active_conn_type}")
        print(f"Healthy:     {colored('Yes', Colors.GREEN) if is_healthy else colored('No', Colors.RED)}")
        print(f"Last Check:  {format_timestamp(last_check_time)}")
        
        # Show recent failures if any
        recent_failures = db_health.get("recent_failures", [])
        if recent_failures:
            print("\nRecent Database Failures:")
            for i, failure in enumerate(recent_failures[-5:]):
                timestamp = format_timestamp(failure.get("timestamp", 0))
                error = failure.get("error", "Unknown error")
                print(f"  {i+1}. [{timestamp}] {colored(error, Colors.RED)}")
        
        # Show API endpoints
        print("\n--- API Endpoints ---")
        api_health = health.get("api", {})
        endpoints = api_health.get("endpoints", {})
        
        if endpoints:
            for endpoint, info in endpoints.items():
                status = info.get("status", "Unknown")
                status_color = get_status_color(status)
                response_time = info.get("response_time", 0)
                print(f"  {endpoint}: {colored(status, status_color)} ({response_time:.1f} ms)")
        else:
            print("  No API endpoints monitored")
    
    # Show alert summary
    alerts = data.get("alerts", [])
    if alerts:
        alert_count = len(alerts)
        critical_count = sum(1 for a in alerts if a.get("severity") == "critical")
        error_count = sum(1 for a in alerts if a.get("severity") == "error")
        warning_count = sum(1 for a in alerts if a.get("severity") == "warning")
        
        print("\n--- Alert Summary ---")
        print(f"Total Alerts: {alert_count}")
        if critical_count > 0:
            print(f"Critical:     {colored(critical_count, Colors.RED)}")
        if error_count > 0:
            print(f"Error:        {colored(error_count, Colors.RED)}")
        if warning_count > 0:
            print(f"Warning:      {colored(warning_count, Colors.YELLOW)}")
        
        # Show recent alerts
        print("\nMost Recent Alert:")
        latest = alerts[-1]
        severity = latest.get("severity", "unknown").upper()
        severity_color = get_status_color(severity)
        message = latest.get("message", "No message")
        source = latest.get("source", "unknown")
        timestamp = format_timestamp(latest.get("timestamp", 0))
        
        print(f"  [{timestamp}] {colored(severity, severity_color)} from {source}")
        print(f"  {message}")
    
    print("\nUse 'python monitoring_cli.py metrics' for detailed performance metrics")
    print("Use 'python monitoring_cli.py alerts' to view all recent alerts")
    
    print()

def show_metrics(args):
    """Show performance metrics."""
    data = get_monitoring_data()
    if not data:
        print("No monitoring data available. Is the monitoring service running?")
        return
    
    print(colored("\n=== CryoProtect Performance Metrics ===", Colors.HEADER + Colors.BOLD))
    print(f"Last updated: {data.get('last_updated', 'Unknown')}")
    print()
    
    # Filter metrics by type if specified
    metrics_type = args.type.lower() if args.type else None
    
    # Database metrics
    if not metrics_type or metrics_type == "database":
        db_metrics = data.get("performance", {}).get("database", {})
        if db_metrics:
            print(colored("--- Database Metrics ---", Colors.BOLD))
            
            # Operation counts
            operations = db_metrics.get("operations", {})
            total_ops = operations.get("total", 0)
            success_ops = operations.get("successful", 0)
            failed_ops = operations.get("failed", 0)
            
            if total_ops > 0:
                success_rate = (success_ops / total_ops) * 100
                success_color = Colors.GREEN if success_rate >= 95 else Colors.YELLOW if success_rate >= 80 else Colors.RED
                
                print(f"Operations:    {total_ops} total, {success_ops} successful, {failed_ops} failed")
                print(f"Success Rate:  {colored(f'{success_rate:.1f}%', success_color)}")
            else:
                print("Operations:    None recorded")
            
            # Timing metrics
            timing = db_metrics.get("timing", {})
            avg_time = timing.get("avg_execution_time", 0) * 1000  # Convert to ms
            min_time = timing.get("min_execution_time", float('inf'))
            if min_time != float('inf'):
                min_time *= 1000  # Convert to ms
            max_time = timing.get("max_execution_time", 0) * 1000  # Convert to ms
            
            time_color = Colors.GREEN if avg_time < 50 else Colors.YELLOW if avg_time < 200 else Colors.RED
            
            print(f"Avg Time:      {colored(f'{avg_time:.2f} ms', time_color)}")
            if min_time != float('inf'):
                print(f"Min Time:      {min_time:.2f} ms")
            print(f"Max Time:      {max_time:.2f} ms")
            
            # Database-specific metrics
            db_stats = db_metrics.get("database", {})
            queries = db_stats.get("query_count", 0)
            transactions = db_stats.get("transaction_count", 0)
            rollbacks = db_stats.get("rollback_count", 0)
            conn_attempts = db_stats.get("connection_attempts", 0)
            conn_failures = db_stats.get("connection_failures", 0)
            
            print(f"Queries:       {queries}")
            print(f"Transactions:  {transactions}")
            if rollbacks > 0:
                print(f"Rollbacks:     {colored(rollbacks, Colors.YELLOW)}")
            
            if conn_failures > 0:
                failure_rate = (conn_failures / conn_attempts) * 100 if conn_attempts > 0 else 0
                failure_color = Colors.GREEN if failure_rate < 5 else Colors.YELLOW if failure_rate < 20 else Colors.RED
                print(f"Connections:   {conn_attempts} attempts, {colored(f'{conn_failures} failures ({failure_rate:.1f}%)', failure_color)}")
            else:
                print(f"Connections:   {conn_attempts} attempts, 0 failures")
            
            # Custom metrics
            custom = db_metrics.get("custom_metrics", {})
            if custom and args.verbose:
                print("\nCustom Metrics:")
                for name, value in custom.items():
                    print(f"  {name}: {value}")
            
            print()
    
    # API metrics
    if not metrics_type or metrics_type == "api":
        api_metrics = data.get("performance", {}).get("api", {})
        if api_metrics:
            print(colored("--- API Metrics ---", Colors.BOLD))
            
            # Operation counts
            operations = api_metrics.get("operations", {})
            total_ops = operations.get("total", 0)
            success_ops = operations.get("successful", 0)
            failed_ops = operations.get("failed", 0)
            
            if total_ops > 0:
                success_rate = (success_ops / total_ops) * 100
                success_color = Colors.GREEN if success_rate >= 95 else Colors.YELLOW if success_rate >= 80 else Colors.RED
                
                print(f"Operations:    {total_ops} total, {success_ops} successful, {failed_ops} failed")
                print(f"Success Rate:  {colored(f'{success_rate:.1f}%', success_color)}")
            else:
                print("Operations:    None recorded")
            
            # Timing metrics
            timing = api_metrics.get("timing", {})
            avg_time = timing.get("avg_execution_time", 0) * 1000  # Convert to ms
            min_time = timing.get("min_execution_time", float('inf'))
            if min_time != float('inf'):
                min_time *= 1000  # Convert to ms
            max_time = timing.get("max_execution_time", 0) * 1000  # Convert to ms
            
            time_color = Colors.GREEN if avg_time < 100 else Colors.YELLOW if avg_time < 500 else Colors.RED
            
            print(f"Avg Time:      {colored(f'{avg_time:.2f} ms', time_color)}")
            if min_time != float('inf'):
                print(f"Min Time:      {min_time:.2f} ms")
            print(f"Max Time:      {max_time:.2f} ms")
            
            # Throughput metrics
            throughput = api_metrics.get("throughput", {})
            items_processed = throughput.get("items_processed", 0)
            items_per_second = throughput.get("items_per_second", 0)
            
            if items_processed > 0:
                print(f"Throughput:    {items_processed} items total, {items_per_second:.2f} items/second")
            
            # Custom metrics
            custom = api_metrics.get("custom_metrics", {})
            if custom and args.verbose:
                print("\nCustom Metrics:")
                for name, value in custom.items():
                    print(f"  {name}: {value}")
            
            print()
    
    # System metrics
    if not metrics_type or metrics_type == "system":
        system_metrics = data.get("performance", {}).get("system", {})
        if system_metrics:
            print(colored("--- System Metrics ---", Colors.BOLD))
            
            # Memory usage
            memory = system_metrics.get("memory", {})
            current_memory = memory.get("current_memory_usage", 0)
            peak_memory = memory.get("peak_memory_usage", 0)
            
            print(f"Current Memory: {format_bytes(current_memory)}")
            print(f"Peak Memory:    {format_bytes(peak_memory)}")
            
            # Custom metrics
            custom = system_metrics.get("custom_metrics", {})
            if custom and args.verbose:
                print("\nCustom Metrics:")
                for name, value in custom.items():
                    print(f"  {name}: {value}")
            
            print()
    
    # Resource usage history
    if (not metrics_type or metrics_type == "resources") and args.verbose:
        resources = data.get("resources", {})
        if resources:
            print(colored("--- Resource Usage History ---", Colors.BOLD))
            
            cpu_history = resources.get("cpu_history", [])
            memory_history = resources.get("memory_history", [])
            disk_history = resources.get("disk_history", [])
            timestamps = resources.get("timestamps", [])
            
            if cpu_history and timestamps:
                # Show the last 5 data points
                num_points = min(5, len(cpu_history))
                print(f"Last {num_points} measurements:")
                
                for i in range(-num_points, 0):
                    timestamp = format_timestamp(timestamps[i]) if i < len(timestamps) else "Unknown"
                    cpu = cpu_history[i] if i < len(cpu_history) else "N/A"
                    memory = memory_history[i] if i < len(memory_history) else "N/A"
                    disk = disk_history[i] if i < len(disk_history) else "N/A"
                    
                    cpu_color = Colors.GREEN if cpu < 70 else Colors.YELLOW if cpu < 90 else Colors.RED
                    memory_color = Colors.GREEN if memory < 70 else Colors.YELLOW if memory < 90 else Colors.RED
                    disk_color = Colors.GREEN if disk < 70 else Colors.YELLOW if disk < 90 else Colors.RED
                    
                    print(f"  [{timestamp}] CPU: {colored(f'{cpu:.1f}%', cpu_color)}, "
                          f"Memory: {colored(f'{memory:.1f}%', memory_color)}, "
                          f"Disk: {colored(f'{disk:.1f}%', disk_color)}")
            
            print()
    
    print("Use 'python monitoring_cli.py metrics --type database|api|system|resources' to filter metrics")
    print("Use 'python monitoring_cli.py metrics --verbose' to see detailed metrics")
    
    print()

def show_progress(args):
    """Show progress of ongoing operations."""
    data = get_monitoring_data()
    if not data:
        print("No monitoring data available. Is the monitoring service running?")
        return
    
    # Get progress trackers if available
    if monitoring_available:
        monitor = MonitoringService.get_instance()
        progress_trackers = monitor.get_progress_trackers()
    else:
        # Try to load from file
        try:
            progress_file = os.path.join("monitoring", "progress_trackers.json")
            if os.path.exists(progress_file):
                with open(progress_file, 'r') as f:
                    progress_trackers = json.load(f)
            else:
                progress_trackers = {}
        except Exception:
            progress_trackers = {}
    
    if not progress_trackers:
        print("\nNo active progress trackers found.")
        return
    
    print(colored("\n=== Active Progress Trackers ===", Colors.HEADER + Colors.BOLD))
    print()
    
    for name, tracker in progress_trackers.items():
        total = tracker.get("total_items", 0)
        processed = tracker.get("processed_items", 0)
        successful = tracker.get("successful_items", 0)
        failed = tracker.get("failed_items", 0)
        
        if total > 0:
            progress_pct = (processed / total) * 100
        else:
            progress_pct = 0
        
        # Determine color based on progress and failures
        if failed > 0:
            failure_rate = (failed / processed) * 100 if processed > 0 else 0
            if failure_rate > 20:
                progress_color = Colors.RED
            elif failure_rate > 5:
                progress_color = Colors.YELLOW
            else:
                progress_color = Colors.BLUE
        else:
            progress_color = Colors.GREEN
        
        print(colored(f"Progress Tracker: {name}", Colors.BOLD))
        print(f"Status:       {processed}/{total} items processed ({progress_pct:.1f}%)")
        
        # Create a simple progress bar
        bar_width = 40
        completed_width = int(bar_width * progress_pct / 100)
        progress_bar = "[" + colored("#" * completed_width, progress_color) + " " * (bar_width - completed_width) + "]"
        print(f"Progress:     {progress_bar}")
        
        # Show success/failure stats
        if processed > 0:
            success_rate = (successful / processed) * 100
            success_color = Colors.GREEN if success_rate > 95 else Colors.YELLOW if success_rate > 80 else Colors.RED
            print(f"Results:      {successful} successful, {failed} failed ({colored(f'{success_rate:.1f}%', success_color)} success rate)")
        
        # Show performance stats
        items_per_second = tracker.get("items_per_second", 0)
        est_remaining = tracker.get("estimated_remaining_time", 0)
        est_completion = tracker.get("estimated_completion", None)
        
        print(f"Speed:        {items_per_second:.1f} items/second")
        print(f"Remaining:    {format_duration(est_remaining)}")
        
        if est_completion:
            try:
                completion_time = datetime.fromisoformat(est_completion)
                print(f"ETA:          {completion_time.strftime('%Y-%m-%d %H:%M:%S')}")
            except (ValueError, TypeError):
                print(f"ETA:          {est_completion}")
        
        # Show more details if verbose
        if args.verbose:
            elapsed = tracker.get("elapsed_time", 0)
            avg_batch_time = tracker.get("avg_batch_time", 0) * 1000  # Convert to ms
            avg_batch_size = tracker.get("avg_batch_size", 0)
            
            print(f"Elapsed:      {format_duration(elapsed)}")
            print(f"Batch Avg:    {avg_batch_size:.1f} items, {avg_batch_time:.1f} ms")
        
        print()
    
    print("Use 'python monitoring_cli.py progress --verbose' to see detailed progress information")
    print()

def show_alerts(args):
    """Show recent alerts."""
    data = get_monitoring_data()
    if not data:
        print("No monitoring data available. Is the monitoring service running?")
        return
    
    alerts = data.get("alerts", [])
    if not alerts:
        print("\nNo alerts found.")
        return
    
    print(colored("\n=== Recent Alerts ===", Colors.HEADER + Colors.BOLD))
    print()
    
    # Filter alerts by source or severity if specified
    if args.source:
        source_filter = args.source.lower()
        alerts = [a for a in alerts if a.get("source", "").lower() == source_filter]
    
    if args.severity:
        severity_filter = args.severity.lower()
        alerts = [a for a in alerts if a.get("severity", "").lower() == severity_filter]
    
    # Limit number of alerts to display
    limit = args.limit if args.limit else 10
    alerts = alerts[-limit:]
    
    for alert in alerts:
        severity = alert.get("severity", "unknown").upper()
        severity_color = get_status_color(severity)
        message = alert.get("message", "No message")
        source = alert.get("source", "unknown")
        timestamp = format_timestamp(alert.get("timestamp", 0))
        alert_type = alert.get("alert_type", "unknown")
        
        print(f"[{timestamp}] {colored(severity, severity_color)} - {source}")
        print(f"Type: {alert_type}")
        print(f"Message: {message}")
        
        # Show context if verbose
        if args.verbose and 'context' in alert and alert['context']:
            context_str = json.dumps(alert['context'], indent=2)
            print(f"Context: {context_str}")
        
        print()
    
    if len(data.get("alerts", [])) > limit:
        print(f"Showing {limit} of {len(data.get('alerts', []))} alerts.")
        print(f"Use 'python monitoring_cli.py alerts --limit N' to view more alerts.")
    
    print("Use 'python monitoring_cli.py alerts --verbose' to see alert context")
    print("Use 'python monitoring_cli.py alerts --source SOURCE' to filter by source")
    print("Use 'python monitoring_cli.py alerts --severity SEVERITY' to filter by severity")
    print()

def show_errors(args):
    """Show recent errors."""
    data = get_monitoring_data()
    if not data:
        print("No monitoring data available. Is the monitoring service running?")
        return
    
    errors = data.get("errors", [])
    if not errors:
        print("\nNo errors found.")
        return
    
    print(colored("\n=== Recent Errors ===", Colors.HEADER + Colors.BOLD))
    print()
    
    # Filter errors by type if specified
    if args.type:
        type_filter = args.type.lower()
        errors = [e for e in errors if e.get("error_type", "").lower() == type_filter]
    
    # Limit number of errors to display
    limit = args.limit if args.limit else 10
    errors = errors[-limit:]
    
    for error in errors:
        error_type = error.get("error_type", "unknown")
        message = error.get("message", "No message")
        timestamp = format_timestamp(error.get("timestamp", 0))
        
        print(f"[{timestamp}] {colored(error_type, Colors.RED)}")
        print(f"Message: {message}")
        
        # Show correlation ID if available
        if 'correlation_id' in error:
            print(f"Correlation ID: {error['correlation_id']}")
        
        # Show stack trace if verbose
        if args.verbose and 'stack_trace' in error and error['stack_trace']:
            print(f"\nStack Trace:")
            print(colored(error['stack_trace'], Colors.YELLOW))
        
        print()
    
    if len(data.get("errors", [])) > limit:
        print(f"Showing {limit} of {len(data.get('errors', []))} errors.")
        print(f"Use 'python monitoring_cli.py errors --limit N' to view more errors.")
    
    print("Use 'python monitoring_cli.py errors --verbose' to see stack traces")
    print("Use 'python monitoring_cli.py errors --type ERROR_TYPE' to filter by error type")
    print()

def silence_alert(args):
    """Silence alerts by type."""
    if not monitoring_available:
        print("Monitoring service is not available. Cannot silence alerts.")
        return
    
    if not args.alert_type:
        print("Error: Alert type must be specified.")
        print("Usage: python monitoring_cli.py silence --alert-type ALERT_TYPE [--duration SECONDS]")
        return
    
    monitor = MonitoringService.get_instance()
    duration = args.duration if args.duration else 3600  # Default: 1 hour
    
    monitor.silence_alert(args.alert_type, duration)
    
    print(f"Alert '{args.alert_type}' silenced for {format_duration(duration)}.")
    print(f"Use 'python monitoring_cli.py unsilence --alert-type {args.alert_type}' to unsilence.")

def unsilence_alert(args):
    """Unsilence alerts by type."""
    if not monitoring_available:
        print("Monitoring service is not available. Cannot unsilence alerts.")
        return
    
    if not args.alert_type:
        print("Error: Alert type must be specified.")
        print("Usage: python monitoring_cli.py unsilence --alert-type ALERT_TYPE")
        return
    
    monitor = MonitoringService.get_instance()
    monitor.unsilence_alert(args.alert_type)
    
    print(f"Alert '{args.alert_type}' unsilenced.")

def generate_report(args):
    """Generate monitoring report."""
    data = get_monitoring_data()
    if not data:
        print("No monitoring data available. Is the monitoring service running?")
        return
    
    # Determine report format
    format_type = args.format.lower() if args.format else "text"
    
    if format_type == "json":
        # Generate JSON report
        report_file = args.output if args.output else "monitoring_report.json"
        with open(report_file, 'w') as f:
            json.dump(data, f, indent=2)
        print(f"Report saved to {report_file}")
        return
    
    # Generate text report
    report = []
    
    # Report header
    report.append(f"CryoProtect Monitoring Report")
    report.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report.append(f"Hostname: {data.get('system_info', {}).get('hostname', 'Unknown')}")
    report.append("")
    
    # System Health
    report.append("=== System Health ===")
    health = data.get("health", {})
    
    # Database status
    db_status = health.get("database", {}).get("status", "Unknown")
    report.append(f"Database: {db_status}")
    
    # API status
    api_status = health.get("api", {}).get("status", "Unknown")
    report.append(f"API: {api_status}")
    
    # System status
    system_status = health.get("system", {}).get("status", "Unknown")
    report.append(f"System: {system_status}")
    
    # Resource Usage
    system_health = health.get("system", {})
    cpu_usage = system_health.get("cpu_usage", 0)
    memory_usage = system_health.get("memory_usage", 0)
    disk_usage = system_health.get("disk_usage", 0)
    
    report.append("")
    report.append("=== Resource Usage ===")
    report.append(f"CPU: {cpu_usage:.1f}%")
    report.append(f"Memory: {memory_usage:.1f}%")
    report.append(f"Disk: {disk_usage:.1f}%")
    
    # Performance Metrics
    report.append("")
    report.append("=== Performance Metrics ===")
    
    # Database metrics
    db_metrics = data.get("performance", {}).get("database", {})
    if db_metrics:
        report.append("Database:")
        operations = db_metrics.get("operations", {})
        total_ops = operations.get("total", 0)
        success_ops = operations.get("successful", 0)
        failed_ops = operations.get("failed", 0)
        
        if total_ops > 0:
            success_rate = (success_ops / total_ops) * 100
            report.append(f"  Operations: {total_ops} total, {success_ops} successful, {failed_ops} failed")
            report.append(f"  Success Rate: {success_rate:.1f}%")
        else:
            report.append("  Operations: None recorded")
        
        timing = db_metrics.get("timing", {})
        avg_time = timing.get("avg_execution_time", 0) * 1000  # Convert to ms
        report.append(f"  Avg Time: {avg_time:.2f} ms")
        
        db_stats = db_metrics.get("database", {})
        queries = db_stats.get("query_count", 0)
        transactions = db_stats.get("transaction_count", 0)
        rollbacks = db_stats.get("rollback_count", 0)
        report.append(f"  Queries: {queries}")
        report.append(f"  Transactions: {transactions}")
        report.append(f"  Rollbacks: {rollbacks}")
    
    # API metrics
    api_metrics = data.get("performance", {}).get("api", {})
    if api_metrics:
        report.append("API:")
        operations = api_metrics.get("operations", {})
        total_ops = operations.get("total", 0)
        success_ops = operations.get("successful", 0)
        failed_ops = operations.get("failed", 0)
        
        if total_ops > 0:
            success_rate = (success_ops / total_ops) * 100
            report.append(f"  Operations: {total_ops} total, {success_ops} successful, {failed_ops} failed")
            report.append(f"  Success Rate: {success_rate:.1f}%")
        else:
            report.append("  Operations: None recorded")
        
        timing = api_metrics.get("timing", {})
        avg_time = timing.get("avg_execution_time", 0) * 1000  # Convert to ms
        report.append(f"  Avg Time: {avg_time:.2f} ms")
    
    # System metrics
    system_metrics = data.get("performance", {}).get("system", {})
    if system_metrics:
        report.append("System:")
        memory = system_metrics.get("memory", {})
        current_memory = memory.get("current_memory_usage", 0)
        peak_memory = memory.get("peak_memory_usage", 0)
        
        report.append(f"  Current Memory: {format_bytes(current_memory)}")
        report.append(f"  Peak Memory: {format_bytes(peak_memory)}")
    
    # Progress Trackers
    if monitoring_available:
        monitor = MonitoringService.get_instance()
        progress_trackers = monitor.get_progress_trackers()
        
        if progress_trackers:
            report.append("")
            report.append("=== Active Progress Trackers ===")
            
            for name, tracker in progress_trackers.items():
                total = tracker.get("total_items", 0)
                processed = tracker.get("processed_items", 0)
                successful = tracker.get("successful_items", 0)
                failed = tracker.get("failed_items", 0)
                
                if total > 0:
                    progress_pct = (processed / total) * 100
                else:
                    progress_pct = 0
                
                items_per_second = tracker.get("items_per_second", 0)
                est_remaining = tracker.get("estimated_remaining_time", 0)
                
                report.append(f"{name}:")
                report.append(f"  Progress: {processed}/{total} items ({progress_pct:.1f}%)")
                report.append(f"  Results: {successful} successful, {failed} failed")
                report.append(f"  Speed: {items_per_second:.1f} items/second")
                report.append(f"  Remaining: {format_duration(est_remaining)}")
    
    # Recent Alerts
    alerts = data.get("alerts", [])
    if alerts:
        report.append("")
        report.append("=== Recent Alerts ===")
        
        # Show last 5 alerts
        for alert in alerts[-5:]:
            severity = alert.get("severity", "unknown").upper()
            message = alert.get("message", "No message")
            source = alert.get("source", "unknown")
            timestamp = format_timestamp(alert.get("timestamp", 0))
            
            report.append(f"[{timestamp}] {severity} - {source}")
            report.append(f"  {message}")
    
    # Recent Errors
    errors = data.get("errors", [])
    if errors:
        report.append("")
        report.append("=== Recent Errors ===")
        
        # Show last 5 errors
        for error in errors[-5:]:
            error_type = error.get("error_type", "unknown")
            message = error.get("message", "No message")
            timestamp = format_timestamp(error.get("timestamp", 0))
            
            report.append(f"[{timestamp}] {error_type}")
            report.append(f"  {message}")
    
    # Save report to file if output specified
    if args.output:
        with open(args.output, 'w') as f:
            f.write("\n".join(report))
        print(f"Report saved to {args.output}")
    else:
        print("\n".join(report))

def open_dashboard(args):
    """Open the monitoring dashboard in a browser."""
    port = args.port if args.port else 5001
    host = args.host if args.host else "localhost"
    
    url = f"http://{host}:{port}/monitoring"
    
    print(f"Opening dashboard at {url}")
    webbrowser.open(url)

def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="CryoProtect Monitoring CLI")
    subparsers = parser.add_subparsers(dest="command", help="Command to execute")
    
    # Status command
    status_parser = subparsers.add_parser("status", help="Show current system health status")
    status_parser.add_argument("-v", "--verbose", action="store_true", help="Show verbose output")
    
    # Metrics command
    metrics_parser = subparsers.add_parser("metrics", help="Show performance metrics")
    metrics_parser.add_argument("-t", "--type", choices=["database", "api", "system", "resources"],
                              help="Type of metrics to show")
    metrics_parser.add_argument("-v", "--verbose", action="store_true", help="Show verbose output")
    
    # Progress command
    progress_parser = subparsers.add_parser("progress", help="Show progress of ongoing operations")
    progress_parser.add_argument("-v", "--verbose", action="store_true", help="Show verbose output")
    
    # Alerts command
    alerts_parser = subparsers.add_parser("alerts", help="Show recent alerts")
    alerts_parser.add_argument("-s", "--source", help="Filter alerts by source")
    alerts_parser.add_argument("--severity", help="Filter alerts by severity")
    alerts_parser.add_argument("-l", "--limit", type=int, help="Limit number of alerts to show")
    alerts_parser.add_argument("-v", "--verbose", action="store_true", help="Show verbose output")
    
    # Errors command
    errors_parser = subparsers.add_parser("errors", help="Show recent errors")
    errors_parser.add_argument("-t", "--type", help="Filter errors by type")
    errors_parser.add_argument("-l", "--limit", type=int, help="Limit number of errors to show")
    errors_parser.add_argument("-v", "--verbose", action="store_true", help="Show verbose output")
    
    # Silence command
    silence_parser = subparsers.add_parser("silence", help="Silence alerts by type")
    silence_parser.add_argument("--alert-type", required=True, help="Alert type to silence")
    silence_parser.add_argument("-d", "--duration", type=int, help="Duration to silence for in seconds")
    
    # Unsilence command
    unsilence_parser = subparsers.add_parser("unsilence", help="Unsilence alerts by type")
    unsilence_parser.add_argument("--alert-type", required=True, help="Alert type to unsilence")
    
    # Report command
    report_parser = subparsers.add_parser("report", help="Generate monitoring report")
    report_parser.add_argument("-o", "--output", help="Output file path")
    report_parser.add_argument("-f", "--format", choices=["text", "json"], default="text",
                             help="Report format")
    
    # Dashboard command
    dashboard_parser = subparsers.add_parser("dashboard", help="Open the monitoring dashboard in a browser")
    dashboard_parser.add_argument("-p", "--port", type=int, help="Dashboard port")
    dashboard_parser.add_argument("--host", help="Dashboard host")
    
    args = parser.parse_args()
    
    # Execute appropriate command
    if args.command == "status":
        show_status(args)
    elif args.command == "metrics":
        show_metrics(args)
    elif args.command == "progress":
        show_progress(args)
    elif args.command == "alerts":
        show_alerts(args)
    elif args.command == "errors":
        show_errors(args)
    elif args.command == "silence":
        silence_alert(args)
    elif args.command == "unsilence":
        unsilence_alert(args)
    elif args.command == "report":
        generate_report(args)
    elif args.command == "dashboard":
        open_dashboard(args)
    else:
        # Show help if no command specified
        parser.print_help()

if __name__ == "__main__":
    main()