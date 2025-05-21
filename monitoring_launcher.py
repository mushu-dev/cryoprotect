#!/usr/bin/env python3
"""
Monitoring Launcher for CryoProtect

This script serves as the main entry point for the CryoProtect monitoring system.
It:
1. Loads configuration from file, environment variables, and command-line arguments
2. Sets up logging based on configuration
3. Initializes the monitoring service
4. Sets up notification integrations
5. Starts the monitoring service and dashboard
6. Provides a graceful shutdown mechanism

Usage:
    python monitoring_launcher.py [--config CONFIG_FILE] [--dashboard-port PORT]
"""

import os
import sys
import time
import signal
import logging
import argparse
import threading
import webbrowser
import json
import socket
from logging.handlers import RotatingFileHandler
from typing import Dict, Any, Optional, List, Union
from datetime import datetime, timedelta

# Import configuration loader
from monitoring_config_loader import MonitoringConfig

# Import monitoring service
try:
    from unified_monitoring import MonitoringService, start_monitoring
    from monitoring_integrations import create_notifiers_from_config
    unified_monitoring_available = True
except ImportError:
    unified_monitoring_available = False

# Import monitoring utilities if available
try:
    from monitoring_utils import ConnectionHealthMonitor, create_monitoring_directory
    monitoring_utils_available = True
except ImportError:
    monitoring_utils_available = False

# Configure default logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('monitoring_launcher')

# Global variables
monitoring_service = None
database_monitor = None
should_exit = threading.Event()

def setup_logging(config):
    """
    Set up logging based on configuration.
    
    Args:
        config: Configuration manager instance
    """
    # Get logging configuration
    log_level_name = config.get("logging.level", "INFO")
    log_file = config.get("logging.file")
    log_format = config.get("logging.format", 
                          '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    max_size_mb = config.get("logging.max_size_mb", 10)
    backup_count = config.get("logging.backup_count", 5)
    
    # Convert level name to level
    log_level = getattr(logging, log_level_name.upper(), logging.INFO)
    
    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(log_level)
    
    # Clear existing handlers
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
    
    # Create formatter
    formatter = logging.Formatter(log_format)
    
    # Add console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)
    
    # Add file handler if configured
    if log_file:
        # Create directory if it doesn't exist
        log_dir = os.path.dirname(log_file)
        if log_dir and not os.path.exists(log_dir):
            os.makedirs(log_dir)
        
        file_handler = RotatingFileHandler(
            log_file,
            maxBytes=max_size_mb * 1024 * 1024,
            backupCount=backup_count
        )
        file_handler.setFormatter(formatter)
        root_logger.addHandler(file_handler)
    
    logger.info(f"Logging configured with level {log_level_name}")
    if log_file:
        logger.info(f"Log file: {log_file}")

def setup_monitoring(config, flask_app=None):
    """
    Set up and start the monitoring service.
    
    Args:
        config: Configuration manager instance
        flask_app: Optional Flask application to monitor
        
    Returns:
        MonitoringService: Initialized monitoring service
    """
    global monitoring_service, database_monitor
    
    if not unified_monitoring_available:
        logger.error("Unified monitoring module not available.")
        return None
    
    # Get configuration values
    dashboard_port = config.get("monitoring.dashboard_port", 5001)
    dashboard_enabled = config.get("monitoring.dashboard_enabled", True)
    
    # Initialize monitoring service
    monitoring_service = start_monitoring(
        app=flask_app,
        dashboard_port=dashboard_port if dashboard_enabled else None
    )
    
    # Update configuration
    monitoring_service.config.update(config.get("monitoring", {}))
    
    # Set up notification integrations
    setup_notifications(config, monitoring_service)
    
    # Open dashboard in browser if configured
    if dashboard_enabled and config.get("monitoring.open_dashboard_on_start", False):
        open_dashboard(dashboard_port)

    # Initialize database connection monitoring if available
    if monitoring_utils_available:
        # Create directory for monitoring data
        create_monitoring_directory()
        
        # Set up database connection monitoring
        db_monitoring_interval = config.get("monitoring.health_check_intervals.database", 60)
        database_monitor = ConnectionHealthMonitor.get_instance()
        database_monitor.health_check_interval = db_monitoring_interval
        database_monitor.start_monitoring()
        logger.info(f"Database connection monitoring started with interval {db_monitoring_interval}s")
    
    logger.info(f"Monitoring service started")
    if dashboard_enabled:
        logger.info(f"Dashboard available at http://localhost:{dashboard_port}/monitoring")
    
    return monitoring_service

def setup_notifications(config, monitoring_service):
    """
    Set up notification integrations.
    
    Args:
        config: Configuration manager instance
        monitoring_service: Monitoring service instance
    """
    # Load notification configuration
    notifications_config = config.get("notifications", {})
    
    # Create notifiers from configuration
    try:
        notifiers = create_notifiers_from_config(notifications_config)
        
        # Register notifiers as alert handlers
        for notifier in notifiers:
            monitoring_service.add_alert_handler(notifier.handle_alert)
            logger.info(f"Registered alert handler: {notifier.name}")
    except Exception as e:
        logger.error(f"Error setting up notifications: {str(e)}")

def open_dashboard(port):
    """
    Open monitoring dashboard in web browser.
    
    Args:
        port: Dashboard port
    """
    url = f"http://localhost:{port}/monitoring"
    
    try:
        # Wait a bit for the service to start
        time.sleep(2)
        webbrowser.open(url)
    except Exception as e:
        logger.error(f"Error opening dashboard in browser: {str(e)}")

def signal_handler(sig, frame):
    """
    Handle interrupt signals.
    
    Args:
        sig: Signal number
        frame: Current stack frame
    """
    global should_exit
    
    logger.info(f"Received signal {sig}, shutting down...")
    should_exit.set()

def shutdown_monitoring():
    """Shut down the monitoring service."""
    global monitoring_service, database_monitor
    
    # Stop database connection monitoring
    if monitoring_utils_available and database_monitor:
        logger.info("Stopping database connection monitoring...")
        database_monitor.stop_monitoring()
        logger.info("Database connection monitoring stopped.")
    
    # Stop monitoring service
    if monitoring_service:
        logger.info("Stopping monitoring service...")
        monitoring_service.stop()
        logger.info("Monitoring service stopped.")

def save_system_info():
    """Save system information for diagnostics."""
    try:
        # Create monitoring directory if it doesn't exist
        os.makedirs("monitoring", exist_ok=True)
        
        # Gather system information
        system_info = {
            "hostname": socket.gethostname(),
            "platform": sys.platform,
            "python_version": sys.version,
            "timestamp": datetime.now().isoformat(),
            "monitoring_modules": {
                "unified_monitoring": unified_monitoring_available,
                "monitoring_utils": monitoring_utils_available,
            }
        }
        
        # Save to file
        with open("monitoring/system_info.json", 'w') as f:
            json.dump(system_info, f, indent=2)
            
        logger.info("System information saved to monitoring/system_info.json")
    except Exception as e:
        logger.warning(f"Could not save system information: {str(e)}")

def setup_prometheus_exporter(config):
    """
    Set up Prometheus metrics exporter if configured.
    
    Args:
        config: Configuration manager instance
    """
    if not config.is_enabled("prometheus"):
        return
    
    prometheus_port = config.get("prometheus.port", 9090)
    
    try:
        from prometheus_client import start_http_server, Gauge, Counter, Summary
        
        # Start Prometheus HTTP server
        start_http_server(prometheus_port)
        logger.info(f"Prometheus metrics available at http://localhost:{prometheus_port}/metrics")
        
        # Create metrics
        HEALTH_STATUS = Gauge('cryoprotect_health_status', 
                             'Health status of CryoProtect components', 
                             ['component'])
        
        CPU_USAGE = Gauge('cryoprotect_cpu_usage', 'CPU usage percentage')
        MEMORY_USAGE = Gauge('cryoprotect_memory_usage', 'Memory usage percentage')
        DISK_USAGE = Gauge('cryoprotect_disk_usage', 'Disk usage percentage')
        
        DB_OPERATIONS = Counter('cryoprotect_db_operations_total', 
                               'Total database operations',
                               ['operation_type', 'status'])
        
        API_REQUEST_DURATION = Summary('cryoprotect_api_request_duration_seconds', 
                                    'API request duration in seconds',
                                    ['endpoint'])
        
        # Start metrics collection thread
        def metrics_collector():
            while not should_exit.is_set():
                if monitoring_service:
                    data = monitoring_service.get_monitoring_data()
                    
                    # Update health metrics
                    health = data.get("health", {})
                    
                    for component in ["database", "api", "system"]:
                        status = health.get(component, {}).get("status", "Unknown")
                        status_value = 1 if status == "OK" else 0
                        HEALTH_STATUS.labels(component=component).set(status_value)
                    
                    # Update resource usage metrics
                    system_health = health.get("system", {})
                    CPU_USAGE.set(system_health.get("cpu_usage", 0))
                    MEMORY_USAGE.set(system_health.get("memory_usage", 0))
                    DISK_USAGE.set(system_health.get("disk_usage", 0))
                    
                    # Wait for next collection interval
                    time.sleep(15)
        
        # Start collector thread
        collector_thread = threading.Thread(
            target=metrics_collector, 
            daemon=True
        )
        collector_thread.start()
        logger.info("Prometheus metrics collection started")
        
    except ImportError:
        logger.warning("Prometheus client library not available. Prometheus exporter disabled.")
    except Exception as e:
        logger.error(f"Error setting up Prometheus exporter: {str(e)}")

def main():
    """Main entry point."""
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="CryoProtect Monitoring Launcher")
    parser.add_argument("--config", help="Path to configuration file")
    parser.add_argument("--dashboard-port", type=int, help="Dashboard port")
    parser.add_argument("--no-dashboard", action="store_true", help="Disable dashboard")
    parser.add_argument("--open-browser", action="store_true", help="Open dashboard in browser")
    parser.add_argument("--system-info", action="store_true", help="Save system information and exit")
    parser.add_argument("--print-config", action="store_true", help="Print configuration and exit")
    
    args, _ = parser.parse_known_args()
    
    # Save system information if requested
    if args.system_info:
        save_system_info()
        return 0
    
    # Load configuration
    config = MonitoringConfig.get_instance(config_file=args.config)
    
    # Print configuration if requested
    if args.print_config:
        print(json.dumps(config.get_all(), indent=2))
        return 0
    
    # Override configuration with command-line arguments
    if args.dashboard_port:
        config.set("monitoring.dashboard_port", args.dashboard_port)
    
    if args.no_dashboard:
        config.set("monitoring.dashboard_enabled", False)
    
    if args.open_browser:
        config.set("monitoring.open_dashboard_on_start", True)
    
    # Set up logging
    setup_logging(config)
    
    # Check if monitoring is enabled
    if not config.is_enabled("monitoring"):
        logger.info("Monitoring is disabled in configuration. Exiting.")
        return 0
    
    # Set up signal handlers
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)
    
    # Set up Prometheus exporter if configured
    setup_prometheus_exporter(config)
    
    # Set up and start monitoring
    app = None
    try:
        from app import app
        logger.info("Found Flask application, integrating with monitoring")
    except ImportError:
        logger.info("No Flask application found, using standalone monitoring")
    
    monitoring = setup_monitoring(config, app)
    
    if not monitoring:
        logger.error("Failed to start monitoring service. Exiting.")
        return 1
    
    # Main loop - keep running until signaled to exit
    logger.info("Monitoring service running. Press Ctrl+C to exit.")
    try:
        while not should_exit.is_set():
            time.sleep(1)
    except KeyboardInterrupt:
        logger.info("Keyboard interrupt received, shutting down...")
    finally:
        shutdown_monitoring()
    
    return 0

if __name__ == "__main__":
    sys.exit(main())