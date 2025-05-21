#!/usr/bin/env python3
"""
CryoProtect Production Test Monitor

This script monitors the real-time metrics and logs during production testing.
It displays a live dashboard of container resource usage, API endpoint performance,
and test progress.
"""

import os
import sys
import time
import json
import signal
import subprocess
import threading
import datetime
from typing import Dict, List, Any, Optional

# Try to import rich for nice terminal formatting
try:
    from rich.console import Console
    from rich.panel import Panel
    from rich.layout import Layout
    from rich.table import Table
    from rich.live import Live
    from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TimeElapsedColumn
    from rich.text import Text
    RICH_AVAILABLE = True
except ImportError:
    RICH_AVAILABLE = False

class ProductionTestMonitor:
    """Monitors production tests and displays real-time metrics."""
    
    def __init__(self):
        """Initialize the monitor."""
        self.running = True
        self.containers = ["cryoprotect-app", "cryoprotect-rdkit", "cryoprotect-test"]
        self.container_stats = {}
        self.api_metrics = {}
        self.test_progress = {
            "status": "Not Started",
            "total": 0,
            "passed": 0,
            "failed": 0,
            "skipped": 0,
            "current_test": "None"
        }
        self.logs = []
        self.max_logs = 20
        
        if RICH_AVAILABLE:
            self.console = Console()
            self.layout = Layout()
            self.setup_layout()
    
    def setup_layout(self):
        """Set up the rich layout."""
        if not RICH_AVAILABLE:
            return
            
        self.layout.split(
            Layout(name="header", size=3),
            Layout(name="main", ratio=1),
            Layout(name="footer", size=3)
        )
        
        self.layout["main"].split_row(
            Layout(name="containers", ratio=1),
            Layout(name="tests", ratio=1)
        )
        
        self.layout["tests"].split(
            Layout(name="progress", size=10),
            Layout(name="logs", ratio=1)
        )
        
    def start(self):
        """Start the monitoring threads."""
        # Register signal handler for clean shutdown
        signal.signal(signal.SIGINT, self.handle_sigint)
        
        # Start monitoring threads
        threading.Thread(target=self.monitor_containers, daemon=True).start()
        threading.Thread(target=self.monitor_api, daemon=True).start()
        threading.Thread(target=self.monitor_tests, daemon=True).start()
        threading.Thread(target=self.monitor_logs, daemon=True).start()
        
        # Start the display
        if RICH_AVAILABLE:
            with Live(self.layout, refresh_per_second=1, console=self.console) as live:
                while self.running:
                    self.update_layout()
                    time.sleep(0.5)
        else:
            # Fallback to simple text display
            while self.running:
                self.display_simple()
                time.sleep(1)
                os.system('clear' if os.name == 'posix' else 'cls')
    
    def handle_sigint(self, sig, frame):
        """Handle Ctrl+C gracefully."""
        print("\nShutting down monitor...")
        self.running = False
        time.sleep(1)
        sys.exit(0)
    
    def monitor_containers(self):
        """Monitor container resource usage."""
        while self.running:
            for container in self.containers:
                try:
                    # Get stats using podman stats
                    cmd = ["podman", "stats", container, "--no-stream", "--format", "json"]
                    process = subprocess.run(cmd, capture_output=True, text=True, check=False)
                    
                    if process.returncode == 0 and process.stdout:
                        try:
                            stats = json.loads(process.stdout)
                            if stats and len(stats) > 0:
                                self.container_stats[container] = {
                                    "cpu": stats[0].get("cpu_percent", "N/A"),
                                    "memory": stats[0].get("mem_usage", "N/A"),
                                    "memory_percent": stats[0].get("mem_percent", "N/A"),
                                    "status": "Running"
                                }
                            else:
                                self.container_stats[container] = {
                                    "status": "No Data",
                                    "cpu": "N/A",
                                    "memory": "N/A",
                                    "memory_percent": "N/A"
                                }
                        except json.JSONDecodeError:
                            self.container_stats[container] = {
                                "status": "Parse Error",
                                "cpu": "N/A",
                                "memory": "N/A",
                                "memory_percent": "N/A"
                            }
                    else:
                        self.container_stats[container] = {
                            "status": "Not Running",
                            "cpu": "N/A",
                            "memory": "N/A",
                            "memory_percent": "N/A"
                        }
                except Exception as e:
                    self.container_stats[container] = {
                        "status": f"Error: {str(e)}",
                        "cpu": "N/A",
                        "memory": "N/A",
                        "memory_percent": "N/A"
                    }
            
            time.sleep(2)
    
    def monitor_api(self):
        """Monitor API endpoint performance."""
        endpoints = [
            {"url": "http://localhost:5001/health", "name": "App Health"},
            {"url": "http://localhost:5002/health", "name": "RDKit Health"},
            {"url": "http://localhost:5001/api/v1/rdkit/check", "name": "RDKit Integration"},
            {"url": "http://localhost:5001/molecule/CCO", "name": "Molecule Properties"}
        ]
        
        import requests
        from requests.exceptions import RequestException
        
        while self.running:
            for endpoint in endpoints:
                try:
                    start_time = time.time()
                    response = requests.get(endpoint["url"], timeout=5)
                    end_time = time.time()
                    
                    if response.status_code == 200:
                        self.api_metrics[endpoint["name"]] = {
                            "status": "UP",
                            "response_time": round((end_time - start_time) * 1000, 2),  # ms
                            "status_code": response.status_code
                        }
                    else:
                        self.api_metrics[endpoint["name"]] = {
                            "status": "ERROR",
                            "response_time": round((end_time - start_time) * 1000, 2),  # ms
                            "status_code": response.status_code
                        }
                except RequestException as e:
                    self.api_metrics[endpoint["name"]] = {
                        "status": "DOWN",
                        "response_time": 0,
                        "status_code": "N/A",
                        "error": str(e)
                    }
                except Exception as e:
                    self.api_metrics[endpoint["name"]] = {
                        "status": "ERROR",
                        "response_time": 0,
                        "status_code": "N/A",
                        "error": str(e)
                    }
            
            time.sleep(5)
    
    def monitor_tests(self):
        """Monitor test progress."""
        # Check for real_data_test_results.json file
        while self.running:
            try:
                if os.path.exists("real_data_test_results.json"):
                    with open("real_data_test_results.json", "r") as f:
                        results = json.load(f)
                    
                    self.test_progress["status"] = results.get("status", "Unknown")
                    self.test_progress["total"] = results.get("total_tests", 0)
                    self.test_progress["passed"] = results.get("passed_tests", 0)
                    self.test_progress["failed"] = results.get("failed_tests", 0)
                    self.test_progress["skipped"] = results.get("skipped_tests", 0)
                    
                    # Get the most recent test case
                    if "test_cases" in results and results["test_cases"]:
                        self.test_progress["current_test"] = results["test_cases"][-1]["name"]
            except Exception as e:
                self.test_progress["status"] = f"Error: {str(e)}"
            
            time.sleep(3)
    
    def monitor_logs(self):
        """Monitor log files."""
        log_file = "real_data_test.log"
        
        while self.running:
            try:
                if os.path.exists(log_file):
                    with open(log_file, "r") as f:
                        lines = f.readlines()
                    
                    # Get the last N lines
                    recent_logs = lines[-self.max_logs:]
                    self.logs = [line.strip() for line in recent_logs]
            except Exception as e:
                self.logs.append(f"Error reading logs: {str(e)}")
            
            time.sleep(2)
    
    def update_layout(self):
        """Update the rich layout with current data."""
        if not RICH_AVAILABLE:
            return
        
        # Header
        header = Panel(
            Text(f"CryoProtect Production Test Monitor - {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", style="bold"),
            style="blue"
        )
        self.layout["header"].update(header)
        
        # Container stats
        container_table = Table(title="Container Status")
        container_table.add_column("Container")
        container_table.add_column("Status")
        container_table.add_column("CPU %")
        container_table.add_column("Memory")
        container_table.add_column("Memory %")
        
        for container, stats in self.container_stats.items():
            status_style = "green" if stats["status"] == "Running" else "red"
            container_table.add_row(
                container,
                Text(stats["status"], style=status_style),
                str(stats["cpu"]),
                str(stats["memory"]),
                str(stats["memory_percent"])
            )
        
        # API metrics
        api_table = Table(title="API Endpoints")
        api_table.add_column("Endpoint")
        api_table.add_column("Status")
        api_table.add_column("Response Time (ms)")
        api_table.add_column("Status Code")
        
        for endpoint, metrics in self.api_metrics.items():
            status_style = "green" if metrics["status"] == "UP" else "red"
            api_table.add_row(
                endpoint,
                Text(metrics["status"], style=status_style),
                str(metrics["response_time"]),
                str(metrics["status_code"])
            )
        
        # Combine container and API tables
        all_metrics = Table.grid()
        all_metrics.add_row(container_table)
        all_metrics.add_row(api_table)
        self.layout["containers"].update(Panel(all_metrics, title="System Metrics"))
        
        # Test progress
        status_style = {
            "Passed": "green",
            "Failed": "red",
            "Running": "yellow",
            "Not Started": "blue"
        }.get(self.test_progress["status"], "white")
        
        progress_table = Table.grid()
        progress_table.add_row(Text(f"Status: {self.test_progress['status']}", style=status_style))
        progress_table.add_row(Text(f"Current Test: {self.test_progress['current_test']}"))
        progress_table.add_row(Text(f"Total: {self.test_progress['total']} | Passed: {self.test_progress['passed']} | Failed: {self.test_progress['failed']} | Skipped: {self.test_progress['skipped']}"))
        
        progress = Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TimeElapsedColumn()
        )
        
        with progress:
            completed = self.test_progress["passed"] + self.test_progress["failed"] + self.test_progress["skipped"]
            total = max(1, self.test_progress["total"]) if self.test_progress["total"] > 0 else 100
            task = progress.add_task("Test Progress", total=total, completed=completed)
            progress_display = progress._live.renderable
        
        progress_panel = Panel(
            progress_table,
            title="Test Progress",
            style="green" if self.test_progress["status"] == "Passed" else "red" if self.test_progress["status"] == "Failed" else "yellow"
        )
        self.layout["progress"].update(progress_panel)
        
        # Logs
        log_panel = Panel("\n".join(self.logs), title="Recent Logs")
        self.layout["logs"].update(log_panel)
        
        # Footer
        footer = Panel(
            Text("Press Ctrl+C to exit", style="bold"),
            style="blue"
        )
        self.layout["footer"].update(footer)
    
    def display_simple(self):
        """Display a simple text-based dashboard."""
        print(f"=== CryoProtect Production Test Monitor - {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} ===\n")
        
        # Container stats
        print("=== Container Status ===")
        for container, stats in self.container_stats.items():
            print(f"{container}: {stats['status']} | CPU: {stats['cpu']} | Memory: {stats['memory']} ({stats['memory_percent']})")
        print()
        
        # API metrics
        print("=== API Endpoints ===")
        for endpoint, metrics in self.api_metrics.items():
            print(f"{endpoint}: {metrics['status']} | Response Time: {metrics['response_time']}ms | Status Code: {metrics['status_code']}")
        print()
        
        # Test progress
        print("=== Test Progress ===")
        print(f"Status: {self.test_progress['status']}")
        print(f"Current Test: {self.test_progress['current_test']}")
        print(f"Total: {self.test_progress['total']} | Passed: {self.test_progress['passed']} | Failed: {self.test_progress['failed']} | Skipped: {self.test_progress['skipped']}")
        print()
        
        # Recent logs
        print("=== Recent Logs ===")
        for log in self.logs[-5:]:
            print(log)

def main():
    """Main function."""
    try:
        # Check if rich is installed, otherwise install it
        if not RICH_AVAILABLE:
            print("Rich library not found. Installing...")
            subprocess.run([sys.executable, "-m", "pip", "install", "rich"], check=True)
            print("Rich installed. Please restart the script.")
            sys.exit(0)
        
        print("Starting CryoProtect Production Test Monitor...")
        monitor = ProductionTestMonitor()
        monitor.start()
    except KeyboardInterrupt:
        print("\nMonitor stopped.")
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()