#!/usr/bin/env python3
"""
Monitoring Utilities for CryoProtect v2 Database Operations

This module provides monitoring capabilities for database operations:
1. Connection health monitoring with automatic fallback
2. Progress tracking with time estimation
3. Performance metrics collection
4. Centralized logging and reporting

Based on specifications in DATABASE_POPULATION_ISSUES.md (Section 5.3)
"""

import os
import time
import json
import logging
import threading
import psutil
import socket
from typing import Dict, Any, List, Optional, Union, Callable, Tuple
from datetime import datetime, timedelta
from contextlib import contextmanager
import uuid

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('monitoring_utils')

# Try to import db_connection_utils, handle import error gracefully
try:
    from db_connection_utils import ConnectionManager, safe_transaction
except ImportError:
    logger.warning("Could not import ConnectionManager from db_connection_utils")
    ConnectionManager = None
    safe_transaction = None


class PerformanceMetrics:
    """
    Collects and stores performance metrics for database operations.
    
    This class provides methods to track various performance metrics such as
    query execution time, memory usage, and throughput. It also provides
    methods to generate reports and alerts based on these metrics.
    """
    
    def __init__(self, name: str, metrics_file: Optional[str] = None):
        """
        Initialize performance metrics collector.
        
        Args:
            name: Name of the metrics collector (for identification)
            metrics_file: Optional path to file for persisting metrics
        """
        self.name = name
        self.metrics_file = metrics_file
        self.start_time = time.time()
        self.lock = threading.RLock()
        self.metrics = {
            "name": name,
            "start_time": self.start_time,
            "operations": {
                "total": 0,
                "successful": 0,
                "failed": 0
            },
            "timing": {
                "total_execution_time": 0,
                "min_execution_time": float('inf'),
                "max_execution_time": 0,
                "avg_execution_time": 0
            },
            "throughput": {
                "items_processed": 0,
                "items_per_second": 0,
                "bytes_processed": 0,
                "bytes_per_second": 0
            },
            "memory": {
                "peak_memory_usage": 0,
                "current_memory_usage": 0
            },
            "database": {
                "connection_attempts": 0,
                "connection_failures": 0,
                "query_count": 0,
                "transaction_count": 0,
                "rollback_count": 0
            },
            "custom_metrics": {},
            "history": []
        }
# Load existing metrics if file exists
        if metrics_file and os.path.exists(metrics_file):
            try:
                with open(metrics_file, 'r') as f:
                    saved_metrics = json.load(f)
                    # Merge saved metrics with current metrics
                    self._merge_metrics(saved_metrics)
                    logger.info(f"Loaded metrics from {metrics_file}")
            except Exception as e:
                logger.warning(f"Could not load metrics from {metrics_file}: {str(e)}")
    
    def _merge_metrics(self, saved_metrics: Dict[str, Any]) -> None:
        """
        Merge saved metrics with current metrics.
        
        Args:
            saved_metrics: Dictionary of saved metrics
        """
        with self.lock:
            # Preserve history
            if "history" in saved_metrics:
                self.metrics["history"] = saved_metrics["history"]
            
            # Preserve operation counts
            if "operations" in saved_metrics:
                for key in self.metrics["operations"]:
                    if key in saved_metrics["operations"]:
                        self.metrics["operations"][key] = saved_metrics["operations"][key]
            
            # Preserve throughput
            if "throughput" in saved_metrics:
                for key in self.metrics["throughput"]:
                    if key in saved_metrics["throughput"]:
                        self.metrics["throughput"][key] = saved_metrics["throughput"][key]
            
            # Preserve database metrics
            if "database" in saved_metrics:
                for key in self.metrics["database"]:
                    if key in saved_metrics["database"]:
                        self.metrics["database"][key] = saved_metrics["database"][key]
            
            # Preserve custom metrics
            if "custom_metrics" in saved_metrics:
                self.metrics["custom_metrics"] = saved_metrics["custom_metrics"]
    
    def record_operation(self, operation_type: str, success: bool, execution_time: float,
                        items_processed: int = 0, bytes_processed: int = 0) -> None:
        """
        Record a database operation.
        
        Args:
            operation_type: Type of operation (e.g., "query", "transaction")
            success: Whether the operation was successful
            execution_time: Time taken to execute the operation (in seconds)
            items_processed: Number of items processed in the operation
            bytes_processed: Number of bytes processed in the operation
        """
        with self.lock:
            # Update operation counts
            self.metrics["operations"]["total"] += 1
            if success:
                self.metrics["operations"]["successful"] += 1
            else:
                self.metrics["operations"]["failed"] += 1
            
            # Update timing metrics
            self.metrics["timing"]["total_execution_time"] += execution_time
            self.metrics["timing"]["min_execution_time"] = min(
                self.metrics["timing"]["min_execution_time"], execution_time
            ) if self.metrics["timing"]["min_execution_time"] != float('inf') else execution_time
            self.metrics["timing"]["max_execution_time"] = max(
                self.metrics["timing"]["max_execution_time"], execution_time
            )
            self.metrics["timing"]["avg_execution_time"] = (
                self.metrics["timing"]["total_execution_time"] / 
                self.metrics["operations"]["total"]
            )
            
            # Update throughput metrics
            self.metrics["throughput"]["items_processed"] += items_processed
            self.metrics["throughput"]["bytes_processed"] += bytes_processed
            
            elapsed = time.time() - self.start_time
            if elapsed > 0:
                self.metrics["throughput"]["items_per_second"] = (
                    self.metrics["throughput"]["items_processed"] / elapsed
                )
                self.metrics["throughput"]["bytes_per_second"] = (
                    self.metrics["throughput"]["bytes_processed"] / elapsed
                )
            
            # Update memory metrics
            process = psutil.Process(os.getpid())
            current_memory = process.memory_info().rss
            self.metrics["memory"]["current_memory_usage"] = current_memory
            self.metrics["memory"]["peak_memory_usage"] = max(
                self.metrics["memory"]["peak_memory_usage"], current_memory
            )
            
            # Update database metrics based on operation type
            if operation_type == "connection":
                self.metrics["database"]["connection_attempts"] += 1
                if not success:
                    self.metrics["database"]["connection_failures"] += 1
            elif operation_type == "query":
                self.metrics["database"]["query_count"] += 1
            elif operation_type == "transaction":
                self.metrics["database"]["transaction_count"] += 1
            elif operation_type == "rollback":
                self.metrics["database"]["rollback_count"] += 1
            
            # Add to history (with limited size to prevent memory issues)
            history_entry = {
                "timestamp": time.time(),
                "operation_type": operation_type,
                "success": success,
                "execution_time": execution_time,
                "items_processed": items_processed
            }
            
            self.metrics["history"].append(history_entry)
            if len(self.metrics["history"]) > 1000:  # Keep last 1000 operations
                self.metrics["history"] = self.metrics["history"][-1000:]
            
            # Save metrics to file if specified
            if self.metrics_file:
                try:
                    with open(self.metrics_file, 'w') as f:
                        json.dump(self.metrics, f, indent=2)
                except Exception as e:
                    logger.warning(f"Could not save metrics to {self.metrics_file}: {str(e)}")
    
    def record_custom_metric(self, metric_name: str, value: Any) -> None:
        """
        Record a custom metric.
        
        Args:
            metric_name: Name of the metric
            value: Value of the metric
        """
        with self.lock:
            self.metrics["custom_metrics"][metric_name] = value
            
            # Save metrics to file if specified
            if self.metrics_file:
                try:
                    with open(self.metrics_file, 'w') as f:
                        json.dump(self.metrics, f, indent=2)
                except Exception as e:
                    logger.warning(f"Could not save metrics to {self.metrics_file}: {str(e)}")
    
    def get_metrics(self) -> Dict[str, Any]:
        """
        Get current metrics.
        
        Returns:
            Dictionary of current metrics
        """
        with self.lock:
            # Update elapsed time and throughput calculations
            elapsed = time.time() - self.start_time
            self.metrics["elapsed_time"] = elapsed
            
            if elapsed > 0:
                self.metrics["throughput"]["items_per_second"] = (
                    self.metrics["throughput"]["items_processed"] / elapsed
                )
                self.metrics["throughput"]["bytes_per_second"] = (
                    self.metrics["throughput"]["bytes_processed"] / elapsed
                )
            
            # Update memory metrics
            process = psutil.Process(os.getpid())
            self.metrics["memory"]["current_memory_usage"] = process.memory_info().rss
            
            return self.metrics.copy()
    
    def generate_report(self) -> str:
        """
        Generate a human-readable report of the metrics.
        
        Returns:
            String containing the report
        """
        metrics = self.get_metrics()
        elapsed = metrics["elapsed_time"]
        
        report = [
            f"Performance Report for {self.name}",
            f"=======================================",
            f"Duration: {timedelta(seconds=int(elapsed))}",
            f"",
            f"Operations:",
            f"  Total: {metrics['operations']['total']}",
            f"  Successful: {metrics['operations']['successful']}",
            f"  Failed: {metrics['operations']['failed']}",
            f"  Success Rate: {metrics['operations']['successful'] / metrics['operations']['total'] * 100:.2f}% (if total > 0)",
            f"",
            f"Timing:",
            f"  Average Execution Time: {metrics['timing']['avg_execution_time'] * 1000:.2f} ms",
            f"  Min Execution Time: {metrics['timing']['min_execution_time'] * 1000:.2f} ms (if not inf)",
            f"  Max Execution Time: {metrics['timing']['max_execution_time'] * 1000:.2f} ms",
            f"",
            f"Throughput:",
            f"  Items Processed: {metrics['throughput']['items_processed']}",
            f"  Items Per Second: {metrics['throughput']['items_per_second']:.2f}",
            f"  Bytes Processed: {metrics['throughput']['bytes_processed']}",
            f"  Bytes Per Second: {metrics['throughput']['bytes_per_second']:.2f}",
            f"",
            f"Memory:",
            f"  Current Memory Usage: {metrics['memory']['current_memory_usage'] / (1024 * 1024):.2f} MB",
            f"  Peak Memory Usage: {metrics['memory']['peak_memory_usage'] / (1024 * 1024):.2f} MB",
            f"",
            f"Database:",
            f"  Connection Attempts: {metrics['database']['connection_attempts']}",
            f"  Connection Failures: {metrics['database']['connection_failures']}",
            f"  Query Count: {metrics['database']['query_count']}",
            f"  Transaction Count: {metrics['database']['transaction_count']}",
            f"  Rollback Count: {metrics['database']['rollback_count']}",
            f"",
            f"Custom Metrics:"
        ]
        
        for name, value in metrics["custom_metrics"].items():
            report.append(f"  {name}: {value}")
        
        return "\n".join(report)
    
    def reset(self) -> None:
        """Reset metrics to initial state."""
        with self.lock:
            self.start_time = time.time()
            self.metrics = {
                "name": self.name,
                "start_time": self.start_time,
                "operations": {
                    "total": 0,
                    "successful": 0,
                    "failed": 0
                },
                "timing": {
                    "total_execution_time": 0,
                    "min_execution_time": float('inf'),
                    "max_execution_time": 0,
                    "avg_execution_time": 0
                },
                "throughput": {
                    "items_processed": 0,
                    "items_per_second": 0,
                    "bytes_processed": 0,
                    "bytes_per_second": 0
                },
                "memory": {
                    "peak_memory_usage": 0,
                    "current_memory_usage": 0
                },
                "database": {
                    "connection_attempts": 0,
                    "connection_failures": 0,
                    "query_count": 0,
                    "transaction_count": 0,
                    "rollback_count": 0
                },
                "custom_metrics": {},
                "history": []
            }
            
            # Save reset metrics to file if specified
            if self.metrics_file:
                try:
                    with open(self.metrics_file, 'w') as f:
                        json.dump(self.metrics, f, indent=2)
                except Exception as e:
                    logger.warning(f"Could not save reset metrics to {self.metrics_file}: {str(e)}")
class ConnectionHealthMonitor:
    """
    Monitors database connection health and provides automatic fallback.
    
    This class integrates with ConnectionManager to monitor connection health
    and provide automatic fallback to alternative connection methods when
    the primary connection fails.
    """
    
    _instance = None
    _lock = threading.Lock()
    
    @classmethod
    def get_instance(cls) -> 'ConnectionHealthMonitor':
        """
        Get singleton instance of ConnectionHealthMonitor.
        
        Returns:
            ConnectionHealthMonitor: Singleton instance
        """
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = ConnectionHealthMonitor()
        return cls._instance
    
    def __init__(self):
        """Initialize connection health monitor."""
        if ConnectionHealthMonitor._instance is not None:
            raise RuntimeError("ConnectionHealthMonitor is a singleton. Use get_instance() instead.")
            
        ConnectionHealthMonitor._instance = self
        
        # Initialize metrics
        self.metrics = PerformanceMetrics("connection_health", "monitoring/connection_health.json")
        
        # Initialize connection manager if available
        self.connection_manager = ConnectionManager.get_instance() if ConnectionManager else None
        
        # Initialize health check interval
        self.health_check_interval = int(os.getenv('CONNECTION_HEALTH_CHECK_INTERVAL', '60'))
        
        # Initialize health check thread
        self.health_check_thread = None
        self.stop_health_check = threading.Event()
        
        # Initialize connection health status
        self.health_status = {
            "last_check_time": 0,
            "is_healthy": False,
            "active_connection_type": None,
            "connection_types": {},
            "recent_failures": [],
            "recent_successes": []
        }
        
        logger.info("ConnectionHealthMonitor initialized")
    
    def start_monitoring(self) -> None:
        """Start connection health monitoring thread."""
        if self.health_check_thread is not None and self.health_check_thread.is_alive():
            logger.warning("Health check thread is already running")
            return
            
        self.stop_health_check.clear()
        self.health_check_thread = threading.Thread(
            target=self._health_check_loop,
            daemon=True
        )
        self.health_check_thread.start()
        logger.info(f"Started connection health monitoring (interval: {self.health_check_interval}s)")
    
    def stop_monitoring(self) -> None:
        """Stop connection health monitoring thread."""
        if self.health_check_thread is None or not self.health_check_thread.is_alive():
            logger.warning("Health check thread is not running")
            return
            
        self.stop_health_check.set()
        self.health_check_thread.join(timeout=5)
        if self.health_check_thread.is_alive():
            logger.warning("Health check thread did not terminate gracefully")
        else:
            logger.info("Stopped connection health monitoring")
    
    def _health_check_loop(self) -> None:
        """Run health check loop in a separate thread."""
        while not self.stop_health_check.is_set():
            try:
                self.check_connection_health()
            except Exception as e:
                logger.error(f"Error in health check loop: {str(e)}")
            
            # Sleep until next check interval or stop event is set
            self.stop_health_check.wait(self.health_check_interval)
    
    def check_connection_health(self) -> Dict[str, Any]:
        """
        Check connection health and update status.
        
        Returns:
            Dictionary with connection health status
        """
        if not self.connection_manager:
            logger.warning("ConnectionManager not available, cannot check connection health")
            self.health_status["is_healthy"] = False
            self.health_status["last_check_time"] = time.time()
            return self.health_status
        
        start_time = time.time()
        
        # Get connection manager stats
        conn_stats = self.connection_manager.get_stats() if hasattr(self.connection_manager, 'get_stats') else {}
        
        # Check circuit breaker states
        circuit_breakers = {}
        if hasattr(self.connection_manager, 'circuit_breakers'):
            for name, breaker in self.connection_manager.circuit_breakers.items():
                circuit_breakers[name] = breaker.get_state()
        
        # Check active connection pool
        active_pool = self.connection_manager.active_pool if hasattr(self.connection_manager, 'active_pool') else None
        
        # Try to get a test connection
        connection = None
        is_healthy = False
        error_message = None
        
        try:
            connection = self.connection_manager.get_connection()
            if connection:
                # Test the connection with a simple query
                if hasattr(connection, 'execute_query'):
                    result = connection.execute_query("SELECT 1 as test")
                    is_healthy = result and len(result) > 0 and result[0].get('test') == 1
                else:
                    # Assume connection is healthy if we got one
                    is_healthy = True
                
                self.metrics.record_operation(
                    operation_type="connection",
                    success=is_healthy,
                    execution_time=time.time() - start_time
                )
                
                if is_healthy:
                    self.health_status["recent_successes"].append({
                        "timestamp": time.time(),
                        "connection_type": active_pool,
                        "latency": time.time() - start_time
                    })
                    # Keep only last 10 successes
                    if len(self.health_status["recent_successes"]) > 10:
                        self.health_status["recent_successes"] = self.health_status["recent_successes"][-10:]
                else:
                    error_message = "Connection test query failed"
                    self.health_status["recent_failures"].append({
                        "timestamp": time.time(),
                        "connection_type": active_pool,
                        "error": error_message
                    })
                    # Keep only last 10 failures
                    if len(self.health_status["recent_failures"]) > 10:
                        self.health_status["recent_failures"] = self.health_status["recent_failures"][-10:]
        except Exception as e:
            is_healthy = False
            error_message = str(e)
            
            self.metrics.record_operation(
                operation_type="connection",
                success=False,
                execution_time=time.time() - start_time
            )
            
            self.health_status["recent_failures"].append({
                "timestamp": time.time(),
                "connection_type": active_pool,
                "error": error_message
            })
            # Keep only last 10 failures
            if len(self.health_status["recent_failures"]) > 10:
                self.health_status["recent_failures"] = self.health_status["recent_failures"][-10:]
        finally:
            if connection and hasattr(connection, 'close'):
                connection.close()
        
        # Update health status
        self.health_status.update({
            "last_check_time": time.time(),
            "is_healthy": is_healthy,
            "active_connection_type": active_pool,
            "connection_types": circuit_breakers,
            "check_duration": time.time() - start_time,
            "error_message": error_message
        })
        
        # Log health status
        if is_healthy:
            logger.info(f"Connection health check passed (type: {active_pool}, "
                       f"latency: {self.health_status['check_duration']:.3f}s)")
        else:
            logger.warning(f"Connection health check failed (type: {active_pool}, "
                          f"error: {error_message})")
        
        return self.health_status
    
    def get_health_status(self) -> Dict[str, Any]:
        """
        Get current connection health status.
        
        Returns:
            Dictionary with connection health status
        """
        # Check if we need to refresh the status
        if time.time() - self.health_status["last_check_time"] > self.health_check_interval:
            return self.check_connection_health()
        
        return self.health_status
class ProgressTracker:
    """
    Tracks progress of batch operations with time estimation.
    
    This class provides methods to track progress of batch operations,
    estimate completion time, and generate progress reports.
    """
    
    def __init__(self, name: str, total_items: int, checkpoint_file: Optional[str] = None):
        """
        Initialize progress tracker.
        
        Args:
            name: Name of the operation being tracked
            total_items: Total number of items to process
            checkpoint_file: Optional path to checkpoint file for resumable operations
        """
        self.name = name
        self.total_items = total_items
        self.checkpoint_file = checkpoint_file
        self.start_time = time.time()
        self.last_update_time = self.start_time
        self.processed_items = 0
        self.successful_items = 0
        self.failed_items = 0
        self.current_batch = 0
        self.total_batches = 0
        self.batch_sizes = []
        self.batch_times = []
        self.lock = threading.RLock()
        
        # Generate a unique ID for this tracking session
        self.tracking_id = str(uuid.uuid4())
        
        # Load checkpoint if exists
        if checkpoint_file and os.path.exists(checkpoint_file):
            try:
                with open(checkpoint_file, 'r') as f:
                    checkpoint = json.load(f)
                    if "processed" in checkpoint:
                        self.processed_items = checkpoint["processed"]
                        logger.info(f"Loaded progress from checkpoint: {self.processed_items}/{self.total_items} items processed")
            except Exception as e:
                logger.warning(f"Could not load checkpoint from {checkpoint_file}: {str(e)}")
        
        # Initialize metrics
        self.metrics = PerformanceMetrics(f"progress_{name}", f"monitoring/progress_{name}.json")
        
        logger.info(f"Progress tracker initialized for {name}: {self.processed_items}/{total_items} items")
    
    def update(self, items_processed: int, successful: int = None, failed: int = None,
              batch_size: int = None, save_checkpoint: bool = True) -> Dict[str, Any]:
        """
        Update progress.
        
        Args:
            items_processed: Number of items processed in this update
            successful: Number of items successfully processed (default: all)
            failed: Number of items that failed processing (default: 0)
            batch_size: Size of the current batch (for batch operations)
            save_checkpoint: Whether to save checkpoint (if checkpoint_file is set)
            
        Returns:
            Dictionary with current progress status
        """
        with self.lock:
            current_time = time.time()
            
            # Set defaults
            if successful is None:
                successful = items_processed
            if failed is None:
                failed = 0
            
            # Update counters
            self.processed_items += items_processed
            self.successful_items += successful
            self.failed_items += failed
            
            # Ensure we don't exceed total items
            self.processed_items = min(self.processed_items, self.total_items)
            
            # Update batch information
            if batch_size is not None:
                self.current_batch += 1
                self.batch_sizes.append(batch_size)
                self.batch_times.append(current_time - self.last_update_time)
            
            # Calculate progress percentage
            progress_pct = (self.processed_items / self.total_items) * 100 if self.total_items > 0 else 0
            
            # Calculate elapsed time
            elapsed = current_time - self.start_time
            
            # Calculate items per second
            items_per_sec = self.processed_items / elapsed if elapsed > 0 else 0
            
            # Calculate estimated time remaining
            remaining_items = self.total_items - self.processed_items
            est_remaining_time = remaining_items / items_per_sec if items_per_sec > 0 else 0
            
            # Calculate estimated completion time
            est_completion = datetime.now() + timedelta(seconds=est_remaining_time)
            
            # Calculate average batch processing time (from last 10 batches)
            recent_batch_times = self.batch_times[-10:] if self.batch_times else [0]
            avg_batch_time = sum(recent_batch_times) / len(recent_batch_times)
            
            # Calculate average batch size (from last 10 batches)
            recent_batch_sizes = self.batch_sizes[-10:] if self.batch_sizes else [0]
            avg_batch_size = sum(recent_batch_sizes) / len(recent_batch_sizes) if recent_batch_sizes else 0
            
            # Calculate items per batch
            items_per_batch = avg_batch_size
            
            # Calculate batches per second
            batches_per_sec = 1 / avg_batch_time if avg_batch_time > 0 else 0
            
            # Update metrics
            self.metrics.record_operation(
                operation_type="batch",
                success=(failed == 0),
                execution_time=current_time - self.last_update_time,
                items_processed=items_processed
            )
            
            # Save checkpoint if requested
            if save_checkpoint and self.checkpoint_file:
                checkpoint = {
                    "position": self.processed_items,
                    "processed": self.processed_items,
                    "successful": self.successful_items,
                    "failed": self.failed_items,
                    "last_updated": datetime.now().isoformat(),
                    "last_batch": self.current_batch
                }
                
                try:
                    with open(self.checkpoint_file, 'w') as f:
                        json.dump(checkpoint, f)
                except Exception as e:
                    logger.warning(f"Could not save checkpoint to {self.checkpoint_file}: {str(e)}")
            
            # Update last update time
            self.last_update_time = current_time
            
            # Create progress status
            status = {
                "name": self.name,
                "tracking_id": self.tracking_id,
                "total_items": self.total_items,
                "processed_items": self.processed_items,
                "successful_items": self.successful_items,
                "failed_items": self.failed_items,
                "progress_percentage": progress_pct,
                "elapsed_time": elapsed,
                "items_per_second": items_per_sec,
                "estimated_remaining_time": est_remaining_time,
                "estimated_completion": est_completion.isoformat(),
                "current_batch": self.current_batch,
                "avg_batch_time": avg_batch_time,
                "avg_batch_size": avg_batch_size,
                "items_per_batch": items_per_batch,
                "batches_per_second": batches_per_sec,
                "last_update_time": current_time
            }
            
            # Log progress
            logger.info(f"Progress [{self.name}]: {self.processed_items}/{self.total_items} items "
                       f"({progress_pct:.1f}%) at {items_per_sec:.1f} items/sec, "
                       f"ETA: {est_completion.strftime('%Y-%m-%d %H:%M:%S')}")
            
            return status
    
    def get_status(self) -> Dict[str, Any]:
        """
        Get current progress status.
        
        Returns:
            Dictionary with current progress status
        """
        with self.lock:
            current_time = time.time()
            
            # Calculate progress percentage
            progress_pct = (self.processed_items / self.total_items) * 100 if self.total_items > 0 else 0
            
            # Calculate elapsed time
            elapsed = current_time - self.start_time
            
            # Calculate items per second
            items_per_sec = self.processed_items / elapsed if elapsed > 0 else 0
            
            # Calculate estimated time remaining
            remaining_items = self.total_items - self.processed_items
            est_remaining_time = remaining_items / items_per_sec if items_per_sec > 0 else 0
            
            # Calculate estimated completion time
            est_completion = datetime.now() + timedelta(seconds=est_remaining_time)
            
            # Calculate average batch processing time (from last 10 batches)
            recent_batch_times = self.batch_times[-10:] if self.batch_times else [0]
            avg_batch_time = sum(recent_batch_times) / len(recent_batch_times)
            
            # Calculate average batch size (from last 10 batches)
            recent_batch_sizes = self.batch_sizes[-10:] if self.batch_sizes else [0]
            avg_batch_size = sum(recent_batch_sizes) / len(recent_batch_sizes) if recent_batch_sizes else 0
            
            # Calculate items per batch
            items_per_batch = avg_batch_size
            
            # Calculate batches per second
            batches_per_sec = 1 / avg_batch_time if avg_batch_time > 0 else 0
            
            # Create progress status
            return {
                "name": self.name,
                "tracking_id": self.tracking_id,
                "total_items": self.total_items,
                "processed_items": self.processed_items,
                "successful_items": self.successful_items,
                "failed_items": self.failed_items,
                "progress_percentage": progress_pct,
                "elapsed_time": elapsed,
                "items_per_second": items_per_sec,
                "estimated_remaining_time": est_remaining_time,
                "estimated_completion": est_completion.isoformat(),
                "current_batch": self.current_batch,
                "avg_batch_time": avg_batch_time,
                "avg_batch_size": avg_batch_size,
                "items_per_batch": items_per_batch,
                "batches_per_second": batches_per_sec
            }
    
    def generate_report(self) -> str:
        """
        Generate a human-readable progress report.
        
        Returns:
            String containing the progress report
        """
        status = self.get_status()
        
        report = [
            f"Progress Report for {self.name}",
            f"=======================================",
            f"Status: {status['processed_items']}/{status['total_items']} items processed "
            f"({status['progress_percentage']:.1f}%)",
            f"",
            f"Time:",
            f"  Elapsed: {timedelta(seconds=int(status['elapsed_time']))}",
            f"  Estimated Remaining: {timedelta(seconds=int(status['estimated_remaining_time']))}",
            f"  Estimated Completion: {datetime.fromisoformat(status['estimated_completion']).strftime('%Y-%m-%d %H:%M:%S')}",
            f"",
            f"Performance:",
            f"  Items Per Second: {status['items_per_second']:.1f}",
            f"  Batches Per Second: {status['batches_per_second']:.2f}",
            f"  Average Batch Size: {status['avg_batch_size']:.1f} items",
            f"  Average Batch Time: {status['avg_batch_time'] * 1000:.1f} ms",
            f"",
            f"Results:",
            f"  Successful: {status['successful_items']} items",
            f"  Failed: {status['failed_items']} items",
            f"  Success Rate: {status['successful_items'] / status['processed_items'] * 100:.1f}% (if processed > 0)"
        ]
        
        return "\n".join(report)
    
    def reset(self) -> None:
        """Reset progress tracker to initial state."""
        with self.lock:
            self.start_time = time.time()
            self.last_update_time = self.start_time
            self.processed_items = 0
            self.successful_items = 0
            self.failed_items = 0
            self.current_batch = 0
            self.batch_sizes = []
            self.batch_times = []
            
            # Generate a new tracking ID
            self.tracking_id = str(uuid.uuid4())
            
            # Reset metrics
            self.metrics.reset()
            
            # Save reset checkpoint if specified
            if self.checkpoint_file:
                checkpoint = {
                    "position": 0,
                    "processed": 0,
                    "successful": 0,
                    "failed": 0,
                    "last_updated": datetime.now().isoformat(),
                    "last_batch": 0
                }
                
                try:
                    with open(self.checkpoint_file, 'w') as f:
                        json.dump(checkpoint, f)
                except Exception as e:
                    logger.warning(f"Could not save reset checkpoint to {self.checkpoint_file}: {str(e)}")
            
            logger.info(f"Progress tracker reset for {self.name}")
# Utility functions for monitoring and metrics collection

def create_monitoring_directory() -> None:
    """
    Create monitoring directory if it doesn't exist.
    
    This function creates the monitoring directory for storing metrics,
    checkpoints, and logs.
    """
    monitoring_dir = "monitoring"
    if not os.path.exists(monitoring_dir):
        try:
            os.makedirs(monitoring_dir)
            logger.info(f"Created monitoring directory: {monitoring_dir}")
        except Exception as e:
            logger.warning(f"Could not create monitoring directory: {str(e)}")


def monitor_connection_health(interval: int = 60) -> None:
    """
    Start connection health monitoring.
    
    This function starts the connection health monitoring thread that
    periodically checks the health of database connections.
    
    Args:
        interval: Health check interval in seconds
    """
    # Create monitoring directory
    create_monitoring_directory()
    
    # Set health check interval
    os.environ['CONNECTION_HEALTH_CHECK_INTERVAL'] = str(interval)
    
    # Start monitoring
    monitor = ConnectionHealthMonitor.get_instance()
    monitor.start_monitoring()
    
    logger.info(f"Started connection health monitoring with interval {interval}s")


def track_progress(name: str, total_items: int, checkpoint_file: Optional[str] = None) -> ProgressTracker:
    """
    Create a progress tracker for batch operations.
    
    This function creates a progress tracker for monitoring batch operations
    and estimating completion time.
    
    Args:
        name: Name of the operation being tracked
        total_items: Total number of items to process
        checkpoint_file: Optional path to checkpoint file for resumable operations
        
    Returns:
        ProgressTracker instance
    """
    # Create monitoring directory
    create_monitoring_directory()
    
    # Create checkpoint directory if needed
    if checkpoint_file:
        checkpoint_dir = os.path.dirname(checkpoint_file)
        if checkpoint_dir and not os.path.exists(checkpoint_dir):
            try:
                os.makedirs(checkpoint_dir)
                logger.info(f"Created checkpoint directory: {checkpoint_dir}")
            except Exception as e:
                logger.warning(f"Could not create checkpoint directory: {str(e)}")
    
    # Create progress tracker
    tracker = ProgressTracker(name, total_items, checkpoint_file)
    
    return tracker


def collect_performance_metrics(name: str, metrics_file: Optional[str] = None) -> PerformanceMetrics:
    """
    Create a performance metrics collector.
    
    This function creates a performance metrics collector for tracking
    database operation performance.
    
    Args:
        name: Name of the metrics collector
        metrics_file: Optional path to file for persisting metrics
        
    Returns:
        PerformanceMetrics instance
    """
    # Create monitoring directory
    create_monitoring_directory()
    
    # Create metrics directory if needed
    if metrics_file:
        metrics_dir = os.path.dirname(metrics_file)
        if metrics_dir and not os.path.exists(metrics_dir):
            try:
                os.makedirs(metrics_dir)
                logger.info(f"Created metrics directory: {metrics_dir}")
            except Exception as e:
                logger.warning(f"Could not create metrics directory: {str(e)}")
    
    # Create performance metrics collector
    metrics = PerformanceMetrics(name, metrics_file)
    
    return metrics


def get_connection_health_status() -> Dict[str, Any]:
    """
    Get current connection health status.
    
    This function returns the current connection health status from the
    ConnectionHealthMonitor.
    
    Returns:
        Dictionary with connection health status
    """
    monitor = ConnectionHealthMonitor.get_instance()
    return monitor.get_health_status()


def generate_monitoring_report() -> str:
    """
    Generate a comprehensive monitoring report.
    
    This function generates a comprehensive report of all monitoring metrics,
    including connection health, progress tracking, and performance metrics.
    
    Returns:
        String containing the monitoring report
    """
    # Create monitoring directory
    create_monitoring_directory()
    
    # Get connection health status
    health_status = get_connection_health_status()
    
    # Build report
    report = [
        "CryoProtect v2 Database Monitoring Report",
        "==========================================",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "Connection Health:",
        f"  Status: {'Healthy' if health_status['is_healthy'] else 'Unhealthy'}",
        f"  Active Connection Type: {health_status['active_connection_type']}",
        f"  Last Check: {datetime.fromtimestamp(health_status['last_check_time']).strftime('%Y-%m-%d %H:%M:%S')}",
        f"  Check Duration: {health_status.get('check_duration', 0) * 1000:.2f} ms",
    ]
    
    if not health_status['is_healthy'] and 'error_message' in health_status:
        report.append(f"  Error: {health_status['error_message']}")
    
    report.append("")
    report.append("Circuit Breaker Status:")
    
    if 'connection_types' in health_status:
        for conn_type, state in health_status['connection_types'].items():
            report.append(f"  {conn_type}: {state['state']}")
    else:
        report.append("  No circuit breaker information available")
    
    report.append("")
    report.append("Recent Connection Failures:")
    
    if 'recent_failures' in health_status and health_status['recent_failures']:
        for failure in health_status['recent_failures'][-5:]:  # Show last 5 failures
            timestamp = datetime.fromtimestamp(failure['timestamp']).strftime('%Y-%m-%d %H:%M:%S')
            report.append(f"  {timestamp} - {failure['connection_type']}: {failure['error']}")
    else:
        report.append("  No recent connection failures")
    
    # Add metrics files information
    report.append("")
    report.append("Available Metrics Files:")
    
    metrics_files = []
    if os.path.exists("monitoring"):
        metrics_files = [f for f in os.listdir("monitoring") if f.endswith(".json")]
    
    if metrics_files:
        for metrics_file in metrics_files:
            file_path = os.path.join("monitoring", metrics_file)
            file_size = os.path.getsize(file_path)
            file_modified = datetime.fromtimestamp(os.path.getmtime(file_path)).strftime('%Y-%m-%d %H:%M:%S')
            report.append(f"  {metrics_file} ({file_size/1024:.1f} KB, last modified: {file_modified})")
    else:
        report.append("  No metrics files available")
    
    return "\n".join(report)


# Context manager for timing and recording database operations
@contextmanager
def timed_operation(metrics: PerformanceMetrics, operation_type: str, items_count: int = 0):
    """
    Context manager for timing and recording database operations.
    
    This context manager measures the execution time of a database operation
    and records it in the provided metrics collector.
    
    Args:
        metrics: PerformanceMetrics instance
        operation_type: Type of operation (e.g., "query", "transaction")
        items_count: Number of items processed in the operation
        
    Yields:
        None
    """
    start_time = time.time()
    success = True
    
    try:
        yield
    except Exception as e:
        success = False
        raise
    finally:
        execution_time = time.time() - start_time
        metrics.record_operation(
            operation_type=operation_type,
            success=success,
            execution_time=execution_time,
            items_processed=items_count
        )