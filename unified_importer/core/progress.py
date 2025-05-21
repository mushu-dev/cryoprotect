"""
Progress tracking and reporting for molecular data imports.

This module provides utilities for tracking import progress, calculating 
estimated completion times, and generating progress reports.
"""

import time
import json
import logging
import datetime
from typing import Dict, Any, Optional, List, Callable
from collections import deque


class ProgressTracker:
    """
    Track progress of molecular data import operations.
    
    This class tracks metrics like completed items, failed items, 
    processing rates, and estimated completion times.
    """
    
    def __init__(
        self,
        total_items: int = 0,
        window_size: int = 100,
        logger: Optional[logging.Logger] = None
    ):
        """
        Initialize progress tracker.
        
        Args:
            total_items: Total number of items to process
            window_size: Size of the window for rate calculations
            logger: Logger instance for progress updates
        """
        self.total_items = total_items
        self.processed_items = 0
        self.successful_items = 0
        self.failed_items = 0
        self.skipped_items = 0
        
        self.start_time = time.time()
        self.last_update_time = self.start_time
        
        # For rate calculations
        self.window_size = window_size
        self.processing_times = deque(maxlen=window_size)
        self.batch_sizes = deque(maxlen=window_size)
        
        # Error tracking
        self.error_counts: Dict[str, int] = {}
        
        # Logger
        self.logger = logger or logging.getLogger(__name__)
        
        # For progress callbacks
        self.progress_callbacks: List[Callable[[Dict[str, Any]], None]] = []
    
    def update(
        self,
        processed: int = 0,
        successful: int = 0,
        failed: int = 0,
        skipped: int = 0,
        errors: Optional[Dict[str, int]] = None
    ) -> None:
        """
        Update progress metrics.
        
        Args:
            processed: Number of processed items in this batch
            successful: Number of successfully processed items in this batch
            failed: Number of failed items in this batch
            skipped: Number of skipped items in this batch
            errors: Dictionary of error counts by type
        """
        current_time = time.time()
        elapsed = current_time - self.last_update_time
        
        # Update counters
        self.processed_items += processed
        self.successful_items += successful
        self.failed_items += failed
        self.skipped_items += skipped
        
        # Update error counts
        if errors:
            for error_type, count in errors.items():
                self.error_counts[error_type] = self.error_counts.get(error_type, 0) + count
        
        # Update rate calculations
        if processed > 0 and elapsed > 0:
            self.processing_times.append(elapsed)
            self.batch_sizes.append(processed)
        
        self.last_update_time = current_time
        
        # Trigger callbacks
        progress_data = self.get_progress_data()
        for callback in self.progress_callbacks:
            try:
                callback(progress_data)
            except Exception as e:
                self.logger.error(f"Error in progress callback: {str(e)}")
    
    def register_callback(self, callback: Callable[[Dict[str, Any]], None]) -> None:
        """
        Register a callback function for progress updates.
        
        Args:
            callback: Function that takes a progress data dictionary
        """
        if callback not in self.progress_callbacks:
            self.progress_callbacks.append(callback)
    
    def unregister_callback(self, callback: Callable[[Dict[str, Any]], None]) -> None:
        """
        Unregister a progress callback function.
        
        Args:
            callback: Previously registered callback function
        """
        if callback in self.progress_callbacks:
            self.progress_callbacks.remove(callback)
    
    def get_processing_rate(self) -> float:
        """
        Calculate the current processing rate (items per second).
        
        Returns:
            Processing rate or 0 if not enough data
        """
        if not self.processing_times or not self.batch_sizes:
            if self.processed_items > 0 and time.time() > self.start_time:
                # Fall back to overall average if no window data
                total_elapsed = time.time() - self.start_time
                return self.processed_items / total_elapsed
            return 0
            
        # Calculate based on recent window
        total_time = sum(self.processing_times)
        total_items = sum(self.batch_sizes)
        
        if total_time > 0:
            return total_items / total_time
        return 0
    
    def get_estimated_completion_time(self) -> Optional[float]:
        """
        Estimate the completion time based on current progress.
        
        Returns:
            Estimated completion time as unix timestamp, or None if can't estimate
        """
        if (self.total_items <= 0 or 
            self.processed_items <= 0 or 
            self.processed_items >= self.total_items):
            return None
            
        rate = self.get_processing_rate()
        if rate <= 0:
            return None
            
        remaining_items = self.total_items - self.processed_items
        estimated_seconds = remaining_items / rate
        
        return time.time() + estimated_seconds
    
    def get_elapsed_time(self) -> float:
        """
        Get the elapsed time since the start of processing.
        
        Returns:
            Elapsed time in seconds
        """
        return time.time() - self.start_time
    
    def get_progress_percentage(self) -> float:
        """
        Calculate the current progress percentage.
        
        Returns:
            Progress as a percentage (0-100)
        """
        if self.total_items <= 0:
            return 0
        return min(100.0, (self.processed_items / self.total_items) * 100)
    
    def get_progress_data(self) -> Dict[str, Any]:
        """
        Get comprehensive progress data.
        
        Returns:
            Dictionary with all progress metrics
        """
        current_time = time.time()
        elapsed = current_time - self.start_time
        
        estimated_completion = self.get_estimated_completion_time()
        eta_string = None
        
        if estimated_completion:
            eta_seconds = estimated_completion - current_time
            eta_string = str(datetime.timedelta(seconds=int(eta_seconds)))
        
        return {
            'timestamp': current_time,
            'elapsed_seconds': elapsed,
            'elapsed_formatted': str(datetime.timedelta(seconds=int(elapsed))),
            'total_items': self.total_items,
            'processed_items': self.processed_items,
            'successful_items': self.successful_items,
            'failed_items': self.failed_items,
            'skipped_items': self.skipped_items,
            'progress_percentage': self.get_progress_percentage(),
            'processing_rate': self.get_processing_rate(),
            'items_per_minute': self.get_processing_rate() * 60,
            'estimated_completion': estimated_completion,
            'estimated_time_remaining': eta_string,
            'error_counts': self.error_counts.copy()
        }
    
    def log_progress(self, level: int = logging.INFO) -> None:
        """
        Log current progress information.
        
        Args:
            level: Logging level to use
        """
        progress_data = self.get_progress_data()
        
        msg = (
            f"Progress: {progress_data['progress_percentage']:.1f}% "
            f"({progress_data['processed_items']}/{progress_data['total_items']})"
        )
        
        if progress_data['estimated_time_remaining']:
            msg += f", ETA: {progress_data['estimated_time_remaining']}"
            
        msg += f", Rate: {progress_data['items_per_minute']:.1f} items/min"
        
        if progress_data['failed_items'] > 0:
            msg += f", Failed: {progress_data['failed_items']}"
            
        if progress_data['skipped_items'] > 0:
            msg += f", Skipped: {progress_data['skipped_items']}"
            
        self.logger.log(level, msg, extra={'structured_data': progress_data})
    
    def generate_report(self) -> Dict[str, Any]:
        """
        Generate a detailed progress report.
        
        Returns:
            Dictionary with detailed report information
        """
        progress_data = self.get_progress_data()
        
        # Add additional report metrics
        progress_data.update({
            'success_rate': (
                (self.successful_items / self.processed_items) * 100
                if self.processed_items > 0 else 0
            ),
            'error_rate': (
                (self.failed_items / self.processed_items) * 100
                if self.processed_items > 0 else 0
            ),
            'skipped_rate': (
                (self.skipped_items / self.processed_items) * 100
                if self.processed_items > 0 else 0
            ),
            'detailed_error_breakdown': self.error_counts
        })
        
        return progress_data
    
    def save_report(self, file_path: str) -> None:
        """
        Save a progress report to a file.
        
        Args:
            file_path: Path to save the report
        """
        report = self.generate_report()
        
        try:
            with open(file_path, 'w') as f:
                json.dump(report, f, indent=2)
                
            self.logger.info(f"Progress report saved to {file_path}")
        except IOError as e:
            self.logger.error(f"Failed to save progress report: {str(e)}")
    
    def reset(self) -> None:
        """Reset all progress metrics."""
        self.processed_items = 0
        self.successful_items = 0
        self.failed_items = 0
        self.skipped_items = 0
        
        self.start_time = time.time()
        self.last_update_time = self.start_time
        
        self.processing_times.clear()
        self.batch_sizes.clear()
        self.error_counts.clear()


class ConsoleProgressReporter:
    """
    Report progress to the console with periodic updates.
    
    This class connects to a ProgressTracker and displays updates
    at regular intervals.
    """
    
    def __init__(
        self,
        tracker: ProgressTracker,
        update_interval: float = 5.0,
        logger: Optional[logging.Logger] = None
    ):
        """
        Initialize console progress reporter.
        
        Args:
            tracker: ProgressTracker instance to report on
            update_interval: Minimum time between updates (seconds)
            logger: Logger instance for progress reports
        """
        self.tracker = tracker
        self.update_interval = update_interval
        self.logger = logger or logging.getLogger(__name__)
        
        self.last_update_time = 0
        
        # Register with tracker
        self.tracker.register_callback(self._progress_callback)
    
    def _progress_callback(self, progress_data: Dict[str, Any]) -> None:
        """
        Process progress updates from the tracker.
        
        Args:
            progress_data: Progress data dictionary
        """
        current_time = time.time()
        
        # Only update console at specified interval
        if current_time - self.last_update_time >= self.update_interval:
            self._display_progress(progress_data)
            self.last_update_time = current_time
    
    def _display_progress(self, progress_data: Dict[str, Any]) -> None:
        """
        Display progress information.
        
        Args:
            progress_data: Progress data dictionary
        """
        # Log the progress at INFO level
        self.tracker.log_progress(logging.INFO)
    
    def stop(self) -> None:
        """Stop the reporter and unregister from tracker."""
        self.tracker.unregister_callback(self._progress_callback)