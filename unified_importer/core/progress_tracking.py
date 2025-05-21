"""
Enhanced progress tracking system for CryoProtect Unified Importer.

This module provides comprehensive progress tracking capabilities with ETA 
calculation, statistics, and multiple reporting formats.
"""

import time
import json
import logging
import datetime
import math
from enum import Enum
from typing import Dict, Any, Optional, List, Union, Callable
from dataclasses import dataclass, field

class ReportFormat(Enum):
    """Available formats for progress reporting."""
    CONSOLE = "console"  # Plain text format for console output
    JSON = "json"        # JSON format for structured data
    CSV = "csv"          # CSV format for data export
    DETAILED = "detailed"  # Detailed text report
    MINIMAL = "minimal"  # Minimal text report with essential info only

@dataclass
class ProgressStats:
    """Statistics for a progress tracking session."""
    total_items: int = 0
    processed_items: int = 0
    success_count: int = 0
    error_count: int = 0
    warning_count: int = 0
    skipped_count: int = 0
    error_categories: Dict[str, int] = field(default_factory=dict)
    
    # Performance metrics
    start_time: float = field(default_factory=time.time)
    end_time: Optional[float] = None
    last_update_time: float = field(default_factory=time.time)
    elapsed_time: float = 0.0
    items_per_second: float = 0.0
    average_process_time: float = 0.0
    recent_process_times: List[float] = field(default_factory=list)
    
    # ETA calculation
    estimated_completion_time: Optional[float] = None
    estimated_time_remaining: float = 0.0
    estimated_completion_datetime: Optional[datetime.datetime] = None
    eta_last_calculated: float = field(default_factory=time.time)
    eta_update_interval: float = 1.0  # Update ETA every second by default
    
    # Adaptive smoothing for ETA
    eta_smoothing_factor: float = 0.2  # Initial smoothing factor
    eta_samples: List[float] = field(default_factory=list)
    eta_max_samples: int = 10
    eta_variance: float = 0.0
    
    # Checkpoints
    checkpoint_times: List[Dict[str, Any]] = field(default_factory=list)
    
    def update_elapsed_time(self) -> None:
        """Update the elapsed time for the progress tracking session."""
        current_time = time.time()
        self.last_update_time = current_time
        if self.end_time:
            self.elapsed_time = self.end_time - self.start_time
        else:
            self.elapsed_time = current_time - self.start_time
    
    def update_processing_rate(self) -> None:
        """Update the processing rate statistics."""
        if self.elapsed_time > 0:
            self.items_per_second = self.processed_items / self.elapsed_time
            
            if self.processed_items > 0:
                self.average_process_time = self.elapsed_time / self.processed_items
    
    def calculate_eta(self) -> None:
        """Calculate the estimated time of completion (ETA)."""
        current_time = time.time()
        
        # Only update ETA if enough time has passed since last calculation
        if current_time - self.eta_last_calculated < self.eta_update_interval:
            return
        
        self.eta_last_calculated = current_time
        
        # Update elapsed time and processing rate
        self.update_elapsed_time()
        self.update_processing_rate()
        
        # Calculate remaining items
        remaining_items = max(0, self.total_items - self.processed_items)
        
        if self.items_per_second > 0 and remaining_items > 0:
            # Simple ETA calculation
            simple_eta = remaining_items / self.items_per_second
            
            # Add to samples for adaptive calculation
            self.eta_samples.append(simple_eta)
            if len(self.eta_samples) > self.eta_max_samples:
                self.eta_samples.pop(0)
            
            # Calculate variance of ETA samples
            if len(self.eta_samples) >= 3:
                mean_eta = sum(self.eta_samples) / len(self.eta_samples)
                squared_diffs = [(eta - mean_eta) ** 2 for eta in self.eta_samples]
                self.eta_variance = sum(squared_diffs) / len(self.eta_samples)
                
                # Adjust smoothing factor based on variance
                # Higher variance = lower smoothing factor (more responsive)
                variance_factor = 1.0 / (1.0 + math.sqrt(self.eta_variance))
                self.eta_smoothing_factor = max(0.1, min(0.9, variance_factor))
            
            # Apply exponential smoothing
            if self.estimated_time_remaining > 0:
                self.estimated_time_remaining = (
                    self.eta_smoothing_factor * simple_eta +
                    (1 - self.eta_smoothing_factor) * self.estimated_time_remaining
                )
            else:
                self.estimated_time_remaining = simple_eta
            
            # Calculate estimated completion time
            self.estimated_completion_time = current_time + self.estimated_time_remaining
            self.estimated_completion_datetime = datetime.datetime.fromtimestamp(
                self.estimated_completion_time
            )
    
    def add_checkpoint(self, name: str, data: Dict[str, Any] = None) -> None:
        """Add a checkpoint to track progress milestones.
        
        Args:
            name: Name of the checkpoint
            data: Additional data to store with the checkpoint
        """
        checkpoint = {
            "name": name,
            "time": time.time(),
            "elapsed": time.time() - self.start_time,
            "processed_items": self.processed_items,
            "success_count": self.success_count,
            "error_count": self.error_count,
            "data": data or {}
        }
        self.checkpoint_times.append(checkpoint)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert the progress stats to a dictionary.
        
        Returns:
            Dictionary with progress statistics
        """
        # Update time-based metrics before returning
        self.update_elapsed_time()
        self.update_processing_rate()
        
        result = {
            "total_items": self.total_items,
            "processed_items": self.processed_items,
            "percent_complete": (self.processed_items / self.total_items * 100) if self.total_items > 0 else 0,
            "success_count": self.success_count,
            "error_count": self.error_count,
            "warning_count": self.warning_count,
            "skipped_count": self.skipped_count,
            "error_categories": self.error_categories,
            "elapsed_time": self.elapsed_time,
            "elapsed_time_formatted": format_time_duration(self.elapsed_time),
            "items_per_second": self.items_per_second,
            "average_process_time": self.average_process_time,
            "start_time": self.start_time,
            "start_time_formatted": format_timestamp(self.start_time),
        }
        
        # Add ETA if available
        if self.estimated_completion_time:
            result.update({
                "estimated_time_remaining": self.estimated_time_remaining,
                "estimated_time_remaining_formatted": format_time_duration(self.estimated_time_remaining),
                "estimated_completion_time": self.estimated_completion_time,
                "estimated_completion_time_formatted": format_timestamp(self.estimated_completion_time),
                "eta_variance": self.eta_variance,
                "eta_smoothing_factor": self.eta_smoothing_factor
            })
        
        # Add end time if finished
        if self.end_time:
            result.update({
                "end_time": self.end_time,
                "end_time_formatted": format_timestamp(self.end_time),
                "total_duration": self.end_time - self.start_time,
                "total_duration_formatted": format_time_duration(self.end_time - self.start_time)
            })
        
        # Add checkpoint info
        if self.checkpoint_times:
            result["checkpoints"] = self.checkpoint_times
            
            # Calculate time between checkpoints
            checkpoint_durations = []
            for i in range(1, len(self.checkpoint_times)):
                prev = self.checkpoint_times[i-1]
                curr = self.checkpoint_times[i]
                duration = curr["time"] - prev["time"]
                checkpoint_durations.append({
                    "from": prev["name"],
                    "to": curr["name"],
                    "duration": duration,
                    "duration_formatted": format_time_duration(duration),
                    "items_processed": curr["processed_items"] - prev["processed_items"]
                })
            
            if checkpoint_durations:
                result["checkpoint_durations"] = checkpoint_durations
        
        return result
    
    def mark_complete(self) -> None:
        """Mark the progress tracking as complete."""
        self.end_time = time.time()
        self.update_elapsed_time()
        self.update_processing_rate()
        self.estimated_time_remaining = 0
        self.estimated_completion_time = self.end_time
        self.estimated_completion_datetime = datetime.datetime.fromtimestamp(self.end_time)
    
    def to_json(self) -> str:
        """Convert the progress stats to a JSON string.
        
        Returns:
            JSON string with progress statistics
        """
        return json.dumps(self.to_dict(), indent=2)
    
    def to_csv(self) -> str:
        """Convert the progress stats to a CSV string.
        
        Returns:
            CSV string with progress statistics
        """
        data = self.to_dict()
        # Flatten the data for CSV format
        flat_data = {
            "total_items": data["total_items"],
            "processed_items": data["processed_items"],
            "percent_complete": data["percent_complete"],
            "success_count": data["success_count"],
            "error_count": data["error_count"],
            "warning_count": data["warning_count"],
            "skipped_count": data["skipped_count"],
            "elapsed_time": data["elapsed_time"],
            "items_per_second": data["items_per_second"],
            "average_process_time": data["average_process_time"]
        }
        
        if "estimated_time_remaining" in data:
            flat_data["estimated_time_remaining"] = data["estimated_time_remaining"]
        
        if "end_time" in data:
            flat_data["total_duration"] = data["total_duration"]
        
        # Convert to CSV
        header = ",".join(flat_data.keys())
        values = ",".join(str(v) for v in flat_data.values())
        return f"{header}\n{values}"
    
    def to_console(self, detailed: bool = False) -> str:
        """Convert the progress stats to a console-friendly string.
        
        Args:
            detailed: Whether to include detailed information
            
        Returns:
            Formatted string for console output
        """
        data = self.to_dict()
        
        if detailed:
            lines = [
                f"Progress: {data['processed_items']}/{data['total_items']} items ({data['percent_complete']:.1f}%)",
                f"Status: {data['success_count']} successful, {data['error_count']} errors, {data['warning_count']} warnings, {data['skipped_count']} skipped",
                f"Time: {data['elapsed_time_formatted']} elapsed"
            ]
            
            if "estimated_time_remaining_formatted" in data:
                lines.append(f"ETA: {data['estimated_time_remaining_formatted']} remaining, completion at {data['estimated_completion_time_formatted']}")
            
            if "end_time_formatted" in data:
                lines.append(f"Completed: {data['end_time_formatted']}, total duration {data['total_duration_formatted']}")
            
            lines.append(f"Performance: {data['items_per_second']:.2f} items/sec, {data['average_process_time']*1000:.2f} ms/item")
            
            if data["error_count"] > 0 and data["error_categories"]:
                error_details = ", ".join(f"{cat}: {count}" for cat, count in data["error_categories"].items())
                lines.append(f"Error categories: {error_details}")
            
            if "checkpoints" in data:
                lines.append("\nCheckpoints:")
                for cp in data["checkpoints"]:
                    cp_time = format_timestamp(cp["time"])
                    cp_elapsed = format_time_duration(cp["elapsed"])
                    lines.append(f"  {cp['name']}: {cp_time} ({cp_elapsed} elapsed)")
            
            if "checkpoint_durations" in data:
                lines.append("\nCheckpoint Durations:")
                for dur in data["checkpoint_durations"]:
                    lines.append(f"  {dur['from']} → {dur['to']}: {dur['duration_formatted']} ({dur['items_processed']} items)")
            
            return "\n".join(lines)
        else:
            # Minimal output
            progress = f"{data['processed_items']}/{data['total_items']} ({data['percent_complete']:.1f}%)"
            status = f"{data['success_count']} ok, {data['error_count']} err"
            
            if "estimated_time_remaining_formatted" in data:
                eta = f"ETA: {data['estimated_time_remaining_formatted']}"
                return f"{progress} - {status} - {eta}"
            elif "end_time" in data:
                duration = f"Done in {data['total_duration_formatted']}"
                return f"{progress} - {status} - {duration}"
            else:
                elapsed = f"Elapsed: {data['elapsed_time_formatted']}"
                return f"{progress} - {status} - {elapsed}"

class ProgressTracker:
    """Enhanced progress tracking with ETA calculation and statistical reporting."""
    
    def __init__(
        self,
        total_items: int = 0,
        description: str = "",
        logger: Optional[logging.Logger] = None,
        update_interval: float = 1.0,
        checkpoint_file: Optional[str] = None,
        report_format: ReportFormat = ReportFormat.CONSOLE,
        on_update: Optional[Callable[[ProgressStats], None]] = None
    ):
        """Initialize the progress tracker.
        
        Args:
            total_items: Total number of items to process
            description: Description of the progress tracking session
            logger: Optional logger instance
            update_interval: Interval for progress updates in seconds
            checkpoint_file: Optional file to save progress checkpoints
            report_format: Format for progress reports
            on_update: Optional callback function for progress updates
        """
        self.description = description
        self.logger = logger or logging.getLogger(__name__)
        self.checkpoint_file = checkpoint_file
        self.report_format = report_format
        self.update_interval = update_interval
        self.on_update = on_update
        
        self.stats = ProgressStats(total_items=total_items)
        self.last_report_time = time.time()
    
    def start(self, total_items: Optional[int] = None) -> None:
        """Start or restart the progress tracking.
        
        Args:
            total_items: Optional updated total items count
        """
        if total_items is not None:
            self.stats.total_items = total_items
        
        self.stats.start_time = time.time()
        self.stats.end_time = None
        self.stats.processed_items = 0
        self.stats.success_count = 0
        self.stats.error_count = 0
        self.stats.warning_count = 0
        self.stats.skipped_count = 0
        self.stats.error_categories = {}
        self.stats.eta_samples = []
        self.stats.checkpoint_times = []
        
        self.last_report_time = time.time()
        
        # Add initial checkpoint
        self.stats.add_checkpoint("start")
        
        # Log the start of tracking
        self.logger.info(
            f"Started progress tracking: {self.description} - "
            f"{self.stats.total_items} items to process"
        )
        
        # Call update callback if provided
        if self.on_update:
            self.on_update(self.stats)
    
    def update(
        self,
        increment: int = 1,
        success: int = 0,
        error: int = 0,
        warning: int = 0,
        skipped: int = 0,
        error_category: Optional[str] = None,
        force_report: bool = False
    ) -> None:
        """Update the progress statistics.
        
        Args:
            increment: Number of items processed in this update
            success: Number of successful items
            error: Number of error items
            warning: Number of warning items
            skipped: Number of skipped items
            error_category: Optional category for errors
            force_report: Whether to force a progress report
        """
        # Update counts
        self.stats.processed_items += increment
        self.stats.success_count += success
        self.stats.error_count += error
        self.stats.warning_count += warning
        self.stats.skipped_count += skipped
        
        # Update error categories
        if error > 0 and error_category:
            self.stats.error_categories[error_category] = (
                self.stats.error_categories.get(error_category, 0) + error
            )
        
        # Calculate ETA
        self.stats.calculate_eta()
        
        # Check if it's time to report progress
        current_time = time.time()
        if force_report or current_time - self.last_report_time >= self.update_interval:
            self.report_progress()
            self.last_report_time = current_time
        
        # Call update callback if provided
        if self.on_update:
            self.on_update(self.stats)
    
    def add_checkpoint(self, name: str, data: Dict[str, Any] = None) -> None:
        """Add a checkpoint to track progress milestones.
        
        Args:
            name: Name of the checkpoint
            data: Additional data to store with the checkpoint
        """
        self.stats.add_checkpoint(name, data)
        
        # Log the checkpoint
        self.logger.info(
            f"Checkpoint '{name}': {self.stats.processed_items}/{self.stats.total_items} items processed "
            f"({self.stats.processed_items / self.stats.total_items * 100:.1f}% complete)"
        )
        
        # Save checkpoint file if specified
        if self.checkpoint_file:
            self._save_checkpoint()
        
        # Call update callback if provided
        if self.on_update:
            self.on_update(self.stats)
    
    def _save_checkpoint(self) -> None:
        """Save the current progress to a checkpoint file."""
        if not self.checkpoint_file:
            return
        
        checkpoint_data = {
            "description": self.description,
            "timestamp": time.time(),
            "stats": self.stats.to_dict()
        }
        
        try:
            with open(self.checkpoint_file, 'w') as f:
                json.dump(checkpoint_data, f, indent=2)
            
            self.logger.debug(f"Saved progress checkpoint to {self.checkpoint_file}")
        except Exception as e:
            self.logger.error(f"Failed to save checkpoint: {e}")
    
    def load_checkpoint(self) -> bool:
        """Load progress from a checkpoint file.
        
        Returns:
            True if checkpoint was loaded successfully, False otherwise
        """
        if not self.checkpoint_file or not os.path.exists(self.checkpoint_file):
            return False
        
        try:
            with open(self.checkpoint_file, 'r') as f:
                checkpoint_data = json.load(f)
            
            # Extract basic statistics
            stats = checkpoint_data.get("stats", {})
            self.description = checkpoint_data.get("description", self.description)
            self.stats.total_items = stats.get("total_items", self.stats.total_items)
            self.stats.processed_items = stats.get("processed_items", 0)
            self.stats.success_count = stats.get("success_count", 0)
            self.stats.error_count = stats.get("error_count", 0)
            self.stats.warning_count = stats.get("warning_count", 0)
            self.stats.skipped_count = stats.get("skipped_count", 0)
            self.stats.error_categories = stats.get("error_categories", {})
            self.stats.checkpoint_times = stats.get("checkpoints", [])
            
            # Add a resume checkpoint
            self.stats.add_checkpoint("resume")
            
            self.logger.info(
                f"Resumed progress from checkpoint: {self.stats.processed_items}/{self.stats.total_items} "
                f"items processed ({self.stats.processed_items / self.stats.total_items * 100:.1f}% complete)"
            )
            
            # Call update callback if provided
            if self.on_update:
                self.on_update(self.stats)
            
            return True
        except Exception as e:
            self.logger.error(f"Failed to load checkpoint: {e}")
            return False
    
    def report_progress(self, format: Optional[ReportFormat] = None) -> str:
        """Generate a progress report.
        
        Args:
            format: Optional format override
            
        Returns:
            Formatted progress report
        """
        report_format = format or self.report_format
        
        # Generate the report in the requested format
        if report_format == ReportFormat.JSON:
            report = self.stats.to_json()
        elif report_format == ReportFormat.CSV:
            report = self.stats.to_csv()
        elif report_format == ReportFormat.DETAILED:
            report = self.stats.to_console(detailed=True)
        elif report_format == ReportFormat.MINIMAL:
            report = self.stats.to_console(detailed=False)
        else:  # Default to console format
            report = self.stats.to_console(detailed=False)
        
        # Log progress at info level
        if report_format in (ReportFormat.CONSOLE, ReportFormat.MINIMAL, ReportFormat.DETAILED):
            self.logger.info(f"Progress: {report}")
        else:
            # For structured formats, log a simple progress message
            self.logger.info(
                f"Progress: {self.stats.processed_items}/{self.stats.total_items} "
                f"({self.stats.processed_items / self.stats.total_items * 100:.1f}%)"
            )
        
        return report
    
    def complete(self) -> None:
        """Mark the progress tracking as complete."""
        self.stats.mark_complete()
        
        # Add final checkpoint
        self.stats.add_checkpoint("complete")
        
        # Generate final report
        report = self.stats.to_console(detailed=True)
        
        # Log completion
        self.logger.info(f"Completed progress tracking: {self.description}")
        self.logger.info(f"Final report:\n{report}")
        
        # Save final checkpoint
        if self.checkpoint_file:
            self._save_checkpoint()
        
        # Call update callback if provided
        if self.on_update:
            self.on_update(self.stats)

def format_time_duration(seconds: float) -> str:
    """Format a time duration in seconds to a human-readable string.
    
    Args:
        seconds: Duration in seconds
        
    Returns:
        Formatted duration string
    """
    if seconds < 0.001:  # Less than 1 millisecond
        return f"{seconds*1000000:.2f}μs"
    elif seconds < 1:  # Less than 1 second
        return f"{seconds*1000:.2f}ms"
    elif seconds < 60:  # Less than 1 minute
        return f"{seconds:.2f}s"
    elif seconds < 3600:  # Less than 1 hour
        minutes = int(seconds / 60)
        secs = seconds % 60
        return f"{minutes}m {secs:.1f}s"
    elif seconds < 86400:  # Less than 1 day
        hours = int(seconds / 3600)
        minutes = int((seconds % 3600) / 60)
        secs = seconds % 60
        return f"{hours}h {minutes}m {secs:.1f}s"
    else:  # Days or more
        days = int(seconds / 86400)
        hours = int((seconds % 86400) / 3600)
        minutes = int((seconds % 3600) / 60)
        return f"{days}d {hours}h {minutes}m"

def format_timestamp(timestamp: float) -> str:
    """Format a timestamp to a human-readable string.
    
    Args:
        timestamp: Unix timestamp
        
    Returns:
        Formatted timestamp string
    """
    dt = datetime.datetime.fromtimestamp(timestamp)
    return dt.strftime("%Y-%m-%d %H:%M:%S")