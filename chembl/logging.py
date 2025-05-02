"""
Centralized, structured logging system for ChEMBL data import.

This module provides a robust logging system for the ChEMBL data import process,
featuring structured JSON logging, automatic log directory creation, and
prevention of duplicate log entries.
"""

import os
import json
import logging
import traceback
import uuid
import socket
import platform
import time
from datetime import datetime
from typing import Dict, Any, Optional, List, Union, Tuple
from pathlib import Path

# Try to import psutil for system metrics, but don't fail if not available
try:
    import psutil
    PSUTIL_AVAILABLE = True
except ImportError:
    PSUTIL_AVAILABLE = False

class ChEMBLLogger:
    """
    Centralized logging system for ChEMBL data import.
    
    Features:
    - Structured JSON logging for progress, errors, and skipped molecules
    - Automatic log directory creation
    - Prevention of duplicate log entries
    - Logging to both file and console
    - System metrics collection (if psutil is available)
    """
    
    def __init__(
        self,
        log_dir: str = "logs",
        progress_log: str = "chembl_progress.jsonl",
        error_log: str = "chembl_errors.jsonl",
        skipped_log: str = "skipped_chembl_molecules.jsonl",
        summary_log: str = "chembl_summary.json",
        general_log: str = "chembl_import.log",
        console_level: int = logging.INFO,
        file_level: int = logging.DEBUG
    ):
        """
        Initialize the logging system.
        
        Args:
            log_dir: Directory to store log files
            progress_log: Filename for progress log
            error_log: Filename for error log
            skipped_log: Filename for skipped molecules log
            summary_log: Filename for summary log
            general_log: Filename for general log
            console_level: Logging level for console output
            file_level: Logging level for file output
        """
        self.log_dir = Path(log_dir)
        self.progress_log_path = self.log_dir / progress_log
        self.error_log_path = self.log_dir / error_log
        self.skipped_log_path = self.log_dir / skipped_log
        self.summary_log_path = self.log_dir / summary_log
        self.general_log_path = self.log_dir / general_log
        
        # Create log directory if it doesn't exist
        self._ensure_log_directory()
        
        # Set up general logger
        self.logger = logging.getLogger("chembl_import")
        self.logger.setLevel(logging.DEBUG)  # Set to lowest level, handlers will filter
        self.logger.propagate = False  # Don't propagate to root logger
        
        # Remove any existing handlers to prevent duplicates
        if self.logger.handlers:
            self.logger.handlers.clear()
        
        # Add console handler
        console_handler = logging.StreamHandler()
        console_handler.setLevel(console_level)
        console_formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
        console_handler.setFormatter(console_formatter)
        self.logger.addHandler(console_handler)
        
        # Add file handler
        file_handler = logging.FileHandler(self.general_log_path)
        file_handler.setLevel(file_level)
        file_formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(name)s - %(message)s")
        file_handler.setFormatter(file_formatter)
        self.logger.addHandler(file_handler)
        
        # Set up progress logger
        self.progress_logger = logging.getLogger("chembl_import.progress")
        self.progress_logger.setLevel(logging.INFO)
        self.progress_logger.propagate = False  # Don't propagate to parent logger
        
        # Remove any existing handlers to prevent duplicates
        if self.progress_logger.handlers:
            self.progress_logger.handlers.clear()
        
        # Add file handler for progress logs
        progress_file_handler = logging.FileHandler(self.progress_log_path)
        progress_file_handler.setFormatter(logging.Formatter('%(message)s'))
        self.progress_logger.addHandler(progress_file_handler)
        
        # Add console handler for progress logs (optional, for debugging)
        progress_console_handler = logging.StreamHandler()
        progress_console_handler.setFormatter(logging.Formatter('%(message)s'))
        self.progress_logger.addHandler(progress_console_handler)
        
        self.logger.info(f"Logging system initialized. Log directory: {self.log_dir}")
        self.logger.info(f"Progress log: {self.progress_log_path}")
        self.logger.info(f"Error log: {self.error_log_path}")
        self.logger.info(f"Skipped molecules log: {self.skipped_log_path}")
    
    def _ensure_log_directory(self) -> None:
        """Create log directory if it doesn't exist."""
        if not self.log_dir.exists():
            self.log_dir.mkdir(parents=True, exist_ok=True)
            print(f"Created log directory: {self.log_dir}")
    
    def log_error(
        self,
        error_type: str,
        message: str,
        context: Dict[str, Any],
        severity: str = "error",
        notify: bool = False
    ) -> str:
        """
        Log a structured error in JSON lines format with full context and traceability.
        
        Args:
            error_type: Type of error (e.g., 'API', 'Database', 'Validation')
            message: Error message
            context: Additional context information
            severity: Error severity level ('error', 'warning', 'critical')
            notify: Whether this error should trigger notifications
            
        Returns:
            str: Unique error ID for reference
        """
        try:
            # Create a safe copy of context with non-serializable objects converted to strings
            safe_context = {}
            for key, value in context.items():
                if key == "exception":
                    # Handle exception separately
                    continue
                try:
                    # Test if value is JSON serializable
                    json.dumps(value)
                    safe_context[key] = value
                except (TypeError, OverflowError):
                    # Convert non-serializable objects to string representation
                    safe_context[key] = str(value)
            
            # Generate a unique ID for this error
            error_id = str(uuid.uuid4())
            
            # Prepare error record with enhanced metadata
            error_record = {
                "error_id": error_id,
                "timestamp": datetime.now().isoformat(),
                "error_type": error_type,
                "severity": severity,
                "message": message,
                "context": safe_context,
                "environment": {
                    "hostname": socket.gethostname(),
                    "platform": platform.platform(),
                    "python_version": platform.python_version(),
                    "process_id": os.getpid()
                },
                "notify": notify
            }
            
            # Add memory usage if psutil is available
            if PSUTIL_AVAILABLE:
                error_record["environment"]["memory_percent"] = psutil.virtual_memory().percent
            
            # Add stack trace if available
            if "exception" in context:
                exc = context["exception"]
                if exc:
                    try:
                        error_record["stack_trace"] = traceback.format_exception(
                            type(exc), exc, exc.__traceback__
                        )
                        # Add exception type for easier filtering
                        error_record["exception_type"] = type(exc).__name__
                    except Exception as trace_error:
                        # Fallback to simple string representation if traceback fails
                        error_record["exception_info"] = str(exc)
                        error_record["exception_type"] = type(exc).__name__
                        self.logger.warning(f"Could not extract stack trace: {str(trace_error)}")
            
            # Write to log file in JSON lines format with atomic operation
            self._write_json_line(self.error_log_path, error_record)
            
            # Also log to standard logger with appropriate severity
            log_method = getattr(self.logger, severity, self.logger.error)
            log_method(f"[{error_id}] {error_type} {severity}: {message}")
            
            return error_id
        
        except Exception as e:
            # Fallback to standard logger if structured logging fails
            self.logger.error(f"Failed to log structured error: {str(e)}")
            self.logger.error(f"Original error - {error_type} {severity}: {message}")
            return str(uuid.uuid4())  # Still return a unique ID
    
    def log_skipped_molecule(
        self,
        chembl_id: str,
        reason: str,
        molecule_data: Optional[Dict[str, Any]] = None,
        category: str = "filter",
        batch_num: Optional[int] = None
    ) -> str:
        """
        Log a skipped molecule in JSON lines format with enhanced traceability and categorization.
        
        Args:
            chembl_id: ChEMBL ID of the skipped molecule
            reason: Reason for skipping
            molecule_data: Additional molecule data if available
            category: Category of skip reason (filter, duplicate, error, validation)
            batch_num: Batch number where the molecule was skipped
            
        Returns:
            str: Unique ID for the skipped record for traceability
        """
        try:
            # Generate a unique ID for this skipped record
            skip_id = str(uuid.uuid4())
            
            # Prepare skipped record with enhanced metadata
            skipped_record = {
                "skip_id": skip_id,
                "timestamp": datetime.now().isoformat(),
                "chembl_id": chembl_id,
                "reason": reason,
                "category": category,
                "batch_num": batch_num,
                "process_id": os.getpid()
            }
            
            # Add molecule data if available
            if molecule_data:
                # Filter out potentially large or sensitive fields
                filtered_data = {k: v for k, v in molecule_data.items()
                               if k not in ["InChI"] and v is not None}
                
                # Add key properties for easier analysis
                key_properties = {}
                for key in ["Name", "Molecular Weight", "LogP", "TPSA", "InChIKey"]:
                    if key in filtered_data:
                        key_properties[key] = filtered_data[key]
                
                skipped_record["key_properties"] = key_properties
                skipped_record["molecule_data"] = filtered_data
            
            # Add stack trace if this was due to an exception
            if molecule_data and "exception" in molecule_data or category == "error":
                exc = molecule_data.get("exception") if molecule_data else None
                if exc:
                    skipped_record["stack_trace"] = traceback.format_exception(
                        type(exc), exc, exc.__traceback__
                    )
                    skipped_record["exception_type"] = type(exc).__name__
            
            # Write to log file in JSON lines format with atomic operation
            self._write_json_line(self.skipped_log_path, skipped_record)
            
            # Also log to standard logger with appropriate category
            self.logger.warning(f"SKIPPED [{category}]: {chembl_id} - {reason} (ID: {skip_id})")
            
            # Return the skip ID for reference
            return skip_id
        
        except Exception as e:
            # Fallback to standard logger if structured logging fails
            self.logger.error(f"Failed to log skipped molecule: {str(e)}")
            self.logger.warning(f"SKIPPED: {chembl_id} - {reason}")
            return f"error_{int(time.time())}"
    
    def log_progress(
        self,
        batch_num: int,
        total_batches: int,
        total_processed: int,
        total_ids: int,
        total_imported: int,
        total_duplicates: int,
        skipped_in_batch: int,
        eta: str,
        memory_info: str = "",
        additional_data: Dict[str, Any] = None
    ) -> Dict[str, Any]:
        """
        Log standardized progress information after each batch with enhanced metrics and traceability.
        
        Args:
            batch_num: Current batch number (0-indexed)
            total_batches: Total number of batches
            total_processed: Total number of molecules processed so far
            total_ids: Total number of molecules to process
            total_imported: Total number of molecules imported so far
            total_duplicates: Total number of duplicates found so far
            skipped_in_batch: Number of molecules skipped in this batch
            eta: Estimated time remaining
            memory_info: Memory usage information (optional)
            additional_data: Any additional data to include in the log (optional)
            
        Returns:
            Dict[str, Any]: The progress record for reference
        """
        # Create standardized progress record with enhanced metrics
        timestamp = datetime.now().isoformat()
        
        # Generate a unique ID for this progress record
        progress_id = str(uuid.uuid4())
        
        # Calculate rates and additional metrics
        elapsed_time = time.time() - (additional_data.get("start_time", 0) if additional_data else 0)
        if elapsed_time > 0 and batch_num > 0:
            molecules_per_second = total_processed / elapsed_time
            molecules_per_batch = total_processed / (batch_num + 1)
            success_rate = (total_imported / total_processed * 100) if total_processed > 0 else 0
        else:
            molecules_per_second = 0
            molecules_per_batch = 0
            success_rate = 0
        
        # Get system metrics
        system_metrics = {"error": "psutil not available"}
        if PSUTIL_AVAILABLE:
            try:
                system_metrics = {
                    "cpu_percent": psutil.cpu_percent(interval=0.1),
                    "memory_percent": psutil.virtual_memory().percent,
                    "disk_usage_percent": psutil.disk_usage('/').percent
                }
            except Exception as e:
                system_metrics = {"error": str(e)}
        
        progress_record = {
            "progress_id": progress_id,
            "timestamp": timestamp,
            "event_type": "batch_progress",
            "batch": {
                "current": batch_num + 1,
                "total": total_batches,
                "percent_complete": round((batch_num + 1) / total_batches * 100, 2) if total_batches > 0 else 0
            },
            "molecules": {
                "processed": total_processed,
                "total": total_ids,
                "percent_complete": round(total_processed / total_ids * 100, 2) if total_ids > 0 else 0,
                "imported": total_imported,
                "duplicates": total_duplicates,
                "skipped_in_batch": skipped_in_batch,
                "total_skipped": additional_data.get("total_skipped", skipped_in_batch) if additional_data else skipped_in_batch
            },
            "performance": {
                "elapsed_seconds": round(elapsed_time, 2),
                "molecules_per_second": round(molecules_per_second, 2),
                "molecules_per_batch": round(molecules_per_batch, 2),
                "success_rate_percent": round(success_rate, 2)
            },
            "system": system_metrics,
            "eta": eta,
            "memory": memory_info.strip() if memory_info else None,
            "process_id": os.getpid()
        }
        
        # Add any additional data
        if additional_data:
            # Don't overwrite existing fields
            for key, value in additional_data.items():
                if key not in progress_record:
                    progress_record[key] = value
        
        # Write to progress log file with atomic operation
        self._write_json_line(self.progress_log_path, progress_record)
        
        # Also log a human-readable message to the standard logger
        self.logger.info(
            f"Batch {batch_num+1}/{total_batches} complete: {total_processed}/{total_ids} ChEMBL IDs processed, "
            f"{total_imported} imported ({round(success_rate, 1)}% success), {total_duplicates} duplicates, "
            f"{skipped_in_batch} skipped in this batch. Rate: {round(molecules_per_second, 2)}/sec. "
            f"ETA: {eta}.{memory_info}"
        )
        
        return progress_record
    
    def write_summary(self, summary_data: Dict[str, Any]) -> None:
        """
        Write a summary report of the import process.
        
        Args:
            summary_data: Summary data to write
        """
        try:
            # Ensure the summary data has a timestamp
            if "timestamp" not in summary_data:
                summary_data["timestamp"] = datetime.now().isoformat()
            
            # Write the summary to a JSON file
            with open(self.summary_log_path, 'w') as f:
                json.dump(summary_data, f, indent=2)
            
            self.logger.info(f"Summary report written to {self.summary_log_path}")
        except Exception as e:
            self.logger.error(f"Failed to write summary report: {str(e)}")
    
    def _write_json_line(self, file_path: Path, data: Dict[str, Any]) -> None:
        """
        Write a JSON line to a file with atomic operation.
        
        Args:
            file_path: Path to the file
            data: Data to write
        """
        try:
            # Create a temporary file
            temp_path = f"{file_path}.tmp"
            
            # Write to the temporary file
            with open(temp_path, "w") as f:
                try:
                    f.write(json.dumps(data) + "\n")
                except (TypeError, OverflowError) as json_error:
                    # If JSON serialization fails, create a simplified record
                    simplified_record = {
                        "timestamp": data.get("timestamp", datetime.now().isoformat()),
                        "error": "JSON serialization failed",
                        "error_details": str(json_error)
                    }
                    f.write(json.dumps(simplified_record) + "\n")
            
            # Append to the log file (atomic operation on most file systems)
            with open(file_path, "a") as f:
                with open(temp_path, "r") as temp_f:
                    f.write(temp_f.read())
            
            # Remove temporary file
            os.remove(temp_path)
        
        except Exception as e:
            # Fall back to direct write if atomic operation fails
            self.logger.error(f"Error writing to temporary log file: {str(e)}")
            try:
                with open(file_path, "a") as f:
                    f.write(json.dumps(data) + "\n")
            except Exception as direct_write_error:
                self.logger.error(f"Failed to write log entry: {str(direct_write_error)}")

# Create a default logger instance
default_logger = ChEMBLLogger()

# Convenience functions that use the default logger
def log_error(error_type: str, message: str, context: Dict[str, Any], 
              severity: str = "error", notify: bool = False) -> str:
    """Convenience function to log an error using the default logger."""
    return default_logger.log_error(error_type, message, context, severity, notify)

def log_skipped_molecule(chembl_id: str, reason: str, 
                         molecule_data: Optional[Dict[str, Any]] = None,
                         category: str = "filter", batch_num: Optional[int] = None) -> str:
    """Convenience function to log a skipped molecule using the default logger."""
    return default_logger.log_skipped_molecule(chembl_id, reason, molecule_data, category, batch_num)

def log_progress(batch_num: int, total_batches: int, total_processed: int,
                total_ids: int, total_imported: int, total_duplicates: int,
                skipped_in_batch: int, eta: str, memory_info: str = "",
                additional_data: Dict[str, Any] = None) -> Dict[str, Any]:
    """Convenience function to log progress using the default logger."""
    return default_logger.log_progress(
        batch_num, total_batches, total_processed, total_ids,
        total_imported, total_duplicates, skipped_in_batch,
        eta, memory_info, additional_data
    )

def write_summary(summary_data: Dict[str, Any]) -> None:
    """Convenience function to write a summary using the default logger."""
    default_logger.write_summary(summary_data)

def get_logger(name: str = "chembl_import") -> logging.Logger:
    """Get a named logger that inherits from the default logger's configuration."""
    return logging.getLogger(name)