#!/usr/bin/env python3
"""
Worker Pool Implementation for Parallel PubChem Data Processing

This module provides a thread-based worker pool implementation for parallel processing
of PubChem data import tasks. It includes:

- WorkerPool: Manages a configurable number of worker threads
- Worker: Processes chunks from the work queue
- ResultCollector: Aggregates results from all workers

The implementation is designed to be thread-safe and efficiently coordinate shared
resources such as rate limiters and circuit breakers.
"""

import os
import time
import logging
import threading
import queue
import random
from typing import List, Dict, Any, Tuple, Optional, Callable, Union, Set
from datetime import datetime
from dataclasses import dataclass, field
import traceback

# Set up logging
logger = logging.getLogger(__name__)

# Default configuration
DEFAULT_NUM_WORKERS = 4
DEFAULT_MAX_QUEUE_SIZE = 100
DEFAULT_RESULT_QUEUE_SIZE = 1000
DEFAULT_WORKER_TIMEOUT = 60  # seconds
DEFAULT_HEALTH_CHECK_INTERVAL = 5  # seconds


@dataclass
class WorkItem:
    """Represents a unit of work to be processed by a worker."""
    id: str
    data: Any
    priority: int = 0
    created_at: float = field(default_factory=time.time)
    attempts: int = 0
    max_attempts: int = 3
    
    def increment_attempt(self) -> None:
        """Increment the attempt counter."""
        self.attempts += 1
    
    def can_retry(self) -> bool:
        """Check if the work item can be retried."""
        return self.attempts < self.max_attempts


@dataclass
class WorkResult:
    """Represents the result of processing a work item."""
    work_item_id: str
    success: bool
    result: Any = None
    error: Optional[Exception] = None
    processing_time: float = 0.0
    worker_id: Optional[int] = None
    completed_at: float = field(default_factory=time.time)