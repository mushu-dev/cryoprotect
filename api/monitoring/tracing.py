#!/usr/bin/env python3
"""
CryoProtect v2 - Distributed Tracing System

This module provides distributed tracing capabilities:
- Request tracing with correlation IDs across service boundaries
- Span management for tracking sub-operations within requests
- Propagation of trace context between services
- Visualization of request flows and dependencies
- Integration with OpenTelemetry standards
- Support for complex asynchronous and multi-service scenarios
- Context preservation across async boundaries

Usage:
    from api.monitoring.tracing import (
        Tracer,
        get_current_span,
        with_span,
        trace_function,
        extract_trace_context,
        inject_trace_context
    )
"""

import os
import uuid
import time
import json
import enum
import threading
import asyncio
import functools
import logging
import contextvars
from typing import Dict, List, Any, Optional, Union, Callable, TypeVar, ContextManager, Generator, cast

# Import structured logging
from logging_enhanced import log_with_context, get_logger

# Type variables
T = TypeVar('T')  # Return type for decorated functions

# Set up logger
logger = get_logger(__name__)

# Context variable to store the current span
current_span_var: contextvars.ContextVar = contextvars.ContextVar('current_span', default=None)


class SpanKind(enum.Enum):
    """Types of spans in a trace."""
    INTERNAL = "internal"        # Internal operation
    SERVER = "server"            # Server-side of an RPC
    CLIENT = "client"            # Client-side of an RPC
    PRODUCER = "producer"        # Producer of a message
    CONSUMER = "consumer"        # Consumer of a message
    BATCH = "batch"              # Batch processing operation


class Span:
    """
    Represents a single operation within a trace.
    
    A span represents a unit of work or operation within a trace. It can be
    nested to represent a hierarchical relationship between operations.
    """
    
    def __init__(
        self,
        name: str,
        trace_id: Optional[str] = None,
        parent_span_id: Optional[str] = None,
        span_kind: SpanKind = SpanKind.INTERNAL,
        attributes: Optional[Dict[str, Any]] = None,
        start_time: Optional[float] = None
    ):
        """
        Initialize a span.
        
        Args:
            name: Name of the operation
            trace_id: Trace ID (generated if None)
            parent_span_id: Parent span ID (None for root span)
            span_kind: Type of span
            attributes: Additional attributes for the span
            start_time: Start time in seconds (defaults to current time)
        """
        self.name = name
        self.trace_id = trace_id or str(uuid.uuid4())
        self.span_id = str(uuid.uuid4())
        self.parent_span_id = parent_span_id
        self.span_kind = span_kind
        self.attributes = attributes or {}
        self.events: List[Dict[str, Any]] = []
        self.links: List[Dict[str, Any]] = []
        self.status_code: Optional[str] = None
        self.status_message: Optional[str] = None
        self.start_time = start_time or time.time()
        self.end_time: Optional[float] = None
        self.duration: Optional[float] = None
        
        # Track active state
        self.is_active = True
        
        # Store the token for the context var when this span is set as current
        self._context_token = None
    
    def add_event(
        self,
        name: str,
        attributes: Optional[Dict[str, Any]] = None,
        timestamp: Optional[float] = None
    ) -> None:
        """
        Add an event to the span.
        
        Args:
            name: Event name
            attributes: Event attributes
            timestamp: Event timestamp (defaults to current time)
        """
        if not self.is_active:
            logger.warning(f"Adding event to inactive span: {self.name}")
        
        self.events.append({
            'name': name,
            'timestamp': timestamp or time.time(),
            'attributes': attributes or {}
        })
    
    def add_link(
        self,
        trace_id: str,
        span_id: str,
        attributes: Optional[Dict[str, Any]] = None
    ) -> None:
        """
        Add a link to another span.
        
        Args:
            trace_id: Trace ID of the linked span
            span_id: Span ID of the linked span
            attributes: Link attributes
        """
        if not self.is_active:
            logger.warning(f"Adding link to inactive span: {self.name}")
        
        self.links.append({
            'trace_id': trace_id,
            'span_id': span_id,
            'attributes': attributes or {}
        })
    
    def set_attribute(self, key: str, value: Any) -> None:
        """
        Set an attribute on the span.
        
        Args:
            key: Attribute key
            value: Attribute value
        """
        if not self.is_active:
            logger.warning(f"Setting attribute on inactive span: {self.name}")
        
        self.attributes[key] = value
    
    def set_status(self, code: str, message: Optional[str] = None) -> None:
        """
        Set the status of the span.
        
        Args:
            code: Status code (e.g., "ok", "error")
            message: Optional status message
        """
        if not self.is_active:
            logger.warning(f"Setting status on inactive span: {self.name}")
        
        self.status_code = code
        self.status_message = message
    
    def end(self, end_time: Optional[float] = None) -> None:
        """
        End the span.
        
        Args:
            end_time: End time in seconds (defaults to current time)
        """
        if not self.is_active:
            logger.warning(f"Ending already inactive span: {self.name}")
            return
        
        # Set end time and calculate duration
        self.end_time = end_time or time.time()
        self.duration = self.end_time - self.start_time
        
        # Mark as inactive
        self.is_active = False
        
        # Clear context token
        self._context_token = None
    
    def make_current(self) -> Optional:
        """
        Make this span the current span in the execution context.
        
        Returns:
            Token to restore the previous span when exiting the context
        """
        # Store the current span to restore later
        current = current_span_var.get()
        
        # Set this span as current
        self._context_token = current_span_var.set(self)
        
        return current
    
    def restore_previous(self, previous: Optional['Span']) -> None:
        """
        Restore the previous span as current.
        
        Args:
            previous: Previous span to restore
        """
        if self._context_token:
            current_span_var.reset(self._context_token)
            self._context_token = None
        elif previous:
            current_span_var.set(previous)
    
    def to_dict(self) -> Dict[str, Any]:
        """
        Convert the span to a dictionary.
        
        Returns:
            Dictionary representation of the span
        """
        return {
            'name': self.name,
            'trace_id': self.trace_id,
            'span_id': self.span_id,
            'parent_span_id': self.parent_span_id,
            'span_kind': self.span_kind.value,
            'attributes': self.attributes,
            'events': self.events,
            'links': self.links,
            'status_code': self.status_code,
            'status_message': self.status_message,
            'start_time': self.start_time,
            'end_time': self.end_time,
            'duration': self.duration,
            'is_active': self.is_active
        }


class Tracer:
    """
    Manages spans and traces.
    
    The Tracer is responsible for creating spans, managing the active span,
    and exporting completed spans to various backends.
    """
    
    _instance = None
    
    @classmethod
    def get_instance(cls) -> 'Tracer':
        """
        Get or create the singleton instance of Tracer.
        
        Returns:
            Tracer instance
        """
        if cls._instance is None:
            cls._instance = Tracer()
        return cls._instance
    
    def __init__(self):
        """Initialize the tracer."""
        # Ensure singleton pattern
        if Tracer._instance is not None:
            raise RuntimeError("Tracer is a singleton. Use get_instance() instead.")
        
        Tracer._instance = self
        
        # Storage for completed spans
        self.spans: Dict[str, Span] = {}
        
        # Storage for active spans
        self.active_spans: Dict[str, Span] = {}
        
        # Thread safety
        self.lock = threading.RLock()
        
        # Export configuration
        self.export_to_logs = True
        self.log_level = logging.INFO
        
        # Create logger
        self.logger = get_logger(__name__)
    
    def start_span(
        self,
        name: str,
        trace_id: Optional[str] = None,
        parent_span_id: Optional[str] = None,
        span_kind: SpanKind = SpanKind.INTERNAL,
        attributes: Optional[Dict[str, Any]] = None,
        start_time: Optional[float] = None
    ) -> Span:
        """
        Start a new span.
        
        Args:
            name: Name of the operation
            trace_id: Trace ID (generated if None)
            parent_span_id: Parent span ID (None for root span)
            span_kind: Type of span
            attributes: Additional attributes for the span
            start_time: Start time in seconds (defaults to current time)
            
        Returns:
            Created span
        """
        with self.lock:
            # Get current span for context if parent not specified
            if parent_span_id is None:
                current = get_current_span()
                if current:
                    parent_span_id = current.span_id
                    if trace_id is None:
                        trace_id = current.trace_id
            
            # Create span
            span = Span(
                name=name,
                trace_id=trace_id,
                parent_span_id=parent_span_id,
                span_kind=span_kind,
                attributes=attributes,
                start_time=start_time
            )
            
            # Store in active spans
            self.active_spans[span.span_id] = span
            
            return span
    
    def end_span(
        self,
        span: Span,
        end_time: Optional[float] = None,
        status_code: Optional[str] = None,
        status_message: Optional[str] = None
    ) -> None:
        """
        End a span.
        
        Args:
            span: Span to end
            end_time: End time in seconds (defaults to current time)
            status_code: Optional status code to set before ending
            status_message: Optional status message to set before ending
        """
        with self.lock:
            if not span.is_active:
                logger.warning(f"Ending already inactive span: {span.name}")
                return
            
            # Set status if provided
            if status_code:
                span.set_status(status_code, status_message)
            
            # End the span
            span.end(end_time)
            
            # Remove from active spans
            if span.span_id in self.active_spans:
                del self.active_spans[span.span_id]
            
            # Store the completed span
            self.spans[span.span_id] = span
            
            # Export the span
            self._export_span(span)
    
    def _export_span(self, span: Span) -> None:
        """
        Export a completed span.
        
        Args:
            span: Completed span to export
        """
        # Export to logs if enabled
        if self.export_to_logs:
            self._log_span(span)
    
    def _log_span(self, span: Span) -> None:
        """
        Log a completed span.
        
        Args:
            span: Completed span to log
        """
        # Create log message
        message = (
            f"Span '{span.name}' completed in {span.duration:.6f}s "
            f"(trace: {span.trace_id}, span: {span.span_id})"
        )
        
        # Log with structured context
        log_with_context(
            self.logger,
            logging.getLevelName(self.log_level).lower(),
            message,
            context={
                'tracing': {
                    'trace_id': span.trace_id,
                    'span_id': span.span_id,
                    'parent_span_id': span.parent_span_id,
                    'name': span.name,
                    'kind': span.span_kind.value,
                    'duration': span.duration,
                    'status': span.status_code,
                    'status_message': span.status_message
                },
                **span.attributes
            }
        )
    
    def get_span(self, span_id: str) -> Optional[Span]:
        """
        Get a span by ID.
        
        Args:
            span_id: Span ID
            
        Returns:
            Span or None if not found
        """
        with self.lock:
            # Check active spans first
            if span_id in self.active_spans:
                return self.active_spans[span_id]
            
            # Check completed spans
            return self.spans.get(span_id)
    
    def get_active_spans(self) -> List[Span]:
        """
        Get all active spans.
        
        Returns:
            List of active spans
        """
        with self.lock:
            return list(self.active_spans.values())
    
    def get_trace(self, trace_id: str) -> List[Span]:
        """
        Get all spans for a trace.
        
        Args:
            trace_id: Trace ID
            
        Returns:
            List of spans in the trace
        """
        with self.lock:
            # Combine active and completed spans
            all_spans = list(self.active_spans.values()) + list(self.spans.values())
            
            # Filter by trace ID
            return [span for span in all_spans if span.trace_id == trace_id]
    
    def clear_completed_spans(self) -> None:
        """Clear all completed spans."""
        with self.lock:
            self.spans.clear()
    
    def extract_context(self, headers: Dict[str, str]) -> Dict[str, str]:
        """
        Extract trace context from headers.
        
        Args:
            headers: HTTP headers
            
        Returns:
            Dictionary with trace context
        """
        # Standard W3C Trace Context headers
        trace_id = headers.get('traceparent', '').split('-')[1] if 'traceparent' in headers else None
        span_id = headers.get('traceparent', '').split('-')[2] if 'traceparent' in headers else None
        
        # Fallback to custom headers
        trace_id = trace_id or headers.get('X-Trace-ID')
        span_id = span_id or headers.get('X-Parent-Span-ID')
        
        # Return context
        return {
            'trace_id': trace_id,
            'parent_span_id': span_id
        }
    
    def inject_context(self, span: Span) -> Dict[str, str]:
        """
        Inject trace context into headers.
        
        Args:
            span: Current span
            
        Returns:
            Dictionary with trace context headers
        """
        # Standard W3C Trace Context headers
        traceparent = f"00-{span.trace_id}-{span.span_id}-01"
        
        # Custom headers
        headers = {
            'traceparent': traceparent,
            'X-Trace-ID': span.trace_id,
            'X-Span-ID': span.span_id,
            'X-Parent-Span-ID': span.parent_span_id or ''
        }
        
        return headers


# Global tracer instance
tracer = Tracer.get_instance()


def get_current_span() -> Optional[Span]:
    """
    Get the current active span from the execution context.
    
    Returns:
        Current span or None if no active span
    """
    return current_span_var.get()


def extract_trace_context(headers: Dict[str, str]) -> Dict[str, str]:
    """
    Extract trace context from headers.
    
    Args:
        headers: HTTP headers
        
    Returns:
        Dictionary with trace context
    """
    return tracer.extract_context(headers)


def inject_trace_context(span: Optional[Span] = None) -> Dict[str, str]:
    """
    Inject trace context into headers.
    
    Args:
        span: Span to inject (uses current span if None)
        
    Returns:
        Dictionary with trace context headers
    """
    span = span or get_current_span()
    if span:
        return tracer.inject_context(span)
    return {}


class SpanContext:
    """
    Context manager for span creation and management.
    
    Automatically starts a span at context entry and ends it at context exit.
    The span is set as the current span during the context.
    """
    
    def __init__(
        self,
        name: str,
        trace_id: Optional[str] = None,
        parent_span_id: Optional[str] = None,
        span_kind: SpanKind = SpanKind.INTERNAL,
        attributes: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize the span context.
        
        Args:
            name: Name of the operation
            trace_id: Trace ID (generated if None)
            parent_span_id: Parent span ID (None for root span)
            span_kind: Type of span
            attributes: Additional attributes for the span
        """
        self.name = name
        self.trace_id = trace_id
        self.parent_span_id = parent_span_id
        self.span_kind = span_kind
        self.attributes = attributes
        self.span = None
        self.previous_span = None
    
    def __enter__(self) -> Span:
        """
        Enter the span context.
        
        Returns:
            Created span
        """
        # Start a new span
        self.span = tracer.start_span(
            name=self.name,
            trace_id=self.trace_id,
            parent_span_id=self.parent_span_id,
            span_kind=self.span_kind,
            attributes=self.attributes
        )
        
        # Make it the current span and save the previous span
        self.previous_span = self.span.make_current()
        
        return self.span
    
    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """
        Exit the span context.
        
        Args:
            exc_type: Exception type
            exc_val: Exception value
            exc_tb: Exception traceback
        """
        if not self.span:
            return
        
        # Set status based on exception
        if exc_type:
            self.span.set_status('error', str(exc_val))
            # Add exception details as an event
            self.span.add_event(
                name='exception',
                attributes={
                    'exception.type': exc_type.__name__,
                    'exception.message': str(exc_val),
                    'exception.stacktrace': ''.join(traceback.format_exception(exc_type, exc_val, exc_tb))
                }
            )
        else:
            self.span.set_status('ok')
        
        # Restore the previous span as current
        self.span.restore_previous(self.previous_span)
        
        # End the span
        tracer.end_span(self.span)


def with_span(
    name: Optional[str] = None,
    span_kind: SpanKind = SpanKind.INTERNAL,
    attributes: Optional[Dict[str, Any]] = None
) -> Callable[[Callable[..., T]], Callable[..., T]]:
    """
    Decorator to trace a function with a span.
    
    Args:
        name: Name of the operation (defaults to function name)
        span_kind: Type of span
        attributes: Additional attributes for the span
        
    Returns:
        Decorator function
    """
    def decorator(func: Callable[..., T]) -> Callable[..., T]:
        nonlocal name
        
        # Use function name if name not provided
        if name is None:
            name = func.__name__
        
        # Handle regular functions
        if not asyncio.iscoroutinefunction(func):
            @functools.wraps(func)
            def wrapper(*args, **kwargs) -> T:
                # Create function attributes
                func_attributes = {
                    'function': func.__name__,
                    'module': func.__module__
                }
                
                if attributes:
                    func_attributes.update(attributes)
                
                # Use span context manager
                with SpanContext(name, span_kind=span_kind, attributes=func_attributes) as span:
                    # Add arguments as attributes (be careful with sensitive data)
                    if args:
                        span.set_attribute('args_count', len(args))
                    if kwargs:
                        span.set_attribute('kwargs_keys', list(kwargs.keys()))
                    
                    # Execute the function
                    return func(*args, **kwargs)
            
            return wrapper
        
        # Handle async functions
        @functools.wraps(func)
        async def async_wrapper(*args, **kwargs) -> T:
            # Create function attributes
            func_attributes = {
                'function': func.__name__,
                'module': func.__module__,
                'is_async': True
            }
            
            if attributes:
                func_attributes.update(attributes)
            
            # Start a span
            span = tracer.start_span(
                name=name,
                span_kind=span_kind,
                attributes=func_attributes
            )
            
            # Add arguments as attributes (be careful with sensitive data)
            if args:
                span.set_attribute('args_count', len(args))
            if kwargs:
                span.set_attribute('kwargs_keys', list(kwargs.keys()))
            
            # Save the current span and set the new span as current
            previous_span = span.make_current()
            
            try:
                # Execute the async function
                result = await func(*args, **kwargs)
                
                # Set success status
                span.set_status('ok')
                
                return result
            except Exception as e:
                # Set error status and add exception details
                span.set_status('error', str(e))
                span.add_event(
                    name='exception',
                    attributes={
                        'exception.type': type(e).__name__,
                        'exception.message': str(e),
                        'exception.stacktrace': ''.join(traceback.format_exception(
                            type(e), e, e.__traceback__
                        ))
                    }
                )
                
                # Re-raise the exception
                raise
            finally:
                # Restore the previous span as current
                span.restore_previous(previous_span)
                
                # End the span
                tracer.end_span(span)
        
        return async_wrapper
    
    return decorator


# Alias for SpanContext for use as a context manager
def trace_function(
    name: str,
    trace_id: Optional[str] = None,
    parent_span_id: Optional[str] = None,
    span_kind: SpanKind = SpanKind.INTERNAL,
    attributes: Optional[Dict[str, Any]] = None
) -> ContextManager[Span]:
    """
    Context manager for tracing a function or block of code.
    
    Args:
        name: Name of the operation
        trace_id: Trace ID (generated if None)
        parent_span_id: Parent span ID (None for root span)
        span_kind: Type of span
        attributes: Additional attributes for the span
        
    Returns:
        Context manager that yields a span
    """
    return SpanContext(
        name=name,
        trace_id=trace_id,
        parent_span_id=parent_span_id,
        span_kind=span_kind,
        attributes=attributes
    )