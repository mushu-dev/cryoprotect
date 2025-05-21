#!/usr/bin/env python3
"""
Cache monitoring utility.

This module provides functions for monitoring the cache status and performance metrics.
It tracks cache hit/miss rates, memory usage, and other relevant statistics.
"""

import os
import time
import logging
import json
from typing import Dict, List, Any, Optional
from datetime import datetime, timedelta

from database.cache import get_cache_stats, get_redis_client
from database.cache_invalidation import invalidation_processor

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Metrics storage
_metrics_history = []
_metrics_start_time = datetime.now()
_max_history_size = 1000  # Maximum number of metric points to store

def get_detailed_cache_metrics() -> Dict[str, Any]:
    """
    Get detailed cache metrics and statistics.
    
    Returns:
        Dictionary with cache metrics
    """
    # Get basic cache stats
    cache_stats = get_cache_stats()
    
    # Get Redis info if available
    try:
        client = get_redis_client()
        info = client.info()
        
        # Calculate hit rate percentage
        hits = info.get('keyspace_hits', 0)
        misses = info.get('keyspace_misses', 0)
        total = hits + misses
        hit_rate_pct = (hits / total * 100) if total > 0 else 0
        
        # Memory metrics
        memory_used = info.get('used_memory', 0)
        memory_peak = info.get('used_memory_peak', 0)
        memory_rss = info.get('used_memory_rss', 0)
        
        # Get invalidation metrics
        try:
            from database.db import execute_query
            invalidation_count = execute_query(
                "SELECT COUNT(*) as count FROM cache_invalidation_events"
            )[0]['count']
            
            pending_count = execute_query(
                "SELECT COUNT(*) as count FROM cache_invalidation_events WHERE processed_at IS NULL"
            )[0]['count']
        except Exception as e:
            logger.warning(f"Error getting invalidation metrics: {e}")
            invalidation_count = 0
            pending_count = 0
        
        # Key metrics by namespace
        namespace_stats = {}
        patterns = [
            'cryoprotect:molecule:*',
            'cryoprotect:property:*',
            'cryoprotect:query:*',
            'cryoprotect:calculation:*',
            'cryoprotect:score:*',
            'cryoprotect:mixture:*'
        ]
        
        for pattern in patterns:
            try:
                key_count = len(list(client.scan_iter(match=pattern)))
                namespace = pattern.split(':')[1] if len(pattern.split(':')) > 1 else pattern
                namespace_stats[namespace] = key_count
            except Exception as e:
                logger.warning(f"Error getting key count for pattern {pattern}: {e}")
        
        # Combine all metrics
        metrics = {
            **cache_stats,
            'hit_rate_pct': round(hit_rate_pct, 2),
            'memory': {
                'used_bytes': memory_used,
                'used_human': info.get('used_memory_human', 'unknown'),
                'peak_bytes': memory_peak,
                'peak_human': info.get('used_memory_peak_human', 'unknown'),
                'rss_bytes': memory_rss,
                'rss_human': info.get('used_memory_rss_human', 'unknown'),
            },
            'uptime': info.get('uptime_in_seconds', 0),
            'clients': info.get('connected_clients', 0),
            'invalidation': {
                'total_events': invalidation_count,
                'pending_events': pending_count,
            },
            'namespaces': namespace_stats,
            'timestamp': datetime.now().isoformat()
        }
        
        return metrics
    except Exception as e:
        logger.warning(f"Error getting detailed cache metrics: {e}")
        
        # Return basic stats if Redis is not available
        return {
            **cache_stats,
            'error': str(e),
            'timestamp': datetime.now().isoformat()
        }

def record_metrics() -> Dict[str, Any]:
    """
    Record current cache metrics to history.
    
    Returns:
        The recorded metrics
    """
    global _metrics_history
    
    # Get current metrics
    metrics = get_detailed_cache_metrics()
    
    # Add to history
    _metrics_history.append(metrics)
    
    # Trim history if needed
    if len(_metrics_history) > _max_history_size:
        _metrics_history = _metrics_history[-_max_history_size:]
    
    return metrics

def get_metrics_history(hours: int = 24) -> List[Dict[str, Any]]:
    """
    Get cache metrics history.
    
    Args:
        hours: Number of hours of history to return
        
    Returns:
        List of metrics snapshots
    """
    if not _metrics_history:
        return []
    
    # Filter by time
    cutoff_time = datetime.now() - timedelta(hours=hours)
    filtered_metrics = []
    
    for metrics in _metrics_history:
        try:
            timestamp = datetime.fromisoformat(metrics['timestamp'])
            if timestamp >= cutoff_time:
                filtered_metrics.append(metrics)
        except (KeyError, ValueError):
            # Skip metrics with missing or invalid timestamp
            pass
    
    return filtered_metrics

def generate_metrics_report() -> str:
    """
    Generate a human-readable metrics report.
    
    Returns:
        Text report of cache metrics
    """
    metrics = get_detailed_cache_metrics()
    
    if 'error' in metrics:
        return f"Cache Monitor Report (ERROR)\n\nError: {metrics['error']}\n"
    
    # Format the report
    report = [
        "Cache Monitor Report",
        f"Timestamp: {metrics['timestamp']}",
        "",
        f"Status: {'Enabled' if metrics['enabled'] else 'Disabled'}",
        f"Total Keys: {metrics['keys']}",
        f"Hit Rate: {metrics.get('hit_rate_pct', 0)}%",
        f"Memory Used: {metrics['memory']['used_human']}",
        f"Memory Peak: {metrics['memory']['peak_human']}",
        "",
        "Keys by Namespace:",
    ]
    
    for namespace, count in metrics.get('namespaces', {}).items():
        report.append(f"- {namespace}: {count}")
    
    report.extend([
        "",
        "Invalidation:",
        f"- Total Events: {metrics['invalidation']['total_events']}",
        f"- Pending Events: {metrics['invalidation']['pending_events']}",
        "",
        f"Redis Uptime: {timedelta(seconds=metrics['uptime'])}",
        f"Connected Clients: {metrics['clients']}",
        "",
        f"Monitor Runtime: {datetime.now() - _metrics_start_time}"
    ])
    
    return "\n".join(report)

def save_metrics_to_file(file_path: str = "cache_metrics.json") -> bool:
    """
    Save current metrics to a JSON file.
    
    Args:
        file_path: Path to the output file
        
    Returns:
        True if successful, False otherwise
    """
    try:
        metrics = get_detailed_cache_metrics()
        
        with open(file_path, 'w') as f:
            json.dump(metrics, f, indent=2)
            
        logger.info(f"Saved cache metrics to {file_path}")
        return True
    except Exception as e:
        logger.error(f"Error saving metrics to file: {e}")
        return False

def process_and_monitor() -> Dict[str, Any]:
    """
    Process invalidation events and monitor cache metrics.
    
    Returns:
        Dictionary with processing and monitoring results
    """
    # Process invalidation events
    from database.cache_invalidation import process_invalidation_events
    processed = process_invalidation_events()
    
    # Record metrics
    metrics = record_metrics()
    
    # Return combined results
    return {
        'processed_events': processed,
        'metrics': metrics,
        'timestamp': datetime.now().isoformat()
    }

if __name__ == "__main__":
    # If run directly, print the metrics report
    print(generate_metrics_report())
    
    # Save metrics to file
    save_metrics_to_file()