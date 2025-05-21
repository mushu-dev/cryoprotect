#!/usr/bin/env python3
"""
Cache invalidation processor script.

This script is responsible for processing cache invalidation events in the database
and keeping the cache synchronized with the database. It's intended to be run
periodically via a scheduler or as a background process.
"""

import os
import time
import signal
import logging
import argparse
import sys
from typing import Dict, List, Any, Optional

# Import cache invalidation module
from database.cache_invalidation import (
    setup_cache_invalidation,
    process_invalidation_events,
    invalidation_processor
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('cache_processor.log')
    ]
)
logger = logging.getLogger(__name__)

# Global flag to control the main loop
running = True

def signal_handler(sig, frame):
    """Handle interrupt signals to gracefully stop the processor."""
    global running
    logger.info(f"Received signal {sig}, shutting down...")
    running = False

def setup_signals():
    """Setup signal handlers for graceful shutdown."""
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Cache invalidation processor')
    parser.add_argument('--daemon', '-d', action='store_true', help='Run as a daemon process')
    parser.add_argument('--interval', '-i', type=int, default=10, help='Polling interval in seconds')
    parser.add_argument('--setup', '-s', action='store_true', help='Setup cache invalidation infrastructure')
    parser.add_argument('--cleanup', '-c', action='store_true', help='Cleanup old invalidation events')
    parser.add_argument('--cleanup-days', type=int, default=7, help='Delete events older than this many days')
    return parser.parse_args()

def run_daemon(interval: int):
    """
    Run the processor as a daemon, continuously processing events.
    
    Args:
        interval: Polling interval in seconds
    """
    logger.info(f"Starting cache invalidation processor (interval: {interval}s)")
    setup_signals()
    
    try:
        # Initial setup
        setup_cache_invalidation()
        
        # Main loop
        while running:
            start_time = time.time()
            
            # Process invalidation events
            processed = process_invalidation_events()
            
            # Periodically clean up old events (every hour)
            if processed > 0 and time.time() % 3600 < interval:
                invalidation_processor.cleanup_old_events()
            
            # Sleep for the interval, accounting for processing time
            elapsed = time.time() - start_time
            sleep_time = max(0.1, interval - elapsed)
            
            if processed > 0:
                logger.info(f"Processed {processed} events in {elapsed:.2f}s, sleeping for {sleep_time:.2f}s")
            
            # Use a small sleep interval to allow faster response to signals
            for _ in range(int(sleep_time * 10)):
                if running:
                    time.sleep(0.1)
                else:
                    break
                    
    except Exception as e:
        logger.error(f"Error in daemon mode: {e}")
        return 1
        
    logger.info("Cache invalidation processor stopped")
    return 0

def main():
    """Main entry point."""
    args = parse_args()
    
    try:
        # Just setup cache invalidation infrastructure
        if args.setup:
            logger.info("Setting up cache invalidation infrastructure...")
            if setup_cache_invalidation():
                logger.info("Cache invalidation infrastructure setup successfully")
                return 0
            else:
                logger.error("Failed to setup cache invalidation infrastructure")
                return 1
        
        # Just cleanup old events
        if args.cleanup:
            logger.info(f"Cleaning up cache invalidation events older than {args.cleanup_days} days...")
            count = invalidation_processor.cleanup_old_events(days=args.cleanup_days)
            logger.info(f"Cleaned up {count} old events")
            return 0
        
        # Run as a daemon
        if args.daemon:
            return run_daemon(args.interval)
        
        # Run once
        logger.info("Processing cache invalidation events...")
        processed = process_invalidation_events()
        logger.info(f"Processed {processed} events")
        return 0
        
    except Exception as e:
        logger.error(f"Error: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())