# Memory Bank Directory

This directory is used by RooFlow to store and manage memory data for the CryoProtect Analyzer application. The Universal Memory Bank (UMB) system allows for persistent storage of important information across sessions.

## Purpose

- Store conversation history
- Cache API responses
- Save user preferences
- Maintain state between sessions
- Optimize API calls by reusing previously fetched data

## Structure

The memory-bank directory will contain JSON files organized by session ID and timestamp. These files should not be manually edited unless absolutely necessary.

## Usage

The memory-bank is automatically managed by the RooFlow system. The application will read from and write to this directory as needed during operation.

## Maintenance

Periodically, old or unused memory files may be archived or removed to prevent the directory from growing too large. This process is typically handled automatically but can be manually triggered if needed.