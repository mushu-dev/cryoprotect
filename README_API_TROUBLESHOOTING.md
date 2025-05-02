# CryoProtect v2 API Troubleshooting Guide

This document provides a comprehensive guide for troubleshooting the CryoProtect v2 API and fixing common issues.

## Table of Contents

1. [Overview](#overview)
2. [Improved Tools](#improved-tools)
3. [Common Issues](#common-issues)
4. [Troubleshooting Workflow](#troubleshooting-workflow)
5. [Monitoring and Diagnostics](#monitoring-and-diagnostics)
6. [Performance Optimization](#performance-optimization)
7. [Recovery Procedures](#recovery-procedures)

## Overview

The CryoProtect v2 API consists of 16 endpoints that provide access to molecule and mixture data, predictions, experiments, and various analysis tools. The API interacts with a Supabase database and uses RDKit for molecular calculations.

Recent fixes have addressed several issues:
- Database table name mismatches (singular vs. plural)
- Endpoint registration issues
- JSON serialization problems
- Error handling improvements

## Improved Tools

We've created several improved tools to help troubleshoot and fix API issues:

### 1. `run_api_fixes_improved.py`

An enhanced version of the original fix script with:
- Real-time progress monitoring
- Timeout handling to prevent hanging
- Detailed performance logging
- Checkpoint system for recovery
- System diagnostics
- Improved error handling

Usage:
```bash
python run_api_fixes_improved.py
```

### 2. `verify_api_standalone.py`

A standalone script to verify API endpoints independently:
- Can test individual endpoints
- Provides detailed error information
- Generates comprehensive reports
- Can be run without applying fixes

Usage:
```bash
# Test all endpoints
python verify_api_standalone.py

# Test a specific endpoint
python verify_api_standalone.py --endpoint /api/v1/mixtures
```

### 3. `fix_database_modular.py`

A modular script to fix database issues:
- Check table existence
- Create missing tables
- Copy data between tables
- Verify table structure

Usage:
```bash
# Run all fixes
python fix_database_modular.py

# Check if required tables exist
python fix_database_modular.py --action check_tables

# Create missing tables
python fix_database_modular.py --action create_tables

# Copy data from singular to plural tables
python fix_database_modular.py --action copy_data

# Verify table structure
python fix_database_modular.py --action verify_structure
```

### 4. `api_monitoring_dashboard.py`

A real-time monitoring dashboard for the API:
- Monitors endpoint health
- Tracks database connection status
- Displays performance metrics
- Logs errors and warnings

Usage:
```bash
python api_monitoring_dashboard.py --port 5001
```

## Common Issues

### 1. Database Configuration Issues

**Symptoms:**
- 500 errors from database-related endpoints
- Missing data in responses
- Database connection errors

**Solutions:**
- Run `fix_database_modular.py` to check and fix database tables
- Verify environment variables in `.env` file
- Check Supabase connection and permissions

### 2. Endpoint Registration Issues

**Symptoms:**
- 404 errors when accessing certain endpoints
- Endpoints not appearing in verification reports

**Solutions:**
- Check endpoint registration in `api/__init__.py`
- Verify URL patterns match the expected format
- Ensure resources are properly imported and registered

### 3. JSON Serialization Issues

**Symptoms:**
- 500 errors when returning data
- Malformed JSON responses
- TypeError or AttributeError in logs

**Solutions:**
- Check response handling in resource classes
- Ensure proper data type conversion
- Add error handling for JSON serialization

### 4. Performance Issues

**Symptoms:**
- Slow response times
- Timeouts
- High CPU/memory usage

**Solutions:**
- Use the monitoring dashboard to identify bottlenecks
- Check database query performance
- Optimize resource-intensive operations

## Troubleshooting Workflow

Follow this workflow to diagnose and fix API issues:

1. **Identify the Issue**
   - Run the monitoring dashboard to get an overview of system health
   - Check logs for error messages
   - Verify which endpoints are failing

2. **Isolate the Problem**
   - Use `verify_api_standalone.py` to test specific endpoints
   - Check database status with `fix_database_modular.py --action check_tables`
   - Review code for the affected endpoints

3. **Apply Fixes**
   - For database issues: Use `fix_database_modular.py`
   - For code issues: Update the relevant files
   - For configuration issues: Check environment variables and settings

4. **Verify Fixes**
   - Run `verify_api_standalone.py` to confirm endpoints are working
   - Use the monitoring dashboard to check system health
   - Test the application to ensure end-to-end functionality

## Monitoring and Diagnostics

### Using the Monitoring Dashboard

The monitoring dashboard provides real-time information about:
- API endpoint status
- Database connection health
- Performance metrics
- Error logs

To interpret the dashboard:
- **Green** status indicates normal operation
- **Yellow** status indicates warnings or potential issues
- **Red** status indicates errors that need immediate attention

### Log Analysis

Important log files:
- `api_fixes_improved.log`: Logs from the improved fix script
- `fix_database_modular.log`: Logs from the database fix script
- `api_verification_standalone.log`: Logs from the verification script
- `api_monitoring.log`: Logs from the monitoring dashboard
- `performance.log`: Performance metrics and timing information

## Performance Optimization

To optimize API performance:

1. **Database Optimization**
   - Ensure proper indexes are in place
   - Optimize queries for frequently accessed data
   - Consider caching for static data

2. **Code Optimization**
   - Minimize database calls
   - Use efficient algorithms for data processing
   - Implement proper error handling to avoid cascading failures

3. **System Resources**
   - Monitor CPU and memory usage
   - Scale resources as needed
   - Consider distributed processing for intensive operations

## Recovery Procedures

If the API becomes unresponsive or encounters critical errors:

1. **Emergency Restart**
   - Stop the running API process
   - Clear any temporary files or locks
   - Restart the API with proper error logging

2. **Database Recovery**
   - Run `fix_database_modular.py` to check and fix database issues
   - Verify data integrity
   - Restore from backup if necessary

3. **Checkpoint Recovery**
   - Check the `checkpoints` directory for the last successful operation
   - Resume from the last checkpoint
   - Skip completed steps to save time

4. **Manual Intervention**
   - If automated recovery fails, manual intervention may be required
   - Check logs for specific error messages
   - Fix issues directly in the code or database

## Conclusion

This troubleshooting guide provides a comprehensive approach to diagnosing and fixing issues with the CryoProtect v2 API. By using the improved tools and following the recommended workflows, you can efficiently identify and resolve problems to maintain a stable and performant API.

For additional assistance, refer to the other documentation files in the project, such as `README_API_FIXES.md` and `README_API.md`.