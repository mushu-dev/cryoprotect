# Enhanced Logging System for CryoProtect v2

This document describes the enhanced, production-ready logging system implemented for CryoProtect v2 as part of PHASE_3.1_STAGE3_TASK_3.1.

## Features

- **JSON-structured logging format** - All logs are formatted as JSON for easy parsing and analysis
- **Correlation IDs** - Each request is assigned a unique correlation ID that is propagated through all logs
- **Contextual information** - Logs include rich contextual information such as user ID, request details, and system information
- **Log rotation and retention policies** - Logs are automatically rotated based on size and retained for a configurable period
- **ELK stack integration** - Logs can be shipped to Elasticsearch for centralized storage and analysis
- **Request/response logging middleware** - All HTTP requests and responses are automatically logged with detailed information

## Components

### 1. Enhanced Logging Module (`logging_enhanced.py`)

The core of the logging system is the `logging_enhanced.py` module, which provides:

- JSON formatting for logs
- Contextual filters to add request and system information
- Handlers for console, file, and Elasticsearch output
- Log rotation based on file size
- Helper functions for logging with context

### 2. ELK Stack Integration

The system includes Docker Compose configuration for the ELK stack:

- **Elasticsearch** - For storing and indexing logs
- **Logstash** - For processing and transforming logs
- **Kibana** - For visualizing and analyzing logs
- **Filebeat** - For shipping logs from files to Elasticsearch

### 3. Request/Response Logging Middleware

The system includes middleware for Flask that:

- Logs all incoming requests with details
- Logs all outgoing responses with status codes and timing information
- Adds correlation IDs to responses for client-side tracking

## Configuration

The logging system can be configured using environment variables:

| Variable | Description | Default |
|----------|-------------|---------|
| LOG_LEVEL | Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL) | INFO |
| LOG_FILE | Path to the log file | logs/cryoprotect.log |
| LOG_TO_FILE | Whether to log to a file (1=yes, 0=no) | 1 |
| LOG_TO_CONSOLE | Whether to log to the console (1=yes, 0=no) | 1 |
| LOG_TO_ELK | Whether to log to Elasticsearch (1=yes, 0=no) | 0 |
| LOG_JSON_FORMAT | Whether to use JSON format (1=yes, 0=no) | 1 |
| LOG_ROTATION_SIZE | Maximum size of log file before rotation (bytes) | 10485760 (10MB) |
| LOG_BACKUP_COUNT | Number of backup log files to keep | 10 |
| ELASTICSEARCH_HOST | Elasticsearch host and port | elasticsearch:9200 |
| ELASTICSEARCH_INDEX | Elasticsearch index name | cryoprotect-logs |
| LOG_CORRELATION_ID_HEADER | HTTP header for correlation ID | X-Correlation-ID |
| LOG_RETENTION_DAYS | Number of days to retain logs | 30 |
| LOG_INCLUDE_CONTEXT | Whether to include context in logs (1=yes, 0=no) | 1 |

## Usage

### Basic Usage

The enhanced logging system is automatically set up when the application starts. No additional configuration is needed for basic usage.

### Logging with Context

To log with additional context:

```python
from logging_enhanced import get_logger, log_with_context

# Get a logger
logger = get_logger(__name__)

# Log with context
log_with_context(
    logger, 'info', 
    "User performed an action", 
    context={
        "user_id": "123",
        "action": "create_experiment",
        "experiment_id": "456"
    }
)
```

### Accessing Logs

Logs are available in multiple locations:

1. **Console** - When `LOG_TO_CONSOLE=1`, logs are output to the console
2. **Log Files** - When `LOG_TO_FILE=1`, logs are written to the configured log file
3. **Kibana** - When using the ELK stack, logs can be viewed and analyzed in Kibana at http://localhost:5601

## ELK Stack Setup

To start the ELK stack:

```bash
# Start the ELK stack
docker-compose -f docker-compose.yml -f docker-compose.elk.yml up -d

# View logs in Kibana
# Open http://localhost:5601 in your browser
```

### Kibana Dashboards

The system includes pre-configured Kibana dashboards for:

- API request monitoring
- Error tracking
- Performance analysis
- User activity monitoring

To import these dashboards, go to Kibana > Management > Saved Objects > Import and select the dashboard JSON files from the `elk/kibana/dashboards` directory.

## Testing

Unit tests for the enhanced logging system are available in `test_logging_enhanced.py`. To run the tests:

```bash
python -m unittest test_logging_enhanced.py
```

## Best Practices

1. **Use correlation IDs** - Always include correlation IDs in logs to track requests across services
2. **Add context** - Include relevant context in logs to make them more useful for debugging
3. **Use structured logging** - Use the `log_with_context` function to add structured data to logs
4. **Monitor Kibana dashboards** - Regularly check Kibana dashboards for errors and performance issues
5. **Set appropriate log levels** - Use appropriate log levels (DEBUG, INFO, WARNING, ERROR) for different types of messages