# CryoProtect v2 API Rate Limiting

This document describes the rate limiting implementation for the CryoProtect v2 API.

## Overview

Rate limiting is a technique used to control the amount of incoming and outgoing traffic to or from a network, application, or service. In the context of the CryoProtect v2 API, rate limiting helps to:

- Protect the API from abuse and DoS attacks
- Ensure fair usage of resources among all users
- Maintain system stability and performance
- Manage infrastructure costs

The CryoProtect v2 API implements a flexible rate limiting system that can be configured based on:

- User identity (authenticated users)
- IP address (unauthenticated users)
- Endpoint-specific limits
- User role-based limits

## Rate Limiting Strategy

The API uses a "fixed-window" rate limiting strategy by default, which counts requests in fixed time windows (e.g., per minute, per hour, per day). This strategy is simple to understand and implement, but can lead to "bursts" of traffic at the edges of time windows.

Alternative strategies that can be configured include:

- **Moving-window**: Tracks requests over a continuously moving time window, providing smoother rate limiting
- **Fixed-window-elastic-expiry**: Similar to fixed-window but with dynamic expiry times

## Rate Limit Headers

When rate limiting is enabled, the API includes the following headers in responses:

| Header | Description |
|--------|-------------|
| `X-RateLimit-Limit` | The maximum number of requests allowed in the current time window |
| `X-RateLimit-Remaining` | The number of requests remaining in the current time window |
| `X-RateLimit-Reset` | The time (in Unix timestamp) when the current rate limit window resets |
| `Retry-After` | The number of seconds to wait before making another request (only included when rate limit is exceeded or close to being exceeded) |
| `X-RateLimit-Policy` | Information about the rate limit policy applied to the request |
| `X-RateLimit-Role` | Role-specific rate limits applied to the request (if applicable) |

## Rate Limit Exceedance

When a rate limit is exceeded, the API responds with:

- HTTP status code `429 Too Many Requests`
- A JSON response body with details about the rate limit exceedance
- A `Retry-After` header indicating when the client should retry the request

Example response:

```json
{
  "status": "error",
  "message": "Rate limit exceeded",
  "details": "5 per minute",
  "retry_after": 30
}
```

## Configuration

Rate limiting is configured in the `rate_limit_config.py` file and can be customized using environment variables. The following settings can be adjusted:

### Global Settings

- `RATE_LIMIT_ENABLED`: Enable or disable rate limiting globally
- `RATE_LIMIT_DEFAULT`: Default rate limits applied to all endpoints
- `RATE_LIMIT_STORAGE_URL`: Storage backend for rate limiting data (memory or Redis)
- `RATE_LIMIT_STRATEGY`: Rate limiting strategy to use
- `RATE_LIMIT_HEADERS_ENABLED`: Enable or disable rate limit headers in responses
- `RATE_LIMIT_RETRY_AFTER`: Default retry-after value for rate limit exceeded responses

### Identification Strategy

- `RATE_LIMIT_BY`: How to identify clients for rate limiting:
  - `"ip"`: Limit by IP address only
  - `"user"`: Limit by user ID only (falls back to IP for unauthenticated requests)
  - `"hybrid"`: Limit by user ID if authenticated, otherwise by IP address (default)

### Endpoint-Specific Limits

The `RATE_LIMIT_ENDPOINTS` dictionary allows setting different rate limits for specific endpoints:

```python
RATE_LIMIT_ENDPOINTS = {
    "/api/v1/auth/login": ["20 per minute", "100 per hour"],
    "/api/v1/molecules": ["50 per minute", "200 per hour"],
    # ...
}
```

### Role-Based Limits

The `RATE_LIMIT_ROLES` dictionary allows setting different rate limits based on user roles:

```python
RATE_LIMIT_ROLES = {
    "admin": ["5000 per day", "500 per hour", "100 per minute"],
    "premium": ["2000 per day", "200 per hour", "40 per minute"],
    "basic": ["1000 per day", "100 per hour", "20 per minute"],
}
```

### Exempt Endpoints

The `RATE_LIMIT_EXEMPT` list allows specifying endpoints that should be exempt from rate limiting:

```python
RATE_LIMIT_EXEMPT = [
    "/api/v1/health",
    "/api/v1/docs",
    # ...
]
```

## Usage in Code

### Applying Rate Limits to Specific Endpoints

You can apply custom rate limits to specific endpoints using the `endpoint_rate_limit` decorator:

```python
from api.rate_limiter import endpoint_rate_limit

@app.route('/api/v1/resource')
@endpoint_rate_limit("5 per minute")
def get_resource():
    # ...
```

### Exempting Endpoints from Rate Limiting

You can exempt specific endpoints from rate limiting using the `exempt_from_rate_limit` decorator:

```python
from api.rate_limiter import exempt_from_rate_limit

@app.route('/api/v1/public-resource')
@exempt_from_rate_limit
def get_public_resource():
    # ...
```

## Environment Variables

Rate limiting can be configured using the following environment variables:

| Variable | Description | Default |
|----------|-------------|---------|
| `RATE_LIMIT_ENABLED` | Enable or disable rate limiting | `true` |
| `RATE_LIMIT_STORAGE_URL` | Storage backend URL (memory:// or redis://) | `memory://` |
| `RATE_LIMIT_STRATEGY` | Rate limiting strategy | `fixed-window` |
| `RATE_LIMIT_BY` | Identification strategy (ip, user, hybrid) | `hybrid` |
| `RATE_LIMIT_HEADERS_ENABLED` | Enable or disable rate limit headers | `true` |
| `RATE_LIMIT_RETRY_AFTER` | Default retry-after value in seconds | `60` |
| `RATE_LIMIT_ROLES` | JSON string of role-based limits | See `rate_limit_config.py` |

Example `.env` file configuration:

```
RATE_LIMIT_ENABLED=true
RATE_LIMIT_STORAGE_URL=redis://localhost:6379/0
RATE_LIMIT_STRATEGY=fixed-window
RATE_LIMIT_BY=hybrid
RATE_LIMIT_HEADERS_ENABLED=true
RATE_LIMIT_RETRY_AFTER=60
RATE_LIMIT_ROLES={"admin":["5000 per day"],"premium":["2000 per day"]}
```

## Distributed Rate Limiting

For production environments with multiple application instances, it's recommended to use Redis as the storage backend for rate limiting. This ensures that rate limits are enforced consistently across all instances.

To enable Redis storage, set the `RATE_LIMIT_STORAGE_URL` environment variable or configuration value:

```
RATE_LIMIT_STORAGE_URL=redis://localhost:6379/0
```

## Monitoring and Debugging

Rate limiting events are logged to the application logs. You can monitor these logs to identify patterns of rate limit exceedance and adjust the configuration as needed.

Log messages include:

- Rate limit configuration at application startup
- Rate limit exceedance events with details about the client and endpoint
- Errors in rate limit header generation

## Best Practices

1. **Start with conservative limits**: Begin with lower limits and increase them as needed based on actual usage patterns.
2. **Communicate limits to users**: Make sure your API documentation clearly communicates the rate limits to users.
3. **Provide backoff guidance**: Include the `Retry-After` header to help clients implement proper backoff strategies.
4. **Monitor rate limit events**: Regularly review logs to identify patterns of rate limit exceedance.
5. **Adjust limits based on user tiers**: Consider implementing different rate limits for different user tiers or subscription levels.

## Testing Rate Limiting

A test endpoint is available to verify that rate limiting is working correctly:

```
GET /api/v1/test/rate-limit
```

This endpoint is limited to 5 requests per minute. You can use it to test the rate limiting functionality:

```bash
# Make multiple requests to test rate limiting
curl -i http://localhost:5000/api/v1/test/rate-limit
```

After 5 requests within a minute, you should receive a 429 Too Many Requests response with appropriate headers and a JSON error message.