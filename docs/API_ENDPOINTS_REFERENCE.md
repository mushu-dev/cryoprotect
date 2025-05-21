# CryoProtect API Endpoints Reference

## Standard Endpoints

### Molecule Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/v1/molecules` | GET | List all molecules with optional filtering and pagination |
| `/api/v1/molecules/{molecule_id}` | GET | Get a specific molecule by ID |
| `/api/v1/molecules/{molecule_id}/calculate-properties` | GET | Calculate properties for a molecule |
| `/api/v1/molecules/{molecule_id}/score` | GET | Get cryoprotection score for a molecule |

### Mixture Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/v1/mixtures` | GET | List all mixtures with optional filtering and pagination |
| `/api/v1/mixtures` | POST | Create a new mixture |
| `/api/v1/mixtures/{mixture_id}` | GET | Get a specific mixture by ID |
| `/api/v1/mixtures/{mixture_id}` | PUT | Update a mixture |
| `/api/v1/mixtures/{mixture_id}` | DELETE | Delete a mixture |
| `/api/v1/mixtures/{mixture_id}/predictions` | GET | List predictions for a mixture |
| `/api/v1/mixtures/{mixture_id}/experiments` | GET | List experiments for a mixture |
| `/api/v1/mixtures/{mixture_id}/compare` | GET | Compare mixture with other mixtures |
| `/api/v1/mixtures/{mixture_id}/score` | GET | Get cryoprotection score for a mixture |

### Batch Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/v1/batch` | POST | Perform operations on multiple molecules or mixtures |

### RDKit Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/v1/rdkit/properties` | POST | Calculate molecular properties using RDKit |
| `/api/v1/rdkit/visualization` | POST | Generate molecular visualizations |
| `/api/v1/rdkit/substructure` | POST | Perform substructure search |
| `/api/v1/rdkit/similarity` | POST | Perform similarity search |

### Scoring Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/v1/scoring/molecules` | POST | Calculate scores for multiple molecules |

### Comparison Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/v1/compare-properties` | POST | Compare properties between molecules or mixtures |

### User Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/v1/user_profile` | GET | Get the current user's profile |
| `/api/v1/user_profile` | PUT | Update the current user's profile |

## Consolidated Molecule Endpoints

### Consolidated Molecule Access

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/v1/consolidated/molecules/{molecule_id}` | GET | Get a molecule with consolidated handling |
| `/api/v1/consolidated/batch` | POST | Batch operations with consolidated handling |
| `/api/v1/molecules/{molecule_id}/primary` | GET | Get the primary molecule for a molecule |
| `/api/v1/consolidated` | GET | List all consolidated molecule relationships |

### Differentiation Group Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/v1/differentiation/groups` | GET | List all differentiation groups |
| `/api/v1/differentiation/groups/{group_id}` | GET | Get details about a differentiation group |
| `/api/v1/molecules/{molecule_id}/differentiation` | GET | Get differentiation information for a molecule |

## Standardized Response Format

All API endpoints return responses in this standardized format:

```json
{
  "status": "success",
  "timestamp": "2023-05-13T12:34:56.789Z",
  "code": 200,
  "message": "Request succeeded",
  "data": {
    // Response data
  }
}
```

Error responses follow a similar format:

```json
{
  "status": "error",
  "timestamp": "2023-05-13T12:34:56.789Z",
  "code": 400,
  "message": "Bad request",
  "errors": [
    {
      "type": "ValidationError",
      "message": "Invalid parameter value",
      "context": "Request validation failed"
    }
  ]
}
```

## Pagination

Endpoints that return lists of items support pagination with these query parameters:

- `page`: Page number (default: 1)
- `per_page`: Number of items per page (default: 20, max: 100)

Example: `/api/v1/molecules?page=2&per_page=50`

Paginated responses include pagination information:

```json
{
  "status": "success",
  "timestamp": "2023-05-13T12:34:56.789Z",
  "code": 200,
  "message": "Request succeeded",
  "data": [
    // Items
  ],
  "pagination": {
    "page": 2,
    "per_page": 50,
    "total_items": 150,
    "total_pages": 3,
    "has_next": false,
    "has_prev": true,
    "links": {
      "self": "/api/v1/molecules?page=2&per_page=50",
      "first": "/api/v1/molecules?page=1&per_page=50",
      "last": "/api/v1/molecules?page=3&per_page=50",
      "prev": "/api/v1/molecules?page=1&per_page=50"
    }
  }
}
```

## Authentication

Protected endpoints require authentication via JWT tokens:

- Include token in the `Authorization` header: `Bearer <token>`
- Token can be obtained from the authentication endpoints
- Tokens expire after a configurable period (default: 1 hour)

## Rate Limiting

API endpoints are rate-limited to prevent abuse:

- Default limit: 120 requests per minute per IP address
- Batch operations count as a single request
- Rate limit headers are included in responses:
  - `X-RateLimit-Limit`: Maximum number of requests allowed in the period
  - `X-RateLimit-Remaining`: Number of requests remaining in the period
  - `X-RateLimit-Reset`: Timestamp when the rate limit resets

## Error Codes

| Code | Description |
|------|-------------|
| 200 | OK - Request succeeded |
| 201 | Created - Resource created successfully |
| 204 | No Content - Request succeeded with no content to return |
| 400 | Bad Request - Invalid request parameters or data |
| 401 | Unauthorized - Authentication required |
| 403 | Forbidden - Permission denied |
| 404 | Not Found - Resource not found |
| 409 | Conflict - Request conflicts with current state of the resource |
| 422 | Unprocessable Entity - Request data validation failed |
| 429 | Too Many Requests - Rate limit exceeded |
| 500 | Internal Server Error - Server encountered an error |
| 503 | Service Unavailable - Service temporarily unavailable |