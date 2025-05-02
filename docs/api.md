# API Documentation

This document describes the API endpoints and usage.

## OpenAPI Specification

[To be added]

## Usage Examples

### Python

#### Get a paginated list of molecules

```python
import requests

url = "https://api.cryoprotect.com/api/molecules?limit=10&offset=0"
headers = {"Authorization": "Bearer <token>"}

response = requests.get(url, headers=headers)

if response.status_code == 200:
    molecules = response.json()
    print(molecules)
else:
    print(f"Error: {response.status_code} - {response.text}")
```

#### Import a molecule from PubChem

```python
import requests
import json

url = "https://api.cryoprotect.com/api/molecules"
headers = {"Authorization": "Bearer <token>", "Content-Type": "application/json"}
data = {"cid": 2244}

response = requests.post(url, headers=headers, data=json.dumps(data))

if response.status_code == 201:
    molecule = response.json()
    print(molecule)
else:
    print(f"Error: {response.status_code} - {response.text}")
```

### JavaScript

#### Get a paginated list of molecules

```javascript
const url = "https://api.cryoprotect.com/api/molecules?limit=10&offset=0";
const headers = {"Authorization": "Bearer <token>"};

fetch(url, {headers: headers})
  .then(response => {
    if (response.ok) {
      return response.json();
    } else {
      throw new Error(`Error: ${response.status} - ${response.statusText}`);
    }
  })
  .then(molecules => console.log(molecules))
  .catch(error => console.error(error));
```

#### Import a molecule from PubChem

```javascript
const url = "https://api.cryoprotect.com/api/molecules";
const headers = {"Authorization": "Bearer <token>", "Content-Type": "application/json"};
const data = {"cid": 2244};

fetch(url, {
  method: "POST",
  headers: headers,
  body: JSON.stringify(data)
})
  .then(response => {
    if (response.status === 201) {
      return response.json();
    } else {
      throw new Error(`Error: ${response.status} - ${response.statusText}`);
    }
  })
  .then(molecule => console.log(molecule))
  .catch(error => console.error(error));
```

## Authentication

The API uses JWT (JSON Web Token) authentication. To access protected endpoints, you need to include a valid JWT token in the `Authorization` header.

To obtain a JWT token, you need to authenticate with the API using your credentials. The specific authentication process may vary depending on the API implementation. Please refer to the authentication documentation for details on how to obtain a JWT token.

Once you have a JWT token, you can include it in the `Authorization` header of your requests as follows:

```
Authorization: Bearer <token>
```

## Rate Limiting

The API uses rate limiting to protect against abuse and ensure the stability of the service. The default rate limits are:

*   200 requests per day
*   50 requests per hour
*   10 requests per minute

If you exceed the rate limit, you will receive a `429 Too Many Requests` error. The response will include a `Retry-After` header indicating how long to wait before making another request.

Example error response:

```json
{
    "status": "error",
    "message": "Rate limit exceeded",
    "details": "10 requests per minute",
    "retry_after": 60
}
```

To avoid exceeding the rate limit, you should:

*   Implement exponential backoff with jitter in your client code.
*   Cache API responses whenever possible.
*   Contact the API administrators if you need higher rate limits.

## Versioning

The API uses semantic versioning. The current API version is `1.0.0`.

The API version is included in the base URL of all API endpoints. For example:

```
https://api.cryoprotect.com/api/v1/molecules
```

In this example, `v1` indicates that the API version is 1.0.0.

To use a specific API version, you should include the version number in the base URL of your requests.