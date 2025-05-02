# CryoProtect v2 API Documentation

This directory contains comprehensive documentation for the CryoProtect v2 API.

## Documentation Overview

The CryoProtect v2 API documentation is provided in OpenAPI/Swagger format, which offers a standardized way to describe RESTful APIs. The documentation includes:

- **Complete API Reference**: All endpoints, methods, parameters, and responses
- **Authentication Details**: How to authenticate with the API
- **Error Handling**: Error codes and recovery procedures
- **Rate Limiting Information**: How rate limiting is implemented and how to handle it
- **Usage Examples**: Examples of how to use the API

## Documentation Features

### OpenAPI/Swagger Specification

The API is documented using the OpenAPI 3.0.2 specification, which provides:

- **Interactive Documentation**: The Swagger UI allows you to explore the API interactively
- **Machine-Readable Format**: The YAML format can be used by tools to generate client libraries
- **Validation**: Request and response validation based on the schema

### Authentication Documentation

The documentation includes detailed information about the authentication process:

- JWT (JSON Web Token) authentication
- Token acquisition and refresh
- Multi-factor authentication (MFA)
- Session management

### Error Handling

Comprehensive error handling documentation is provided:

- Standard HTTP status codes
- Detailed error messages and codes
- Recovery procedures for different error scenarios

### Rate Limiting

The documentation includes information about rate limiting:

- Rate limit headers
- Rate limit thresholds
- How to handle rate limit exceeded errors

## Using the Documentation

### Swagger UI

You can access the interactive Swagger UI by running the CryoProtect v2 application and navigating to:

```
http://localhost:5000/swagger-ui/
```

This interface allows you to:

1. Browse all available endpoints
2. See request and response schemas
3. Try out API calls directly from the browser
4. View authentication requirements

### OpenAPI YAML File

The `openapi.yaml` file in this directory contains the complete API specification. You can use this file with various tools:

- **Code Generation**: Generate client libraries in various languages
- **API Testing**: Use with tools like Postman or Insomnia
- **Documentation Generation**: Create custom documentation

## Key Endpoints

The API is organized into several categories:

### Authentication

- `POST /auth/login`: Authenticate a user and get an access token
- `POST /auth/logout`: Log out the current user
- `POST /auth/refresh`: Refresh an access token
- `GET /auth/validate`: Validate an access token
- `POST /auth/register`: Register a new user
- `POST /auth/reset-password`: Request a password reset
- `POST /auth/update-password`: Update a user's password
- `POST /auth/mfa/initiate`: Initiate multi-factor authentication
- `POST /auth/mfa/verify`: Verify a multi-factor authentication challenge

### Molecules

- `GET /api/v1/molecules`: Get a list of molecules
- `POST /api/v1/molecules`: Create a new molecule
- `GET /api/v1/molecules/{molecule_id}`: Get a specific molecule
- `PUT /api/v1/molecules/{molecule_id}`: Update a molecule
- `DELETE /api/v1/molecules/{molecule_id}`: Delete a molecule
- `POST /api/v1/molecules/{molecule_id}/calculate-properties`: Calculate properties for a molecule
- `GET /api/v1/molecules/{molecule_id}/score`: Get the score for a molecule

### Mixtures

- `GET /api/v1/mixtures`: Get a list of mixtures
- `POST /api/v1/mixtures`: Create a new mixture
- `GET /api/v1/mixtures/{mixture_id}`: Get a specific mixture
- `PUT /api/v1/mixtures/{mixture_id}`: Update a mixture
- `DELETE /api/v1/mixtures/{mixture_id}`: Delete a mixture
- `GET /api/v1/mixtures/{mixture_id}/score`: Get the score for a mixture
- `GET /api/v1/mixtures/{mixture_id}/predictions`: Get predictions for a mixture
- `POST /api/v1/mixtures/{mixture_id}/predictions`: Create a prediction for a mixture

### Comparisons and Batch Operations

- `POST /api/v1/compare-properties`: Compare properties of molecules and/or mixtures
- `POST /api/v1/batch`: Perform batch operations on multiple items

### RDKit Integration

- `POST /api/v1/rdkit/properties`: Calculate properties for a molecule using RDKit
- `POST /api/v1/rdkit/visualization`: Generate a visualization of a molecule
- `POST /api/v1/rdkit/substructure`: Search for molecules containing a specific substructure
- `POST /api/v1/rdkit/similarity`: Search for molecules similar to a reference molecule
- `POST /api/v1/rdkit-enhanced/batch-calculate`: Calculate properties for multiple molecules in batch

### Scoring

- `POST /api/v1/scoring/molecules`: Score one or more molecules

### Export

- `POST /api/v1/export`: Export data in various formats

### Dashboard and User Profile

- `GET /api/v1/dashboard`: Get dashboard data
- `GET /api/v1/user_profile`: Get the current user's profile
- `PUT /api/v1/user_profile`: Update the current user's profile

## Integration Examples

For detailed examples of how to integrate with the CryoProtect v2 API, see the [API Usage Examples](../developer/api-usage-examples.md) document, which includes code samples in JavaScript and Python.

## Security Considerations

When using the API, be aware of the following security considerations:

1. **Authentication**: Always use HTTPS for production environments
2. **Token Storage**: Store tokens securely and refresh them as needed
3. **Rate Limiting**: Implement proper backoff strategies for rate limit errors
4. **Error Handling**: Handle errors gracefully in your application

## Further Resources

- [Developer Guide](../developer-guide.md): General guide for developers
- [Technical Documentation](../technical-documentation.md): Technical details about the system
- [User Guide](../user-guide.md): Guide for end users