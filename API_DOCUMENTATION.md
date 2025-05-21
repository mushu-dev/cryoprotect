# CryoProtect API Documentation

This document provides information about the available API endpoints for frontend integration with the CryoProtect application.

## Base URL

The API is accessible at:
- Development: `/api/v1`
- Production: `https://api.cryoprotect.app/v1` (via Vercel)

## Authentication

Most endpoints require authentication using JWT tokens. Authentication can be handled in two ways:

1. **Authorization Header**: Include a Bearer token in the request header
   ```
   Authorization: Bearer <your_jwt_token>
   ```

2. **HTTP-Only Cookies**: For browser-based clients, authentication can be handled via HTTP-only cookies set during the login process.

## Core Resources

### 1. Molecules

#### 1.1 List Molecules

```
GET /molecules
```

**Query Parameters:**
- `limit` (optional): Maximum number of results to return (default: 20)
- `offset` (optional): Pagination offset (default: 0)
- `sort` (optional): Field to sort by (default: 'name')
- `order` (optional): Sort order ('asc' or 'desc', default: 'asc')
- `filter` (optional): JSON string with filter criteria

**Response:**
```json
{
  "status": "success",
  "data": {
    "molecules": [
      {
        "id": "uuid",
        "name": "string",
        "inchikey": "string",
        "smiles": "string",
        "molecular_formula": "string",
        "molecular_weight": "number",
        "is_consolidated": "boolean",
        "molecule_status": "string",
        "primary_molecule_id": "uuid",
        "created_at": "timestamp",
        "updated_at": "timestamp"
      }
    ],
    "count": "number",
    "total": "number"
  }
}
```

#### 1.2 Get Molecule

```
GET /molecules/{molecule_id}
```

**Response:**
```json
{
  "status": "success",
  "data": {
    "id": "uuid",
    "name": "string",
    "inchikey": "string",
    "smiles": "string",
    "molecular_formula": "string",
    "molecular_weight": "number",
    "is_consolidated": "boolean",
    "molecule_status": "string",
    "primary_molecule_id": "uuid",
    "created_at": "timestamp",
    "updated_at": "timestamp",
    "properties": [
      {
        "id": "uuid",
        "property_name": "string",
        "property_type": "string",
        "numeric_value": "number",
        "text_value": "string",
        "unit": "string",
        "source": "string"
      }
    ]
  }
}
```

#### 1.3 Create Molecule

```
POST /molecules
```

**Request Body:**
```json
{
  "name": "string",
  "inchikey": "string",
  "smiles": "string",
  "molecular_formula": "string",
  "molecular_weight": "number",
  "source": "string"
}
```

**Response:** Same as Get Molecule

#### 1.4 Update Molecule

```
PUT /molecules/{molecule_id}
```

**Request Body:**
```json
{
  "name": "string",
  "smiles": "string",
  "molecular_formula": "string",
  "molecular_weight": "number"
}
```

**Response:** Same as Get Molecule

#### 1.5 Delete Molecule

```
DELETE /molecules/{molecule_id}
```

**Response:**
```json
{
  "status": "success",
  "message": "Molecule deleted successfully"
}
```

### 2. Consolidated Molecules

#### 2.1 Get Consolidated Molecule

```
GET /consolidated-molecules/{molecule_id}
```

**Response:**
```json
{
  "status": "success",
  "data": {
    "id": "uuid",
    "name": "string",
    "inchikey": "string",
    "smiles": "string",
    "molecular_formula": "string",
    "molecular_weight": "number",
    "is_consolidated": "boolean",
    "molecule_status": "string",
    "primary_molecule_id": "uuid",
    "primary_molecule_name": "string",
    "primary_note": "string",
    "audit_history": [
      {
        "id": "uuid",
        "operation": "string",
        "timestamp": "timestamp",
        "user_id": "uuid",
        "old_value": "json",
        "new_value": "json"
      }
    ],
    "duplicate_molecules": [
      {
        "id": "uuid",
        "name": "string"
      }
    ],
    "duplicate_count": "number"
  }
}
```

#### 2.2 Find Potential Duplicate Molecules

```
GET /molecule-consolidation?inchikey={inchikey}
```

**Query Parameters:**
- `inchikey` (required): InChIKey to search for
- `limit` (optional): Maximum number of results to return (default: 100)

**Response:**
```json
{
  "status": "success",
  "data": {
    "molecules": [
      {
        "id": "uuid",
        "name": "string",
        "inchikey": "string",
        "smiles": "string",
        "molecular_formula": "string",
        "molecular_weight": "number",
        "is_consolidated": "boolean",
        "molecule_status": "string",
        "primary_molecule_id": "uuid"
      }
    ],
    "count": "number",
    "grouped_by_status": {
      "primary": [],
      "duplicate": [],
      "original": []
    },
    "inchikey": "string"
  }
}
```

#### 2.3 Update Consolidated Molecule Relationship

```
PUT /consolidated-molecules/{molecule_id}
```

**Request Body:**
```json
{
  "primary_molecule_id": "uuid"
}
```

**Response:**
```json
{
  "status": "success",
  "data": {
    "message": "Molecule successfully consolidated",
    "molecule": {
      "id": "uuid",
      "name": "string",
      "inchikey": "string",
      "smiles": "string",
      "molecular_formula": "string",
      "molecular_weight": "number",
      "is_consolidated": "boolean",
      "molecule_status": "string",
      "primary_molecule_id": "uuid",
      "primary_molecule_name": "string"
    },
    "audit_id": "uuid"
  }
}
```

#### 2.4 Batch Consolidate Molecules

```
POST /molecule-consolidation
```

**Request Body:**
```json
{
  "primary_molecule_id": "uuid",
  "duplicate_molecule_ids": ["uuid", "uuid", ...]
}
```

**Response:**
```json
{
  "status": "success",
  "data": {
    "message": "Successfully consolidated X molecules",
    "primary_molecule_id": "uuid",
    "primary_molecule_name": "string",
    "updated_molecule_ids": ["uuid", "uuid", ...],
    "failed_molecule_ids": ["uuid", "uuid", ...]
  }
}
```

#### 2.5 Migrate Molecule Properties

```
POST /molecule-property-migration
```

**Request Body:**
```json
{
  "source_molecule_id": "uuid",
  "target_molecule_id": "uuid",
  "property_ids": ["uuid", "uuid", ...] // Optional
}
```

**Response:**
```json
{
  "status": "success",
  "data": {
    "message": "Successfully migrated X properties",
    "source_molecule_id": "uuid",
    "target_molecule_id": "uuid",
    "migrated_properties": [
      {
        "source_property_id": "uuid",
        "target_property_id": "uuid",
        "property_name": "string"
      }
    ],
    "skipped_properties": ["uuid", "uuid", ...],
    "migrated_count": "number",
    "skipped_count": "number"
  }
}
```

### 3. Molecular Properties

#### 3.1 List Molecule Properties

```
GET /molecules/{molecule_id}/properties
```

**Response:**
```json
{
  "status": "success",
  "data": {
    "properties": [
      {
        "id": "uuid",
        "molecule_id": "uuid",
        "property_type_id": "uuid",
        "property_name": "string",
        "property_type": "string",
        "numeric_value": "number",
        "text_value": "string",
        "unit": "string",
        "source": "string",
        "created_at": "timestamp",
        "updated_at": "timestamp"
      }
    ],
    "count": "number"
  }
}
```

#### 3.2 Add Molecule Property

```
POST /molecules/{molecule_id}/properties
```

**Request Body:**
```json
{
  "property_type_id": "uuid",
  "property_name": "string",
  "property_type": "string",
  "numeric_value": "number",
  "text_value": "string",
  "unit": "string",
  "source": "string"
}
```

**Response:**
```json
{
  "status": "success",
  "data": {
    "id": "uuid",
    "molecule_id": "uuid",
    "property_type_id": "uuid",
    "property_name": "string",
    "property_type": "string",
    "numeric_value": "number",
    "text_value": "string",
    "unit": "string",
    "source": "string",
    "created_at": "timestamp",
    "updated_at": "timestamp"
  }
}
```

#### 3.3 Update Molecule Property

```
PUT /molecules/{molecule_id}/properties/{property_id}
```

**Request Body:**
```json
{
  "numeric_value": "number",
  "text_value": "string",
  "unit": "string",
  "source": "string"
}
```

**Response:** Same as Add Molecule Property

#### 3.4 Delete Molecule Property

```
DELETE /molecules/{molecule_id}/properties/{property_id}
```

**Response:**
```json
{
  "status": "success",
  "message": "Property deleted successfully"
}
```

### 4. Mixtures

#### 4.1 List Mixtures

```
GET /mixtures
```

**Query Parameters:**
- `limit` (optional): Maximum number of results to return (default: 20)
- `offset` (optional): Pagination offset (default: 0)
- `sort` (optional): Field to sort by (default: 'name')
- `order` (optional): Sort order ('asc' or 'desc', default: 'asc')

**Response:**
```json
{
  "status": "success",
  "data": {
    "mixtures": [
      {
        "id": "uuid",
        "name": "string",
        "description": "string",
        "created_at": "timestamp",
        "updated_at": "timestamp",
        "component_count": "number"
      }
    ],
    "count": "number",
    "total": "number"
  }
}
```

#### 4.2 Get Mixture

```
GET /mixtures/{mixture_id}
```

**Response:**
```json
{
  "status": "success",
  "data": {
    "id": "uuid",
    "name": "string",
    "description": "string",
    "created_at": "timestamp",
    "updated_at": "timestamp",
    "components": [
      {
        "id": "uuid",
        "mixture_id": "uuid",
        "molecule_id": "uuid",
        "molecule_name": "string",
        "concentration": "number",
        "concentration_unit": "string",
        "role": "string"
      }
    ],
    "properties": [
      {
        "id": "uuid",
        "property_name": "string",
        "property_type": "string",
        "numeric_value": "number",
        "text_value": "string",
        "unit": "string",
        "source": "string"
      }
    ]
  }
}
```

#### 4.3 Create Mixture

```
POST /mixtures
```

**Request Body:**
```json
{
  "name": "string",
  "description": "string",
  "components": [
    {
      "molecule_id": "uuid",
      "concentration": "number",
      "concentration_unit": "string",
      "role": "string"
    }
  ]
}
```

**Response:** Same as Get Mixture

#### 4.4 Update Mixture

```
PUT /mixtures/{mixture_id}
```

**Request Body:**
```json
{
  "name": "string",
  "description": "string"
}
```

**Response:** Same as Get Mixture

#### 4.5 Delete Mixture

```
DELETE /mixtures/{mixture_id}
```

**Response:**
```json
{
  "status": "success",
  "message": "Mixture deleted successfully"
}
```

#### 4.6 Add Mixture Component

```
POST /mixtures/{mixture_id}/components
```

**Request Body:**
```json
{
  "molecule_id": "uuid",
  "concentration": "number",
  "concentration_unit": "string",
  "role": "string"
}
```

**Response:**
```json
{
  "status": "success",
  "data": {
    "id": "uuid",
    "mixture_id": "uuid",
    "molecule_id": "uuid",
    "molecule_name": "string",
    "concentration": "number",
    "concentration_unit": "string",
    "role": "string"
  }
}
```

#### 4.7 Update Mixture Component

```
PUT /mixtures/{mixture_id}/components/{component_id}
```

**Request Body:**
```json
{
  "concentration": "number",
  "concentration_unit": "string",
  "role": "string"
}
```

**Response:** Same as Add Mixture Component

#### 4.8 Delete Mixture Component

```
DELETE /mixtures/{mixture_id}/components/{component_id}
```

**Response:**
```json
{
  "status": "success",
  "message": "Component deleted successfully"
}
```

### 5. RDKit Integration

#### 5.1 Get Molecule Properties from RDKit

```
GET /rdkit/properties?smiles={smiles}
```

**Query Parameters:**
- `smiles` (required): SMILES notation of the molecule

**Response:**
```json
{
  "status": "success",
  "data": {
    "properties": {
      "molecular_weight": "number",
      "logp": "number",
      "num_atoms": "number",
      "num_bonds": "number",
      "num_rings": "number",
      "tpsa": "number",
      "h_bond_donors": "number",
      "h_bond_acceptors": "number",
      "rotatable_bonds": "number"
    }
  }
}
```

#### 5.2 Get Molecule Visualization

```
GET /rdkit/visualization?smiles={smiles}&format={format}
```

**Query Parameters:**
- `smiles` (required): SMILES notation of the molecule
- `format` (optional): Output format ('svg', 'png', base64', default: 'svg')
- `width` (optional): Image width in pixels (default: 400)
- `height` (optional): Image height in pixels (default: 300)

**Response:**
```json
{
  "status": "success",
  "data": {
    "visualization": "string", // SVG, base64-encoded PNG, or URL
    "format": "string"
  }
}
```

#### 5.3 Search by Substructure

```
GET /rdkit/substructure?smarts={smarts}&limit={limit}
```

**Query Parameters:**
- `smarts` (required): SMARTS pattern to search for
- `limit` (optional): Maximum number of results to return (default: 50)

**Response:**
```json
{
  "status": "success",
  "data": {
    "matches": [
      {
        "id": "uuid",
        "name": "string",
        "smiles": "string",
        "similarity": "number"
      }
    ],
    "count": "number"
  }
}
```

#### 5.4 Search by Similarity

```
GET /rdkit/similarity?smiles={smiles}&threshold={threshold}&limit={limit}
```

**Query Parameters:**
- `smiles` (required): SMILES notation of the reference molecule
- `threshold` (optional): Similarity threshold (0-1, default: 0.7)
- `limit` (optional): Maximum number of results to return (default: 50)

**Response:**
```json
{
  "status": "success",
  "data": {
    "matches": [
      {
        "id": "uuid",
        "name": "string",
        "smiles": "string",
        "similarity": "number"
      }
    ],
    "count": "number"
  }
}
```

### 6. Authentication

#### 6.1 Login

```
POST /auth/login
```

**Request Body:**
```json
{
  "email": "string",
  "password": "string"
}
```

**Response:**
```json
{
  "status": "success",
  "data": {
    "message": "Authentication successful",
    "user": {
      "id": "uuid",
      "email": "string",
      "role": "string",
      "roles": ["string"],
      "permissions": ["string"]
    }
  }
}
```

#### 6.2 Logout

```
POST /auth/logout
```

**Response:**
```json
{
  "status": "success",
  "message": "Logout successful"
}
```

#### 6.3 Refresh Token

```
POST /auth/refresh
```

**Request Body:**
```json
{
  "refresh_token": "string" // Optional if using cookies
}
```

**Response:**
```json
{
  "status": "success",
  "data": {
    "message": "Session refreshed",
    "token": {
      "access_token": "string",
      "refresh_token": "string" // Only if not using cookies
    }
  }
}
```

#### 6.4 Register

```
POST /auth/register
```

**Request Body:**
```json
{
  "email": "string",
  "password": "string"
}
```

**Response:**
```json
{
  "status": "success",
  "data": {
    "message": "Registration successful",
    "user": {
      "id": "uuid",
      "email": "string"
    }
  }
}
```

#### 6.5 Reset Password

```
POST /auth/reset-password
```

**Request Body:**
```json
{
  "email": "string"
}
```

**Response:**
```json
{
  "status": "success",
  "message": "Password reset email sent"
}
```

### 7. Dashboard Data

#### 7.1 Get Dashboard Stats

```
GET /dashboard/stats
```

**Response:**
```json
{
  "status": "success",
  "data": {
    "molecules": {
      "total": "number",
      "consolidated": "number",
      "unique": "number"
    },
    "mixtures": {
      "total": "number",
      "public": "number",
      "private": "number"
    },
    "experiments": {
      "total": "number",
      "successful": "number",
      "failed": "number"
    },
    "predictions": {
      "total": "number",
      "high_confidence": "number",
      "medium_confidence": "number",
      "low_confidence": "number"
    }
  }
}
```

#### 7.2 Get Recent Activity

```
GET /dashboard/activity
```

**Response:**
```json
{
  "status": "success",
  "data": {
    "activities": [
      {
        "id": "uuid",
        "action": "string",
        "entity_type": "string",
        "entity_id": "uuid",
        "entity_name": "string",
        "timestamp": "timestamp",
        "user_id": "uuid",
        "user_email": "string"
      }
    ]
  }
}
```

## Response Format Standards

All endpoints follow a consistent response format:

### Success Response

```json
{
  "status": "success",
  "data": {
    // Endpoint-specific data
  },
  // Optional pagination info for list endpoints
  "pagination": {
    "limit": "number",
    "offset": "number",
    "total": "number"
  }
}
```

### Error Response

```json
{
  "status": "error",
  "error": {
    "code": "string",
    "message": "string",
    "details": "object" // Optional additional error details
  }
}
```

## Common Status Codes

- `200 OK`: Request succeeded
- `201 Created`: Resource created successfully
- `400 Bad Request`: Invalid input data
- `401 Unauthorized`: Authentication required
- `403 Forbidden`: Insufficient permissions
- `404 Not Found`: Resource not found
- `409 Conflict`: Resource conflict (e.g., duplicate data)
- `422 Unprocessable Entity`: Validation error
- `429 Too Many Requests`: Rate limit exceeded
- `500 Internal Server Error`: Server-side error

## API Rate Limits

- Standard user: 60 requests per minute
- Premium user: 300 requests per minute
- API key: 1000 requests per minute