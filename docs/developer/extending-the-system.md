# Extending the System

This guide provides information on how to extend the CryoProtect system with new features and functionality.

## Overview

The CryoProtect system has been designed with extensibility in mind. You can extend the system in various ways:

1. Adding new tables to the database
2. Extending existing functionality
3. Creating new API endpoints
4. Adding new molecular properties
5. Integrating with external systems

This guide provides detailed instructions for each of these extension points.

## Adding New Tables

### Database Schema Considerations

When adding new tables to the CryoProtect database, follow these guidelines:

1. **Use Plural Table Names**: All table names should be plural (e.g., `experiments`, not `experiment`)
2. **Use UUID Primary Keys**: All tables should have a UUID primary key with DEFAULT gen_random_uuid()
3. **Include Standard Fields**: Include standard fields like created_at, updated_at, and created_by
4. **Add Foreign Key Constraints**: Ensure proper relationships with foreign key constraints
5. **Enable Row Level Security**: Enable RLS and add appropriate policies

### Example: Adding a New Table

Here's an example of adding a new `protocols` table:

```sql
-- Create the table
CREATE TABLE IF NOT EXISTS public.protocols (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    name VARCHAR(255) NOT NULL,
    description TEXT,
    steps JSONB,
    is_public BOOLEAN DEFAULT false,
    project_id UUID REFERENCES public.projects(id) ON DELETE CASCADE,
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    created_by UUID REFERENCES auth.users(id)
);

-- Create indexes
CREATE INDEX IF NOT EXISTS idx_protocols_project_id ON public.protocols(project_id);
CREATE INDEX IF NOT EXISTS idx_protocols_created_by ON public.protocols(created_by);

-- Enable RLS
ALTER TABLE public.protocols ENABLE ROW LEVEL SECURITY;

-- Add RLS policies
CREATE POLICY "Public read access on protocols" 
ON public.protocols 
FOR SELECT 
USING (is_public = true);

CREATE POLICY "Owner full access on protocols" 
ON public.protocols 
USING (created_by = auth.uid());

CREATE POLICY "Team member access on protocols" 
ON public.protocols 
FOR SELECT 
USING (
  EXISTS (
    SELECT 1 FROM projects
    JOIN teams ON teams.id = projects.team_id
    JOIN user_profile ON user_profile.team_id = teams.id
    WHERE projects.id = protocols.project_id
    AND user_profile.auth_user_id = auth.uid()
  )
);

CREATE POLICY "Service role bypass on protocols" 
ON public.protocols 
USING (auth.role() = 'service_role');
```

### Python Model for New Table

Create a Python model for the new table in the `models` directory:

```python
# models/protocol.py
from datetime import datetime
from uuid import UUID, uuid4
from typing import Optional, List, Dict, Any
from pydantic import BaseModel, Field

class Protocol(BaseModel):
    id: UUID = Field(default_factory=uuid4)
    name: str
    description: Optional[str] = None
    steps: Optional[Dict[str, Any]] = None
    is_public: bool = False
    project_id: UUID
    created_at: datetime = Field(default_factory=datetime.now)
    updated_at: datetime = Field(default_factory=datetime.now)
    created_by: Optional[UUID] = None
    
    class Config:
        orm_mode = True
```

## Extending Existing Functionality

### Adding Fields to Existing Tables

To add fields to existing tables:

1. Create a migration script in the `migrations` directory
2. Update the corresponding Python model
3. Update any affected API endpoints
4. Update the UI to display/edit the new fields

Example migration script:

```sql
-- Add a new field to the molecules table
ALTER TABLE public.molecules
ADD COLUMN IF NOT EXISTS toxicity_level VARCHAR(50);

-- Add an index for the new field
CREATE INDEX IF NOT EXISTS idx_molecules_toxicity_level
ON public.molecules(toxicity_level);
```

### Extending Existing Classes

To extend existing classes with new functionality:

1. Create a subclass that inherits from the existing class
2. Override methods as needed
3. Add new methods for the additional functionality

Example:

```python
# models/enhanced_molecule.py
from models.molecule import Molecule
from typing import Optional, List

class EnhancedMolecule(Molecule):
    toxicity_level: Optional[str] = None
    
    def calculate_toxicity(self) -> str:
        """Calculate toxicity level based on molecular properties."""
        # Implementation here
        return "Low"  # Example return value
```

## Creating New API Endpoints

### Adding New Resources

To add new API endpoints:

1. Create a new resource class in `api/resources.py`
2. Register the resource in `api/__init__.py`
3. Implement the required HTTP methods (GET, POST, PUT, DELETE)

Example resource class:

```python
# api/resources.py
from flask_restful import Resource, reqparse
from flask import request
from models.protocol import Protocol
from utils.database import supabase

class ProtocolResource(Resource):
    def get(self, protocol_id=None):
        """Get a protocol or list of protocols."""
        if protocol_id:
            # Get a specific protocol
            result = supabase.table("protocols").select("*").eq("id", protocol_id).execute()
            if not result.data:
                return {"error": "Protocol not found"}, 404
            return result.data[0], 200
        else:
            # Get all protocols
            result = supabase.table("protocols").select("*").execute()
            return result.data, 200
    
    def post(self):
        """Create a new protocol."""
        data = request.get_json()
        # Validate data
        # ...
        
        # Insert into database
        result = supabase.table("protocols").insert(data).execute()
        return result.data[0], 201
    
    def put(self, protocol_id):
        """Update a protocol."""
        if not protocol_id:
            return {"error": "Protocol ID is required"}, 400
            
        data = request.get_json()
        # Validate data
        # ...
        
        # Update in database
        result = supabase.table("protocols").update(data).eq("id", protocol_id).execute()
        if not result.data:
            return {"error": "Protocol not found"}, 404
        return result.data[0], 200
    
    def delete(self, protocol_id):
        """Delete a protocol."""
        if not protocol_id:
            return {"error": "Protocol ID is required"}, 400
            
        # Delete from database
        result = supabase.table("protocols").delete().eq("id", protocol_id).execute()
        if not result.data:
            return {"error": "Protocol not found"}, 404
        return {"message": "Protocol deleted"}, 200
```

Register the resource:

```python
# api/__init__.py
from flask_restful import Api
from .resources import ProtocolResource

# ... existing code ...

# Register resources
api.add_resource(ProtocolResource, '/protocols', '/protocols/<string:protocol_id>')
```

## Adding New Molecular Properties

### Implementing New Property Calculations

To add new molecular property calculations:

1. Create a new function in the appropriate module
2. Register the property in the property_types table
3. Update the UI to display the new property

Example property calculation:

```python
# models/molecular_properties.py
from rdkit import Chem
from rdkit.Chem import Descriptors

def calculate_toxicity(mol_smiles):
    """Calculate toxicity level based on molecular properties."""
    mol = Chem.MolFromSmiles(mol_smiles)
    if not mol:
        return None
        
    # Example: Use logP as a simple toxicity indicator
    logp = Descriptors.MolLogP(mol)
    
    if logp < 0:
        return "Low"
    elif logp < 3:
        return "Medium"
    else:
        return "High"
```

Register the property:

```sql
INSERT INTO public.property_types (name, data_type, description)
VALUES ('toxicity_level', 'string', 'Toxicity level based on molecular properties')
ON CONFLICT (name) DO NOTHING;
```

## Integrating with External Systems

### API Integration

To integrate with external APIs:

1. Create a client class for the external API
2. Implement methods to interact with the API
3. Handle authentication and error cases

Example external API client:

```python
# utils/external_api.py
import requests
from typing import Dict, Any, Optional

class ExternalAPIClient:
    def __init__(self, base_url, api_key):
        self.base_url = base_url
        self.api_key = api_key
        self.session = requests.Session()
        self.session.headers.update({
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json"
        })
    
    def get_data(self, endpoint, params=None):
        """Get data from the external API."""
        url = f"{self.base_url}/{endpoint}"
        response = self.session.get(url, params=params)
        response.raise_for_status()
        return response.json()
    
    def post_data(self, endpoint, data):
        """Post data to the external API."""
        url = f"{self.base_url}/{endpoint}"
        response = self.session.post(url, json=data)
        response.raise_for_status()
        return response.json()
```

### Database Integration

To integrate with external databases:

1. Create a connection class for the external database
2. Implement methods to query and update data
3. Handle connection pooling and error cases

Example external database connection:

```python
# utils/external_db.py
import psycopg2
from psycopg2.pool import SimpleConnectionPool
from contextlib import contextmanager

class ExternalDBConnection:
    def __init__(self, dbname, user, password, host, port=5432, pool_size=5):
        self.pool = SimpleConnectionPool(
            1, pool_size,
            dbname=dbname,
            user=user,
            password=password,
            host=host,
            port=port
        )
    
    @contextmanager
    def get_connection(self):
        """Get a connection from the pool."""
        conn = self.pool.getconn()
        try:
            yield conn
        finally:
            self.pool.putconn(conn)
    
    def query(self, sql, params=None):
        """Execute a query and return the results."""
        with self.get_connection() as conn:
            with conn.cursor() as cur:
                cur.execute(sql, params or ())
                return cur.fetchall()
    
    def execute(self, sql, params=None):
        """Execute a statement and return the number of affected rows."""
        with self.get_connection() as conn:
            with conn.cursor() as cur:
                cur.execute(sql, params or ())
                conn.commit()
                return cur.rowcount
```

## Best Practices for Extensions

When extending the CryoProtect system, follow these best practices:

1. **Maintain Consistency**: Follow the existing patterns and conventions
2. **Write Tests**: Create tests for your extensions
3. **Document Changes**: Update documentation to reflect your changes
4. **Consider Security**: Ensure your extensions maintain the security model
5. **Handle Errors**: Implement proper error handling
6. **Validate Input**: Validate all user input
7. **Use Transactions**: Use database transactions for multi-step operations
8. **Follow 3NF**: Ensure your database changes follow Third Normal Form
9. **Add Indexes**: Add appropriate indexes for performance
10. **Enable RLS**: Enable Row Level Security on new tables

## Example: Complete Extension

Here's a complete example of adding a new feature to track protocol execution results:

1. **Database Schema**:
   ```sql
   -- Create the protocol_results table
   CREATE TABLE IF NOT EXISTS public.protocol_results (
       id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
       protocol_id UUID NOT NULL REFERENCES public.protocols(id) ON DELETE CASCADE,
       experiment_id UUID NOT NULL REFERENCES public.experiments(id) ON DELETE CASCADE,
       result_data JSONB,
       status VARCHAR(50) NOT NULL,
       notes TEXT,
       is_public BOOLEAN DEFAULT false,
       created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
       updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
       created_by UUID REFERENCES auth.users(id),
       UNIQUE(protocol_id, experiment_id)
   );
   
   -- Create indexes
   CREATE INDEX IF NOT EXISTS idx_protocol_results_protocol_id ON public.protocol_results(protocol_id);
   CREATE INDEX IF NOT EXISTS idx_protocol_results_experiment_id ON public.protocol_results(experiment_id);
   CREATE INDEX IF NOT EXISTS idx_protocol_results_status ON public.protocol_results(status);
   
   -- Enable RLS
   ALTER TABLE public.protocol_results ENABLE ROW LEVEL SECURITY;
   
   -- Add RLS policies
   CREATE POLICY "Public read access on protocol_results" 
   ON public.protocol_results 
   FOR SELECT 
   USING (is_public = true);
   
   CREATE POLICY "Owner full access on protocol_results" 
   ON public.protocol_results 
   USING (created_by = auth.uid());
   
   CREATE POLICY "Team member access on protocol_results" 
   ON public.protocol_results 
   FOR SELECT 
   USING (
     EXISTS (
       SELECT 1 FROM protocols
       JOIN projects ON projects.id = protocols.project_id
       JOIN teams ON teams.id = projects.team_id
       JOIN user_profile ON user_profile.team_id = teams.id
       WHERE protocols.id = protocol_results.protocol_id
       AND user_profile.auth_user_id = auth.uid()
     )
   );
   
   CREATE POLICY "Service role bypass on protocol_results" 
   ON public.protocol_results 
   USING (auth.role() = 'service_role');
   ```

2. **Python Model**:
   ```python
   # models/protocol_result.py
   from datetime import datetime
   from uuid import UUID, uuid4
   from typing import Optional, Dict, Any
   from pydantic import BaseModel, Field
   
   class ProtocolResult(BaseModel):
       id: UUID = Field(default_factory=uuid4)
       protocol_id: UUID
       experiment_id: UUID
       result_data: Optional[Dict[str, Any]] = None
       status: str
       notes: Optional[str] = None
       is_public: bool = False
       created_at: datetime = Field(default_factory=datetime.now)
       updated_at: datetime = Field(default_factory=datetime.now)
       created_by: Optional[UUID] = None
       
       class Config:
           orm_mode = True
   ```

3. **API Resource**:
   ```python
   # api/resources.py
   from flask_restful import Resource
   from flask import request
   from models.protocol_result import ProtocolResult
   from utils.database import supabase
   
   class ProtocolResultResource(Resource):
       def get(self, result_id=None):
           """Get a protocol result or list of results."""
           if result_id:
               # Get a specific result
               result = supabase.table("protocol_results").select("*").eq("id", result_id).execute()
               if not result.data:
                   return {"error": "Protocol result not found"}, 404
               return result.data[0], 200
           else:
               # Get all results
               result = supabase.table("protocol_results").select("*").execute()
               return result.data, 200
       
       def post(self):
           """Create a new protocol result."""
           data = request.get_json()
           # Validate data
           # ...
           
           # Insert into database
           result = supabase.table("protocol_results").insert(data).execute()
           return result.data[0], 201
       
       def put(self, result_id):
           """Update a protocol result."""
           if not result_id:
               return {"error": "Result ID is required"}, 400
               
           data = request.get_json()
           # Validate data
           # ...
           
           # Update in database
           result = supabase.table("protocol_results").update(data).eq("id", result_id).execute()
           if not result.data:
               return {"error": "Protocol result not found"}, 404
           return result.data[0], 200
       
       def delete(self, result_id):
           """Delete a protocol result."""
           if not result_id:
               return {"error": "Result ID is required"}, 400
               
           # Delete from database
           result = supabase.table("protocol_results").delete().eq("id", result_id).execute()
           if not result.data:
               return {"error": "Protocol result not found"}, 404
           return {"message": "Protocol result deleted"}, 200
   ```

4. **Register the Resource**:
   ```python
   # api/__init__.py
   from flask_restful import Api
   from .resources import ProtocolResultResource
   
   # ... existing code ...
   
   # Register resources
   api.add_resource(ProtocolResultResource, '/protocol-results', '/protocol-results/<string:result_id>')
   ```

5. **UI Integration**:
   - Create HTML templates for viewing and editing protocol results
   - Add JavaScript for handling form submission and data display
   - Update navigation to include the new feature

By following these guidelines, you can extend the CryoProtect system with new features while maintaining consistency with the existing architecture and design principles.