# Phase 3 Planning: Further Data Quality Enhancements

Based on the successful implementation of Phase 2 molecule consolidation, here's a plan for Phase 3 of the data quality enhancement project.

## Phase 3 Objectives

1. **Complete the Molecule Consolidation**
   - Implement the "differentiate" strategy for molecules with different structures
   - Fix prediction migration issues in complex consolidation
   - Create the necessary database triggers to maintain consolidation integrity

2. **Application Integration**
   - Update API endpoints to handle consolidated molecules
   - Enhance the user interface to show consolidation information
   - Add endpoints to manage and modify consolidation relationships

3. **Further Data Enrichment**
   - Complete missing property values for all molecules
   - Enhance property verification and validation
   - Add additional data sources and cross-references

## Implementation Plan

### 1. Complete the "Differentiate" Strategy

Molecules that were identified as similar but having different structures need special handling:

```python
def process_differentiate_group(conn, group):
    """
    Process a differentiate group.
    These are molecules that appeared similar but have different structures.
    We'll add descriptive properties to clearly differentiate them.
    """
    for molecule in group['molecules']:
        # Add properties to clarify the differences
        properties = {
            'differentiation_group': group['group_id'],
            'chemical_identity': molecule.get('molecular_formula', 'Unknown'),
            'structure_note': f"SMILES: {molecule.get('smiles', 'Unknown')}"
        }
        
        # Update the molecule with differentiation properties
        with conn.cursor() as cursor:
            cursor.execute("""
                UPDATE molecules
                SET properties = properties || %s::jsonb
                WHERE id = %s
            """, (json.dumps(properties), molecule['id']))
```

### 2. Fix Prediction Migration

The current complex consolidation script has issues with prediction migration. The fix would be:

```python
def migrate_predictions(cursor, source_id, target_id, consolidation_log):
    """
    Migrate predictions from source to target molecule.
    """
    # Get all predictions for the source molecule
    cursor.execute("""
        SELECT * FROM predictions
        WHERE molecule_id = %s
    """, (source_id,))
    
    predictions = cursor.fetchall()
    
    for prediction in predictions:
        property_type_id = prediction.get('property_type_id')
        
        # Check if a prediction for this property already exists for the target
        cursor.execute("""
            SELECT 1 FROM predictions
            WHERE molecule_id = %s AND property_type_id = %s
        """, (target_id, property_type_id))
        
        if cursor.fetchone():
            continue
        
        # Migrate the prediction to the target molecule
        new_prediction_id = str(uuid.uuid4())
        
        # Handle results field properly - convert to JSON if needed
        results = prediction.get('results')
        results_json = None
        
        if results:
            if isinstance(results, dict) or isinstance(results, list):
                results_json = Json(results)
            else:
                try:
                    results_dict = json.loads(results)
                    results_json = Json(results_dict)
                except:
                    # If parsing fails, store as is
                    results_json = Json({})
        
        # Similar handling for modification_history
        mod_history = prediction.get('modification_history')
        mod_history_json = None
        
        if mod_history:
            if isinstance(mod_history, dict):
                mod_history_dict = mod_history.copy()
            else:
                try:
                    mod_history_dict = json.loads(mod_history)
                except:
                    mod_history_dict = {}
                    
            mod_history_dict["consolidated_from"] = str(source_id)
            mod_history_json = Json(mod_history_dict)
        
        # Insert the new prediction
        cursor.execute("""
            INSERT INTO predictions (
                id, molecule_id, mixture_id, property_type_id, calculation_method_id,
                numeric_value, text_value, boolean_value, confidence,
                created_at, updated_at, created_by, data_source, version,
                modification_history, name, description, results
            )
            VALUES (
                %s, %s, %s, %s, %s,
                %s, %s, %s, %s,
                NOW(), NOW(), %s, %s, %s,
                %s, %s, %s, %s
            )
            RETURNING id
        """, (
            new_prediction_id, target_id, prediction.get('mixture_id'),
            prediction.get('property_type_id'), prediction.get('calculation_method_id'),
            prediction.get('numeric_value'), prediction.get('text_value'),
            prediction.get('boolean_value'), prediction.get('confidence'),
            prediction.get('created_by'), 
            f"Consolidated from {source_id}" if not prediction.get('data_source') else prediction.get('data_source'), 
            prediction.get('version'),
            mod_history_json, prediction.get('name'), prediction.get('description'),
            results_json
        ))
```

### 3. Database Triggers for Consolidation Integrity

To ensure that new relationships are created with primary molecules, we can use database triggers:

```sql
-- Function to redirect to primary molecule
CREATE OR REPLACE FUNCTION redirect_to_primary_molecule()
RETURNS TRIGGER AS $$
DECLARE
    primary_id UUID;
BEGIN
    -- Get the primary molecule ID
    SELECT get_primary_molecule_id(NEW.molecule_id) INTO primary_id;
    
    -- If it's different from the provided ID, use the primary
    IF primary_id != NEW.molecule_id THEN
        NEW.molecule_id := primary_id;
    END IF;
    
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

-- Trigger for mixture components
CREATE TRIGGER use_primary_molecule_in_mixture_components
BEFORE INSERT OR UPDATE ON mixture_components
FOR EACH ROW EXECUTE FUNCTION redirect_to_primary_molecule();

-- Trigger for predictions
CREATE TRIGGER use_primary_molecule_in_predictions
BEFORE INSERT OR UPDATE ON predictions
FOR EACH ROW EXECUTE FUNCTION redirect_to_primary_molecule();

-- Trigger for molecular properties
CREATE TRIGGER use_primary_molecule_in_molecular_properties
BEFORE INSERT OR UPDATE ON molecular_properties
FOR EACH ROW EXECUTE FUNCTION redirect_to_primary_molecule();
```

### 4. API Enhancements

The API should be enhanced to work with consolidated molecules:

```python
@app.route('/api/v1/molecules/<uuid:molecule_id>', methods=['GET'])
def get_molecule(molecule_id):
    """Get a molecule by ID, handling consolidation."""
    conn = get_db_connection()
    
    try:
        # Check if molecule is consolidated
        cursor = conn.cursor(cursor_factory=RealDictCursor)
        cursor.execute("""
            SELECT properties->>'consolidated_to' AS primary_id
            FROM molecules
            WHERE id = %s
        """, (str(molecule_id),))
        
        result = cursor.fetchone()
        
        if not result:
            return jsonify({'error': 'Molecule not found'}), 404
            
        primary_id = result.get('primary_id')
        
        # If consolidated, include consolidation info
        if primary_id:
            cursor.execute("""
                SELECT * FROM molecules
                WHERE id = %s
            """, (primary_id,))
            
            primary_molecule = cursor.fetchone()
            
            cursor.execute("""
                SELECT * FROM molecules
                WHERE id = %s
            """, (str(molecule_id),))
            
            secondary_molecule = cursor.fetchone()
            
            return jsonify({
                'molecule': secondary_molecule,
                'is_consolidated': True,
                'primary_molecule': primary_molecule
            })
        
        # Not consolidated, return as is
        cursor.execute("""
            SELECT * FROM molecules
            WHERE id = %s
        """, (str(molecule_id),))
        
        molecule = cursor.fetchone()
        
        # Check if it has secondary molecules
        cursor.execute("""
            SELECT * FROM molecules
            WHERE properties->>'consolidated_to' = %s
        """, (str(molecule_id),))
        
        secondary_molecules = cursor.fetchall()
        
        return jsonify({
            'molecule': molecule,
            'is_consolidated': False,
            'secondary_molecules': secondary_molecules
        })
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500
    
    finally:
        conn.close()
```

### 5. UI Enhancements

The UI should be enhanced to show consolidation information:

```javascript
// Example React component for displaying molecule info with consolidation
function MoleculeDetail({ moleculeId }) {
  const [molecule, setMolecule] = useState(null);
  const [isLoading, setIsLoading] = useState(true);
  
  useEffect(() => {
    async function fetchData() {
      try {
        const response = await fetch(`/api/v1/molecules/${moleculeId}`);
        const data = await response.json();
        setMolecule(data);
        setIsLoading(false);
      } catch (error) {
        console.error("Error fetching molecule:", error);
        setIsLoading(false);
      }
    }
    
    fetchData();
  }, [moleculeId]);
  
  if (isLoading) return <div>Loading...</div>;
  if (!molecule) return <div>Molecule not found</div>;
  
  return (
    <div className="molecule-detail">
      <h2>{molecule.molecule.name}</h2>
      
      {molecule.is_consolidated && (
        <div className="consolidation-notice">
          <p>This molecule has been consolidated with:</p>
          <h3>{molecule.primary_molecule.name}</h3>
          <p>All properties and relationships will use the primary molecule.</p>
          <Link to={`/molecules/${molecule.primary_molecule.id}`}>
            View Primary Molecule
          </Link>
        </div>
      )}
      
      {!molecule.is_consolidated && molecule.secondary_molecules.length > 0 && (
        <div className="primary-notice">
          <p>This is a primary molecule with {molecule.secondary_molecules.length} consolidated variants:</p>
          <ul>
            {molecule.secondary_molecules.map(sec => (
              <li key={sec.id}>
                <Link to={`/molecules/${sec.id}`}>{sec.name}</Link>
              </li>
            ))}
          </ul>
        </div>
      )}
      
      {/* Rest of molecule detail display */}
    </div>
  );
}
```

## Timeline and Resources

### Week 1-2: Complete Consolidation
- Implement differentiate strategy
- Fix prediction migration issues
- Set up database triggers

### Week 3-4: Application Integration
- Update API endpoints
- Create new consolidation management endpoints
- Implement UI changes

### Week 5-6: Data Enrichment
- Complete missing properties
- Add property validation
- Implement cross-references

### Resources Needed
- Back-end developer for database and API work
- Front-end developer for UI integration
- Database administrator for trigger implementation
- Computational chemist for property validation

## Success Criteria

1. All duplicate molecules properly consolidated
2. API endpoints correctly handle consolidated molecules
3. UI clearly shows consolidation relationships
4. Database triggers maintain data integrity
5. Applications continue to function with the consolidated data model

## Monitoring and Validation

After implementation, we should:

1. Set up database monitoring to detect any issues with the triggers
2. Create periodic validation checks to ensure consolidation integrity
3. Analyze application logs to detect any issues with the API endpoints
4. Gather user feedback on the UI changes

This Phase 3 plan builds on the successful work in Phase 2 and completes the molecule data quality enhancement project.