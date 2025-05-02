# Roo Agent Database Population Directive

## Task Overview

Your mission is to implement a simplified, verifiable database population approach for the CryoProtect v2 project. This directive provides structured guidance for implementing the database population plan using Supabase MCP tools.

## Task Context

The CryoProtect v2 project requires a populated database with molecular data for cryoprotectant compounds. Previous attempts faced challenges with connection stability, transaction management, and property handling. We are shifting to a simplified approach focused on verifiable progress.

## Success Criteria

1. **Molecule Count**: At least 5,000 molecules in database
2. **Reference Compounds**: All 9 reference compounds with complete properties
3. **Property Completeness**: At least 90% of molecules have all required properties
4. **Query Performance**: Average query time under 50ms

## Task Structure

This task is designed to be executed in discrete phases, each with clear verification steps:

### Phase 1: Database Connection and Schema Validation
- **Purpose**: Ensure the database is accessible and properly structured
- **Actions**: 
  - Use MCP tools to connect to Supabase
  - Verify schema existence (tables, columns)
  - Ensure property types are defined
- **Verification**: All required tables exist and are properly structured
- **Time Budget**: 1 hour

### Phase 2: Reference Compound Population
- **Purpose**: Create a baseline set of critical compounds with complete properties
- **Actions**:
  - Insert 9 reference compounds (details in MCP_DATABASE_POPULATION_GUIDE.md)
  - Create all required properties for each compound
  - Verify property completeness
- **Verification**: All 9 reference compounds exist with complete properties
- **Time Budget**: 2 hours

### Phase 3: ChEMBL Data Import
- **Purpose**: Populate the database with molecular data from ChEMBL
- **Actions**:
  - Import compounds similar to reference compounds
  - Ensure property completeness during import
  - Focus on quality over quantity
- **Verification**: At least 5,000 molecules with complete properties
- **Time Budget**: 6 hours

### Phase 4: Performance Optimization
- **Purpose**: Ensure database queries meet performance requirements
- **Actions**:
  - Add indexes for common query patterns
  - Run database statistics update
  - Measure query performance
- **Verification**: Query time under 50ms for standard operations
- **Time Budget**: 2 hours

## Implementation Guidance

### Code Structure
- Implement each phase as a separate function
- Create clear validation functions for each verification step
- Use checkpointing to track progress
- Log all operations for debugging purposes

### Error Handling
- Implement retry logic for API requests
- Use transactions for database operations
- Log detailed error information
- Create recovery mechanisms for interrupted processes

### Using MCP Tools
- Always use `mcp__supabase__execute_sql` for database operations
- Begin by identifying the project ID with `mcp__supabase__list_projects`
- Use parameterized queries to prevent SQL injection
- Keep transactions small and focused

## Implementation Resources

Refer to these supporting documents for detailed implementation:

1. **MCP_DATABASE_POPULATION_GUIDE.md**: Comprehensive guide with SQL examples
2. **DATABASE_POPULATION_SIMPLIFIED_APPROACH.md**: Strategy overview
3. **MINIMAL_POPULATION_WORKFLOW_GUIDE.md**: Practical tips and workflows
4. **simplified_database_population.py**: Implementation example with detailed functions

## Sample Implementation (Phase 2 Example)

```python
# Example implementation of Reference Compounds population

def populate_reference_compounds(project_id):
    """
    Populate the 9 reference compounds with guaranteed properties.
    
    Args:
        project_id: The Supabase project ID
        
    Returns:
        dict: Statistics about the operation
    """
    # Reference compounds data
    reference_compounds = [
        {
            "name": "Glycerol", 
            "chembl_id": "CHEMBL388978", 
            "smiles": "C(C(CO)O)O", 
            "inchi": "InChI=1S/C3H8O3/c4-1-3(6)2-5/h3-6H,1-2H2", 
            "inchikey": "PEDCQBHIVMGVHV-UHFFFAOYSA-N", 
            "formula": "C3H8O3",
            "properties": {
                "logP": -1.76,
                "h_bond_donors": 3,
                "h_bond_acceptors": 3
            }
        },
        # Add other 8 reference compounds here...
    ]
    
    results = {
        "inserted": 0,
        "properties_added": 0,
        "errors": []
    }
    
    for compound in reference_compounds:
        try:
            # Insert molecule
            molecule_sql = """
            INSERT INTO molecules 
                (name, chembl_id, smiles, inchi, inchikey, formula, data_source)
            VALUES 
                (%s, %s, %s, %s, %s, %s, 'reference')
            ON CONFLICT (chembl_id) DO UPDATE SET
                name = EXCLUDED.name,
                smiles = EXCLUDED.smiles,
                inchi = EXCLUDED.inchi,
                inchikey = EXCLUDED.inchikey,
                formula = EXCLUDED.formula,
                updated_at = NOW()
            RETURNING id;
            """
            
            molecule_result = mcp__supabase__execute_sql(
                project_id=project_id,
                query=molecule_sql,
                params=[
                    compound["name"], 
                    compound["chembl_id"], 
                    compound["smiles"],
                    compound["inchi"], 
                    compound["inchikey"], 
                    compound["formula"]
                ]
            )
            
            if molecule_result and len(molecule_result) > 0:
                molecule_id = molecule_result[0]["id"]
                results["inserted"] += 1
                
                # Insert properties
                for prop_name, value in compound["properties"].items():
                    # Get or create property type
                    prop_type_sql = """
                    INSERT INTO property_types (name, data_type)
                    VALUES (%s, 'numeric')
                    ON CONFLICT (name) DO UPDATE SET updated_at = NOW()
                    RETURNING id;
                    """
                    
                    prop_type_result = mcp__supabase__execute_sql(
                        project_id=project_id,
                        query=prop_type_sql,
                        params=[prop_name]
                    )
                    
                    if prop_type_result and len(prop_type_result) > 0:
                        property_type_id = prop_type_result[0]["id"]
                        
                        # Insert property value
                        property_sql = """
                        INSERT INTO molecular_properties (molecule_id, property_type_id, numeric_value)
                        VALUES (%s, %s, %s)
                        ON CONFLICT (molecule_id, property_type_id) DO UPDATE SET
                            numeric_value = EXCLUDED.numeric_value,
                            updated_at = NOW();
                        """
                        
                        mcp__supabase__execute_sql(
                            project_id=project_id,
                            query=property_sql,
                            params=[molecule_id, property_type_id, value]
                        )
                        
                        results["properties_added"] += 1
                
        except Exception as e:
            results["errors"].append({
                "compound": compound["chembl_id"],
                "error": str(e)
            })
    
    return results
```

## Sample Implementation (Phase 3 Example - ChEMBL Import)

```python
def import_chembl_data(project_id, limit=5000, batch_size=50):
    """
    Import chemical data from ChEMBL focusing on cryoprotectant compounds.
    
    Args:
        project_id: The Supabase project ID
        limit: Maximum number of compounds to import
        batch_size: Number of compounds to process in each batch
        
    Returns:
        dict: Statistics about the import operation
    """
    # Define the ChEMBL search terms for cryoprotectants
    search_terms = [
        "cryoprotectant", 
        "cryoprotection", 
        "cryopreservation", 
        "antifreeze", 
        "freeze protection",
        "polyol", 
        "glycol"
    ]
    
    results = {
        "molecules_imported": 0,
        "properties_imported": 0,
        "errors": [],
        "batches_processed": 0
    }
    
    # Iterate through search terms
    for term in search_terms:
        if results["molecules_imported"] >= limit:
            break
            
        try:
            # Search for compounds matching the term
            search_url = f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_properties__full_mwt__lte=500&limit=1000&q={term}"
            response = requests.get(search_url)
            
            if response.status_code != 200:
                results["errors"].append({
                    "term": term,
                    "error": f"API request failed: {response.status_code}"
                })
                continue
                
            data = response.json()
            molecules = data.get("molecules", [])
            
            # Process in batches
            for i in range(0, len(molecules), batch_size):
                if results["molecules_imported"] >= limit:
                    break
                    
                batch = molecules[i:i+batch_size]
                batch_results = process_chembl_batch(project_id, batch)
                
                # Update results
                results["molecules_imported"] += batch_results["molecules_imported"]
                results["properties_imported"] += batch_results["properties_imported"]
                results["errors"].extend(batch_results["errors"])
                results["batches_processed"] += 1
                
                # Create checkpoint
                create_checkpoint(project_id, {
                    "term": term,
                    "batch": i // batch_size,
                    "molecules_imported": results["molecules_imported"],
                    "properties_imported": results["properties_imported"]
                })
                
                # Log progress
                print(f"Progress: {results['molecules_imported']}/{limit} molecules imported")
                
        except Exception as e:
            results["errors"].append({
                "term": term,
                "error": str(e)
            })
    
    return results

def process_chembl_batch(project_id, batch):
    """
    Process a batch of ChEMBL molecules and import into database.
    
    Args:
        project_id: The Supabase project ID
        batch: List of molecule data from ChEMBL API
        
    Returns:
        dict: Statistics about the batch processing
    """
    results = {
        "molecules_imported": 0,
        "properties_imported": 0,
        "errors": []
    }
    
    for molecule_data in batch:
        try:
            # Skip if missing required fields
            if not all(field in molecule_data for field in ["molecule_chembl_id", "molecule_structures"]):
                continue
                
            if not molecule_data.get("molecule_structures", {}).get("canonical_smiles"):
                continue
                
            # Extract basic properties
            molecule = {
                "chembl_id": molecule_data["molecule_chembl_id"],
                "name": molecule_data.get("pref_name", molecule_data["molecule_chembl_id"]),
                "smiles": molecule_data["molecule_structures"]["canonical_smiles"],
                "inchi": molecule_data["molecule_structures"].get("standard_inchi", ""),
                "inchikey": molecule_data["molecule_structures"].get("standard_inchi_key", ""),
                "formula": molecule_data.get("molecule_properties", {}).get("full_molformula", ""),
                "molecular_weight": molecule_data.get("molecule_properties", {}).get("full_mwt"),
                "data_source": "ChEMBL"
            }
            
            # Insert molecule
            molecule_sql = """
            INSERT INTO molecules 
                (chembl_id, name, smiles, inchi, inchikey, formula, molecular_weight, data_source)
            VALUES 
                (%s, %s, %s, %s, %s, %s, %s, %s)
            ON CONFLICT (chembl_id) DO UPDATE SET
                name = EXCLUDED.name,
                smiles = EXCLUDED.smiles,
                inchi = EXCLUDED.inchi,
                inchikey = EXCLUDED.inchikey,
                formula = EXCLUDED.formula,
                molecular_weight = EXCLUDED.molecular_weight,
                updated_at = NOW()
            RETURNING id;
            """
            
            molecule_result = mcp__supabase__execute_sql(
                project_id=project_id,
                query=molecule_sql,
                params=[
                    molecule["chembl_id"],
                    molecule["name"],
                    molecule["smiles"],
                    molecule["inchi"],
                    molecule["inchikey"],
                    molecule["formula"],
                    molecule["molecular_weight"],
                    molecule["data_source"]
                ]
            )
            
            if not molecule_result or len(molecule_result) == 0:
                continue
                
            molecule_id = molecule_result[0]["id"]
            results["molecules_imported"] += 1
            
            # Extract and insert properties
            properties = {}
            mol_props = molecule_data.get("molecule_properties", {})
            
            # Map ChEMBL properties to our property names
            if "alogp" in mol_props:
                properties["logP"] = mol_props["alogp"]
            if "hba" in mol_props:
                properties["h_bond_acceptors"] = mol_props["hba"]
            if "hbd" in mol_props:
                properties["h_bond_donors"] = mol_props["hbd"]
                
            # Add properties
            for prop_name, value in properties.items():
                try:
                    # Get or create property type
                    prop_type_sql = """
                    INSERT INTO property_types (name, data_type)
                    VALUES (%s, 'numeric')
                    ON CONFLICT (name) DO UPDATE SET updated_at = NOW()
                    RETURNING id;
                    """
                    
                    prop_type_result = mcp__supabase__execute_sql(
                        project_id=project_id,
                        query=prop_type_sql,
                        params=[prop_name]
                    )
                    
                    if prop_type_result and len(prop_type_result) > 0:
                        property_type_id = prop_type_result[0]["id"]
                        
                        # Insert property value
                        property_sql = """
                        INSERT INTO molecular_properties (molecule_id, property_type_id, numeric_value)
                        VALUES (%s, %s, %s)
                        ON CONFLICT (molecule_id, property_type_id) DO UPDATE SET
                            numeric_value = EXCLUDED.numeric_value,
                            updated_at = NOW();
                        """
                        
                        mcp__supabase__execute_sql(
                            project_id=project_id,
                            query=property_sql,
                            params=[molecule_id, property_type_id, value]
                        )
                        
                        results["properties_imported"] += 1
                        
                except Exception as e:
                    results["errors"].append({
                        "molecule_id": molecule_id,
                        "property": prop_name,
                        "error": str(e)
                    })
                    
        except Exception as e:
            results["errors"].append({
                "molecule": molecule_data.get("molecule_chembl_id", "unknown"),
                "error": str(e)
            })
    
    return results

def create_checkpoint(project_id, data):
    """
    Create a checkpoint to track import progress.
    
    Args:
        project_id: The Supabase project ID
        data: Checkpoint data to store
    """
    checkpoint_sql = """
    INSERT INTO import_checkpoints 
        (import_type, checkpoint_data, created_at)
    VALUES 
        ('chembl', %s, NOW())
    RETURNING id;
    """
    
    mcp__supabase__execute_sql(
        project_id=project_id,
        query=checkpoint_sql,
        params=[json.dumps(data)]
    )
```

## Reporting Requirements

Provide a detailed report after each phase including:

1. **Progress Statistics**:
   - Number of molecules successfully added
   - Number of properties populated
   - Percentage of completion
   
2. **Verification Results**:
   - Results of all verification queries
   - Status of each success criterion
   
3. **Technical Details**:
   - Any modifications to the approach
   - Performance metrics
   - Issues encountered and solutions

4. **Next Steps**:
   - Recommended improvements
   - Areas requiring further attention

## Integration Points

This task interacts with these other components:

1. **Database Schema**: Defined in migrations/001_initial_schema.sql
2. **Property Utilities**: See property_utils.py for reference
3. **Verification Scripts**: verify_imported_data.py contains verification logic

## Final Deliverables

1. Complete implementation of all 4 phases
2. Verification that all success criteria are met
3. Documentation of the implementation approach
4. SQL scripts for any schema modifications
5. Performance optimization recommendations

## Timeframe

Complete this task within 11 hours total (spread across the 4 phases).

## Communication Protocol

- Report progress after each phase completion
- Flag any blocking issues immediately
- Document any deviations from the plan with justification
- Provide regular status updates (hourly) during implementation

## Approval Criteria

Your implementation will be considered successful when:
1. All 4 success criteria are fully met
2. The implementation follows the simplified approach
3. No critical bugs or issues are present
4. The performance meets or exceeds requirements