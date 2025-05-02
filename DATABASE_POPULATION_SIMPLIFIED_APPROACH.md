# Database Population Simplified Approach

## Core Philosophy: Simplify, Verify, Then Optimize

Based on metacognitive analysis, we're shifting from complex optimizations to ensuring basic functionality first, with clear verification at each step.

## Phase 1: Minimum Viable Database Population

### Step 1: Focus on Reference Compounds Only
```python
def populate_reference_compounds():
    """
    Populate ONLY the 9 reference compounds with guaranteed properties.
    This is our simplest test of the end-to-end pipeline.
    """
    # Hard-coded reference compounds with complete properties
    reference_compounds = [
        {
            "chembl_id": "CHEMBL388978",
            "name": "Glycerol",
            "smiles": "C(C(CO)O)O",
            "inchi": "InChI=1S/C3H8O3/c4-1-3(6)2-5/h3-6H,1-2H2",
            "inchikey": "PEDCQBHIVMGVHV-UHFFFAOYSA-N",
            "formula": "C3H8O3",
            "properties": {
                "logP": -1.76,
                "h_bond_donors": 3,
                "h_bond_acceptors": 3,
                "molecular_weight": 92.09
            }
        },
        {
            "chembl_id": "CHEMBL1098659",
            "name": "DMSO",
            "smiles": "CS(=O)C",
            "inchi": "InChI=1S/C2H6OS/c1-4(2)3/h1-2H3",
            "inchikey": "IAZDPXIOMUYVGZ-UHFFFAOYSA-N",
            "formula": "C2H6OS",
            "properties": {
                "logP": -1.35,
                "h_bond_donors": 0,
                "h_bond_acceptors": 1,
                "molecular_weight": 78.13
            }
        },
        {
            "chembl_id": "CHEMBL66195",
            "name": "beta-Alanine",
            "smiles": "C(CC(=O)O)N",
            "inchi": "InChI=1S/C3H7NO2/c4-2-1-3(5)6/h1-2,4H2,(H,5,6)",
            "inchikey": "UCMIRNVEIXFBKS-UHFFFAOYSA-N",
            "formula": "C3H7NO2",
            "properties": {
                "logP": -3.17,
                "h_bond_donors": 2,
                "h_bond_acceptors": 3,
                "molecular_weight": 89.09
            }
        },
        {
            "chembl_id": "CHEMBL500033",
            "name": "tert-Butanol",
            "smiles": "CC(C)(C)O",
            "inchi": "InChI=1S/C4H10O/c1-4(2,3)5/h5H,1-3H3",
            "inchikey": "DKGRVHCAPJKBHC-UHFFFAOYSA-N",
            "formula": "C4H10O",
            "properties": {
                "logP": 0.35,
                "h_bond_donors": 1,
                "h_bond_acceptors": 1,
                "molecular_weight": 74.12
            }
        },
        {
            "chembl_id": "CHEMBL1487",
            "name": "Urea",
            "smiles": "C(=O)(N)N",
            "inchi": "InChI=1S/CH4N2O/c2-1(3)4/h(H4,2,3,4)",
            "inchikey": "XSQUKJJJFZCRTK-UHFFFAOYSA-N",
            "formula": "CH4N2O",
            "properties": {
                "logP": -2.11,
                "h_bond_donors": 2,
                "h_bond_acceptors": 1,
                "molecular_weight": 60.06
            }
        },
        {
            "chembl_id": "CHEMBL6196",
            "name": "Ethylene glycol",
            "smiles": "C(CO)O",
            "inchi": "InChI=1S/C2H6O2/c3-1-2-4/h3-4H,1-2H2",
            "inchikey": "LYCAIKOWRPUZTN-UHFFFAOYSA-N",
            "formula": "C2H6O2",
            "properties": {
                "logP": -1.36,
                "h_bond_donors": 2,
                "h_bond_acceptors": 2,
                "molecular_weight": 62.07
            }
        },
        {
            "chembl_id": "CHEMBL967",
            "name": "Propylene glycol",
            "smiles": "CC(CO)O",
            "inchi": "InChI=1S/C3H8O2/c1-3(5)2-4/h3-5H,2H2,1H3",
            "inchikey": "DNIAPMSPPWPWGF-UHFFFAOYSA-N",
            "formula": "C3H8O2",
            "properties": {
                "logP": -0.92,
                "h_bond_donors": 2,
                "h_bond_acceptors": 2,
                "molecular_weight": 76.09
            }
        },
        {
            "chembl_id": "CHEMBL262548",
            "name": "Trehalose",
            "smiles": "C(C1C(C(C(C(O1)OC2C(C(C(C(O2)CO)O)O)O)O)O)O)O",
            "inchi": "InChI=1S/C12H22O11/c13-1-4-7(16)8(17)9(18)11(21-4)23-12-10(19)6(15)5(14)3(2-13)22-12/h3-19H,1-2H2/t3-,4-,5-,6-,7-,8-,9-,10-,11-,12-/m1/s1",
            "inchikey": "LFQSCWFLJHTTHZ-LUBWMEBZSA-N",
            "formula": "C12H22O11",
            "properties": {
                "logP": -3.77,
                "h_bond_donors": 8,
                "h_bond_acceptors": 11,
                "molecular_weight": 342.30
            }
        },
        {
            "chembl_id": "CHEMBL6752",
            "name": "Glycine",
            "smiles": "C(C(=O)O)N",
            "inchi": "InChI=1S/C2H5NO2/c3-1-2(4)5/h1,3H2,(H,4,5)",
            "inchikey": "DHMQDGOQFOQNFH-UHFFFAOYSA-N",
            "formula": "C2H5NO2",
            "properties": {
                "logP": -3.21,
                "h_bond_donors": 1,
                "h_bond_acceptors": 3,
                "molecular_weight": 75.07
            }
        }
    ]
    
    # Use the simplest possible insertion method
    success_count = 0
    for compound in reference_compounds:
        try:
            # 1. Insert basic molecule
            molecule_id = execute_query("""
                INSERT INTO molecules (name, chembl_id, smiles, inchi, inchikey, formula, data_source)
                VALUES (%(name)s, %(chembl_id)s, %(smiles)s, %(inchi)s, %(inchikey)s, %(formula)s, 'reference')
                ON CONFLICT (chembl_id) DO UPDATE SET
                    name = EXCLUDED.name,
                    smiles = EXCLUDED.smiles,
                    inchi = EXCLUDED.inchi,
                    inchikey = EXCLUDED.inchikey,
                    formula = EXCLUDED.formula,
                    updated_at = NOW()
                RETURNING id
            """, compound, fetch_one=True)['id']
            
            # 2. Insert each property one by one
            for prop_name, value in compound['properties'].items():
                # Get or create property type
                prop_type_id = execute_query("""
                    INSERT INTO property_types (name, data_type)
                    VALUES (%(name)s, %(data_type)s)
                    ON CONFLICT (name) DO UPDATE SET
                        updated_at = NOW()
                    RETURNING id
                """, {
                    'name': prop_name,
                    'data_type': 'numeric'
                }, fetch_one=True)['id']
                
                # Insert property value
                execute_query("""
                    INSERT INTO molecular_properties (molecule_id, property_type_id, numeric_value)
                    VALUES (%(molecule_id)s, %(property_type_id)s, %(value)s)
                    ON CONFLICT (molecule_id, property_type_id) DO UPDATE SET
                        numeric_value = %(value)s,
                        updated_at = NOW()
                """, {
                    'molecule_id': molecule_id,
                    'property_type_id': prop_type_id,
                    'value': value
                })
            
            # If we get here, all properties were inserted successfully
            success_count += 1
            print(f"Successfully populated reference compound: {compound['name']}")
            
        except Exception as e:
            print(f"Error populating reference compound {compound['name']}: {str(e)}")
    
    print(f"Successfully populated {success_count}/9 reference compounds")
    return success_count
```

### Step 2: Verify Reference Compounds Only
```python
def verify_reference_compounds_only():
    """
    Verify ONLY that reference compounds exist with required properties.
    This is our simplest verification step.
    """
    # Get all reference compounds
    ref_compounds = execute_query("""
        SELECT id, name, chembl_id FROM molecules
        WHERE chembl_id IN (
            'CHEMBL388978', 'CHEMBL1098659', 'CHEMBL66195', 'CHEMBL500033',
            'CHEMBL1487', 'CHEMBL6196', 'CHEMBL967', 'CHEMBL262548', 'CHEMBL6752'
        )
    """)
    
    found_count = len(ref_compounds)
    print(f"Found {found_count}/9 reference compounds")
    
    # For each found compound, check critical properties
    complete_count = 0
    for compound in ref_compounds:
        required_props = ['logP', 'h_bond_donors', 'h_bond_acceptors']
        missing_props = []
        
        for prop in required_props:
            # Check if property exists
            result = execute_query("""
                SELECT mp.numeric_value
                FROM molecular_properties mp
                JOIN property_types pt ON mp.property_type_id = pt.id
                WHERE mp.molecule_id = %(molecule_id)s AND pt.name = %(prop_name)s
            """, {
                'molecule_id': compound['id'],
                'prop_name': prop
            }, fetch_one=True)
            
            if not result or result['numeric_value'] is None:
                missing_props.append(prop)
        
        if not missing_props:
            complete_count += 1
            print(f"✓ {compound['name']} has all required properties")
        else:
            print(f"✗ {compound['name']} missing properties: {', '.join(missing_props)}")
    
    print(f"Found {complete_count}/{found_count} complete reference compounds")
    
    # Return simple success/failure
    return complete_count == 9
```

## Phase 2: Expand to PubChem Data

### Step 1: Simplified PubChem Import
```python
def import_pubchem_minimal(limit=500):
    """
    Import a minimal set of molecules from PubChem with critical properties.
    Focus on reliability rather than quantity.
    """
    # Get CID list from file
    with open("CID-Synonym-curated", "r") as f:
        cids = [int(line.strip().split("\t")[0]) for line in f 
                if line.strip() and line.strip().split("\t")[0].isdigit()][:limit]
    
    print(f"Starting minimal import of {len(cids)} PubChem compounds")
    
    # Process in small batches of 10
    batch_size = 10
    success_count = 0
    
    for i in range(0, len(cids), batch_size):
        batch = cids[i:i+batch_size]
        print(f"Processing batch {i//batch_size + 1}/{(len(cids) + batch_size - 1)//batch_size}")
        
        # Process each compound in batch
        batch_success = 0
        for cid in batch:
            try:
                # 1. Fetch compound data from PubChem
                data = fetch_pubchem_data(cid)
                if not data:
                    print(f"No data found for CID {cid}")
                    continue
                
                # 2. Insert molecule
                molecule_id = insert_pubchem_molecule(data)
                if not molecule_id:
                    continue
                
                # 3. Insert properties
                property_success = insert_pubchem_properties(molecule_id, data)
                if property_success:
                    batch_success += 1
                    
            except Exception as e:
                print(f"Error processing CID {cid}: {str(e)}")
        
        # Update success count and log progress
        success_count += batch_success
        print(f"Batch complete. Total successful imports: {success_count}/{i+len(batch)}")
        
        # Simple file-based checkpoint
        with open("pubchem_import_checkpoint.txt", "w") as f:
            f.write(f"{i+batch_size}\n{success_count}")
            
        # Small delay to avoid overwhelming API
        time.sleep(1)
    
    print(f"Completed PubChem import: {success_count}/{len(cids)} compounds")
    return success_count

def fetch_pubchem_data(cid):
    """
    Fetch compound data from PubChem with simple retry mechanism.
    Focus on the most crucial properties for cryoprotectants.
    """
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/"
        "MolecularFormula,MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,"
        "IsomericSMILES,InChI,InChIKey,IUPACName,Title/JSON"
    )
    
    for attempt in range(3):  # Simple retry mechanism
        try:
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                data = response.json()
                return data["PropertyTable"]["Properties"][0]
            
            # Handle rate limiting
            if response.status_code == 429:
                wait_time = (attempt + 1) * 2  # Simple exponential backoff
                print(f"Rate limited. Waiting {wait_time} seconds...")
                time.sleep(wait_time)
                continue
                
            print(f"Failed to fetch CID {cid}: HTTP {response.status_code}")
            return None
            
        except Exception as e:
            print(f"Error fetching CID {cid}, attempt {attempt+1}: {str(e)}")
            time.sleep(1)
    
    return None

def insert_pubchem_molecule(data):
    """
    Insert molecule data from PubChem into database.
    Simple, straightforward implementation with clear logging.
    """
    try:
        # Check for required fields
        if not all(key in data for key in ["IsomericSMILES", "InChI", "InChIKey", "MolecularFormula"]):
            print(f"Missing required fields for CID {data['CID']}")
            return None
        
        # Prepare data
        molecule_data = {
            "pubchem_cid": data["CID"],
            "name": data.get("Title", f"CID {data['CID']}"),
            "smiles": data["IsomericSMILES"],
            "inchi": data["InChI"],
            "inchikey": data["InChIKey"],
            "formula": data["MolecularFormula"],
            "molecular_weight": data.get("MolecularWeight"),
            "data_source": "PubChem"
        }
        
        # Execute insert
        result = execute_query("""
            INSERT INTO molecules 
                (pubchem_cid, name, smiles, inchi, inchikey, formula, molecular_weight, data_source)
            VALUES 
                (%(pubchem_cid)s, %(name)s, %(smiles)s, %(inchi)s, %(inchikey)s, 
                 %(formula)s, %(molecular_weight)s, %(data_source)s)
            ON CONFLICT (pubchem_cid) DO UPDATE SET 
                name = EXCLUDED.name,
                smiles = EXCLUDED.smiles,
                inchi = EXCLUDED.inchi, 
                inchikey = EXCLUDED.inchikey,
                formula = EXCLUDED.formula,
                molecular_weight = EXCLUDED.molecular_weight,
                updated_at = NOW()
            RETURNING id
        """, molecule_data, fetch_one=True)
        
        if result and 'id' in result:
            print(f"Successfully inserted molecule for CID {data['CID']}")
            return result['id']
        else:
            print(f"Failed to insert molecule for CID {data['CID']}")
            return None
            
    except Exception as e:
        print(f"Error inserting molecule for CID {data['CID']}: {str(e)}")
        return None

def insert_pubchem_properties(molecule_id, data):
    """
    Insert properties for a PubChem molecule.
    Focus on getting the 3 critical properties correctly.
    """
    try:
        # Extract the 3 critical properties
        properties = {}
        
        # These are the most important properties for our verification criteria
        if "XLogP" in data:
            try:
                properties["logP"] = float(data["XLogP"])
            except (ValueError, TypeError):
                print(f"Invalid XLogP value: {data['XLogP']}")
        
        if "HBondDonorCount" in data:
            try:
                properties["h_bond_donors"] = int(data["HBondDonorCount"])
            except (ValueError, TypeError):
                print(f"Invalid HBondDonorCount value: {data['HBondDonorCount']}")
        
        if "HBondAcceptorCount" in data:
            try:
                properties["h_bond_acceptors"] = int(data["HBondAcceptorCount"])
            except (ValueError, TypeError):
                print(f"Invalid HBondAcceptorCount value: {data['HBondAcceptorCount']}")
        
        # Check if we have all required properties
        if not all(prop in properties for prop in ["logP", "h_bond_donors", "h_bond_acceptors"]):
            print(f"Missing critical properties for CID {data['CID']}")
            return False
        
        # Insert properties one by one for clarity
        for prop_name, value in properties.items():
            # Get or create property type
            prop_type_id = execute_query("""
                INSERT INTO property_types (name, data_type)
                VALUES (%(name)s, %(data_type)s)
                ON CONFLICT (name) DO UPDATE SET
                    updated_at = NOW()
                RETURNING id
            """, {
                'name': prop_name,
                'data_type': 'numeric'
            }, fetch_one=True)['id']
            
            # Insert property value
            execute_query("""
                INSERT INTO molecular_properties (molecule_id, property_type_id, numeric_value)
                VALUES (%(molecule_id)s, %(property_type_id)s, %(value)s)
                ON CONFLICT (molecule_id, property_type_id) DO UPDATE SET
                    numeric_value = %(value)s,
                    updated_at = NOW()
            """, {
                'molecule_id': molecule_id,
                'property_type_id': prop_type_id,
                'value': value
            })
        
        print(f"Successfully inserted properties for molecule {molecule_id}")
        return True
        
    except Exception as e:
        print(f"Error inserting properties for molecule {molecule_id}: {str(e)}")
        return False
```

### Step 2: Verify PubChem Molecules
```python
def verify_pubchem_molecules():
    """
    Verify that PubChem molecules exist with required properties.
    Focus on property quality rather than pure quantity.
    """
    # Count PubChem molecules
    pubchem_count = execute_query("""
        SELECT COUNT(*) as count FROM molecules 
        WHERE pubchem_cid IS NOT NULL
    """, fetch_one=True)['count']
    
    # Count PubChem molecules with all critical properties
    complete_count = execute_query("""
        SELECT COUNT(DISTINCT m.id) as count
        FROM molecules m
        JOIN molecular_properties mp_logp ON m.id = mp_logp.molecule_id
        JOIN property_types pt_logp ON mp_logp.property_type_id = pt_logp.id AND pt_logp.name = 'logP'
        JOIN molecular_properties mp_hbd ON m.id = mp_hbd.molecule_id
        JOIN property_types pt_hbd ON mp_hbd.property_type_id = pt_hbd.id AND pt_hbd.name = 'h_bond_donors'
        JOIN molecular_properties mp_hba ON m.id = mp_hba.molecule_id
        JOIN property_types pt_hba ON mp_hba.property_type_id = pt_hba.id AND pt_hba.name = 'h_bond_acceptors'
        WHERE m.pubchem_cid IS NOT NULL
    """, fetch_one=True)['count']
    
    completion_percent = (complete_count / pubchem_count * 100) if pubchem_count > 0 else 0
    print(f"PubChem molecules: {pubchem_count} total, {complete_count} complete ({completion_percent:.1f}%)")
    
    # Check if we meet requirements
    return pubchem_count >= 500 and completion_percent >= 90
```

## Phase 3: Expand to ChEMBL Data

### Step 1: Simplified ChEMBL Import
```python
def import_chembl_minimal(limit=500):
    """
    Import a minimal set of molecules from ChEMBL with critical properties.
    Focus on reliability rather than quantity.
    """
    # Similar implementation to PubChem import but for ChEMBL
    # Get reference compound ChEMBL IDs as a starting point
    reference_ids = [
        "CHEMBL388978", "CHEMBL1098659", "CHEMBL66195", "CHEMBL500033",
        "CHEMBL1487", "CHEMBL6196", "CHEMBL967", "CHEMBL262548", "CHEMBL6752"
    ]
    
    # Find similar compounds (simplified approach)
    all_compounds = []
    
    # First add our reference IDs
    for ref_id in reference_ids:
        all_compounds.append(ref_id)
    
    # Then find similar compounds by API call
    # (Simplified implementation)
    
    # Deduplicate
    all_compounds = list(set(all_compounds))[:limit]
    
    print(f"Starting import of {len(all_compounds)} ChEMBL compounds")
    
    # Process each compound
    success_count = 0
    for i, chembl_id in enumerate(all_compounds):
        try:
            print(f"Processing compound {i+1}/{len(all_compounds)}: {chembl_id}")
            
            # 1. Fetch compound data from ChEMBL
            data = fetch_chembl_data(chembl_id)
            if not data:
                continue
            
            # 2. Insert molecule
            molecule_id = insert_chembl_molecule(data)
            if not molecule_id:
                continue
            
            # 3. Insert properties
            if insert_chembl_properties(molecule_id, data):
                success_count += 1
                
            # Simple checkpoint every 10 compounds
            if (i + 1) % 10 == 0:
                print(f"Progress: {i+1}/{len(all_compounds)} compounds processed, {success_count} successful")
                
        except Exception as e:
            print(f"Error processing ChEMBL ID {chembl_id}: {str(e)}")
    
    print(f"Completed ChEMBL import: {success_count}/{len(all_compounds)} compounds")
    return success_count
```

### Step 2: Verify Combined Data
```python
def verify_all_molecules():
    """
    Verify all molecules against requirements.
    """
    # Count all molecules
    total_count = execute_query("""
        SELECT COUNT(*) as count FROM molecules
    """, fetch_one=True)['count']
    
    # Count molecules with required properties
    complete_count = execute_query("""
        SELECT COUNT(DISTINCT m.id) as count
        FROM molecules m
        JOIN molecular_properties mp_logp ON m.id = mp_logp.molecule_id
        JOIN property_types pt_logp ON mp_logp.property_type_id = pt_logp.id AND pt_logp.name = 'logP'
        JOIN molecular_properties mp_hbd ON m.id = mp_hbd.molecule_id
        JOIN property_types pt_hbd ON mp_hbd.property_type_id = pt_hbd.id AND pt_hbd.name = 'h_bond_donors'
        JOIN molecular_properties mp_hba ON m.id = mp_hba.molecule_id
        JOIN property_types pt_hba ON mp_hba.property_type_id = pt_hba.id AND pt_hba.name = 'h_bond_acceptors'
    """, fetch_one=True)['count']
    
    # Calculate completion percentage
    completion_percent = (complete_count / total_count * 100) if total_count > 0 else 0
    
    # Print results
    print(f"Total molecules: {total_count}")
    print(f"Complete molecules: {complete_count} ({completion_percent:.1f}%)")
    
    # Check all requirements
    ref_complete = verify_reference_compounds_only()
    pubchem_complete = verify_pubchem_molecules()
    target_met = total_count >= 5000
    completion_met = completion_percent >= 90
    
    print(f"Reference compounds complete: {ref_complete}")
    print(f"PubChem verification: {pubchem_complete}")
    print(f"Target molecule count (5000+): {target_met}")
    print(f"Property completion target (90%+): {completion_met}")
    
    # All requirements met?
    return ref_complete and pubchem_complete and target_met and completion_met
```

## Phase 4: Only Then Optimize Performance

```python
def add_basic_indexes():
    """
    Add simple indexes to improve query performance.
    Only run after basic functionality is verified.
    """
    execute_query("""
        -- Add indexes for commonly queried fields
        CREATE INDEX IF NOT EXISTS idx_molecules_pubchem_cid ON molecules (pubchem_cid);
        CREATE INDEX IF NOT EXISTS idx_molecules_chembl_id ON molecules (chembl_id);
        CREATE INDEX IF NOT EXISTS idx_molecules_inchikey ON molecules (inchikey);
        CREATE INDEX IF NOT EXISTS idx_molecular_properties_molecule_id ON molecular_properties (molecule_id);
        CREATE INDEX IF NOT EXISTS idx_molecular_properties_property_type_id ON molecular_properties (property_type_id);
        CREATE INDEX IF NOT EXISTS idx_molecular_properties_numeric_value ON molecular_properties (numeric_value);
        CREATE INDEX IF NOT EXISTS idx_property_types_name ON property_types (name);
        
        -- Add compound index for frequent join
        CREATE INDEX IF NOT EXISTS idx_molecular_properties_molecule_property ON molecular_properties (molecule_id, property_type_id);
    """)
    
    # Run ANALYZE to update statistics
    execute_query("ANALYZE")
    
    print("Added basic performance indexes")
```

## Utility Functions

```python
def execute_query(query, params=None, fetch_one=False, return_data=True):
    """
    Simple, consistent query execution function.
    Abstracts away database connection details.
    """
    try:
        # Get connection
        from database.connection import get_db_connection
        conn = get_db_connection()
        
        # Execute query
        cursor = conn.cursor()
        cursor.execute(query, params)
        
        # Return results if needed
        if return_data:
            if fetch_one:
                return cursor.fetchone()
            else:
                return cursor.fetchall()
        
        # Commit if not a SELECT query
        if not query.strip().upper().startswith("SELECT"):
            conn.commit()
            
        return True
        
    except Exception as e:
        print(f"Error executing query: {str(e)}")
        print(f"Query: {query}")
        print(f"Params: {params}")
        raise
```