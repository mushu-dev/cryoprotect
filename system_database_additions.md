# Molecular Filtering and Search Mechanisms

## Molecular Selection and Filtering

CryoProtect implements sophisticated molecular filtering mechanisms to ensure scientific relevance of compounds for cryoprotectant research. This section documents the filtering approaches, search algorithms, and data flow patterns used throughout the system.

### Filtering Criteria for Cryoprotectant Relevance

Molecules are filtered using multiple scientifically-validated criteria to identify potentially useful cryoprotectants:

1. **Physicochemical Property Filters**:
   - **Molecular Weight**: 32-500 Da (optimum range for cryoprotectant function)
   - **logP (Octanol-water partition coefficient)**: -3.0 to 1.5 (water solubility with moderate membrane permeability)
   - **Hydrogen Bond Donors**: 1-6 (ability to form stabilizing interactions with biomolecules)
   - **Hydrogen Bond Acceptors**: 2-10 (ability to disrupt water crystallization)
   - **Rotatable Bonds**: 0-8 (conformational flexibility for interacting with water)
   - **Polar Surface Area**: 20-140 Å² (balance between hydrophilicity and cell permeability)

2. **Functional Group Filters**:
   - **Hydroxyl Groups (-OH)**: Minimum of 1 (essential for hydrogen bonding with water)
   - **Amine Groups (-NH₂)**: 0-3 (may enhance interaction with biological macromolecules)
   - **Toxic or Reactive Group Exclusion**: Elimination of molecules containing reactive functionalities (e.g., acyl halides, epoxides, aldehydes) or known toxic groups

3. **Structural Class Filters**:
   - **Polyols**: Multi-hydroxyl compounds like glycerol, propylene glycol
   - **Sugars and Sugar Alcohols**: Glucose, mannitol, trehalose
   - **Amides**: Formamide, acetamide
   - **Sulfoxides**: DMSO and similar compounds
   - **Small Polymers**: PEG variants with appropriate molecular weights

### Property-Based Filtering Implementation

Property-based filtering is implemented via a configuration-driven system that translates scientific criteria into SQL queries or in-memory filters:

```json
{
  "cryoprotectant_filter_profile": {
    "name": "standard_cryoprotectant_profile",
    "description": "Standard filtering criteria for potential cryoprotectants",
    "version": "2.0",
    "numeric_filters": [
      {
        "property": "molecular_weight",
        "min": 32.0,
        "max": 500.0,
        "unit": "Da",
        "importance": "critical"
      },
      {
        "property": "logP",
        "min": -3.0,
        "max": 1.5,
        "unit": "",
        "importance": "high"
      },
      {
        "property": "hydrogen_bond_donors",
        "min": 1,
        "max": 6,
        "unit": "count",
        "importance": "high"
      },
      {
        "property": "hydrogen_bond_acceptors",
        "min": 2,
        "max": 10,
        "unit": "count",
        "importance": "high"
      },
      {
        "property": "rotatable_bonds",
        "min": 0,
        "max": 8,
        "unit": "count",
        "importance": "medium"
      },
      {
        "property": "polar_surface_area",
        "min": 20.0,
        "max": 140.0,
        "unit": "Å²",
        "importance": "medium"
      }
    ],
    "structural_filters": [
      {
        "type": "functional_group",
        "pattern": "hydroxyl",
        "smarts": "[OX2H]",
        "min_count": 1,
        "importance": "critical"
      },
      {
        "type": "functional_group",
        "pattern": "amine_primary",
        "smarts": "[NX3;H2]",
        "min_count": 0,
        "max_count": 3,
        "importance": "medium"
      },
      {
        "type": "exclusion",
        "pattern": "acyl_halide",
        "smarts": "[CX3](=[OX1])[F,Cl,Br,I]",
        "importance": "critical"
      }
    ],
    "class_filters": [
      {
        "class": "polyol",
        "description": "Multiple hydroxyl functional groups",
        "smarts": "[OX2H][CX4][CX4][OX2H]",
        "score_bonus": 1.5
      },
      {
        "class": "sugar",
        "description": "Monosaccharide or disaccharide structures",
        "smarts": "[OX2H]C1[OX2][CX4][CX4][CX4]1",
        "score_bonus": 1.3
      }
    ]
  }
}
```

## Structure-Based Similarity Search

CryoProtect implements multiple molecular similarity algorithms to identify structurally related compounds and evaluate potential cryoprotectant efficacy.

### Fingerprint-Based Similarity Methods

1. **Morgan Fingerprints (ECFP)**:
   - Implementation: RDKit's Morgan fingerprints with radius 2 (ECFP4 equivalent)
   - Bit length: 2048 bits
   - Application: Fast screening of large databases for structural similarity
   - SQL Implementation:

```sql
CREATE OR REPLACE FUNCTION morganbv_similarity(smiles1 TEXT, smiles2 TEXT)
RETURNS FLOAT AS $$
  from rdkit import Chem
  from rdkit.Chem import AllChem
  import rdkit.DataStructs

  # Generate molecules from SMILES
  mol1 = Chem.MolFromSmiles(smiles1)
  mol2 = Chem.MolFromSmiles(smiles2)
  
  if mol1 is None or mol2 is None:
    return 0.0
    
  # Generate Morgan fingerprints (ECFP4)
  fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
  fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
  
  # Calculate Tanimoto similarity
  similarity = rdkit.DataStructs.TanimotoSimilarity(fp1, fp2)
  return similarity
$$ LANGUAGE plpython3u;
```

2. **Topological Fingerprints**:
   - Implementation: RDKit's RDKFingerprint 
   - Bit length: 2048 bits
   - Application: Complementary to Morgan fingerprints, captures different structural features

3. **MACCS Keys**:
   - Implementation: 166-bit structural keys based on MACCS definitions
   - Application: Finds molecules with similar functional group patterns
   - Performance: Lightweight, good for rapid screening

### 3D Similarity Methods

For compounds with experimental or predicted 3D structures:

1. **Shape-Based Similarity**:
   - Implementation: RDKit's Shape Protocols with Gaussian overlay
   - Features: Conformer generation with ETKDG algorithm
   - Application: Identifies molecules with similar 3D shapes that may have similar binding interactions

2. **Pharmacophore-Based Similarity**:
   - Implementation: Custom 3D pharmacophore models representing cryoprotectant features
   - Features: Maps hydrogen-bond donors/acceptors and hydrophobic regions
   - Application: Identifies molecules with potentially similar interaction patterns with water and biological macromolecules

### Database Implementation of Structural Similarity Searching

Structural similarity is implemented via custom PostgreSQL functions and indexes:

```sql
-- Create RDKit extension and necessary functions
CREATE EXTENSION IF NOT EXISTS rdkit;

-- Create index for substructure and similarity searching
CREATE INDEX idx_molecules_rdkit_mol
ON molecules
USING gist(mol);

-- Similarity search function using Morgan fingerprints
CREATE OR REPLACE FUNCTION similar_molecules(query_smiles TEXT, similarity_threshold FLOAT DEFAULT 0.7, limit_count INT DEFAULT 100)
RETURNS TABLE (
  id UUID,
  name TEXT,
  smiles TEXT,
  similarity FLOAT
) AS $$
BEGIN
  RETURN QUERY
  SELECT 
    m.id,
    m.name,
    m.smiles,
    morganbv_similarity(m.smiles, query_smiles) AS similarity
  FROM 
    molecules m
  WHERE 
    morganbv_similarity(m.smiles, query_smiles) >= similarity_threshold
  ORDER BY 
    similarity DESC
  LIMIT limit_count;
END;
$$ LANGUAGE plpgsql;
```

## Property-Based Search Mechanisms

### Primary Cryoprotectant Properties

CryoProtect allows searching and filtering based on key properties relevant to cryoprotection:

1. **Glass Transition Temperature (Tg)**:
   - Range filtering: e.g., Tg > -40°C for practical vitrification
   - Database table: `molecular_properties` with property_type 'glass_transition_temperature'

2. **Cytotoxicity Metrics**:
   - Range filtering: LC50, IC50 values above thresholds for safety
   - Database table: `toxicity_data` with specific cell lines and conditions

3. **Ice Crystal Inhibition**:
   - Range filtering: Quantitative measurements of ice crystal growth inhibition
   - Database table: `molecular_properties` with property_type 'ice_inhibition_factor'

4. **Cell Membrane Permeability**:
   - Range filtering: Permeability coefficients for cell penetration
   - Database table: `molecular_properties` with property_type 'membrane_permeability'

5. **Viscosity Impact**:
   - Range filtering: Solution viscosity at various concentrations
   - Database table: `molecular_properties` with property_type 'viscosity_factor'

### Advanced Property Search Implementation

Advanced property-based searching is implemented with complex SQL queries using parameterized filters:

```sql
-- Example query for retrieving molecules with specific cryoprotectant properties
CREATE OR REPLACE FUNCTION find_cryoprotectant_candidates(
  min_tg FLOAT DEFAULT NULL,
  max_toxicity FLOAT DEFAULT NULL,
  min_permeability FLOAT DEFAULT NULL,
  min_ice_inhibition FLOAT DEFAULT NULL
)
RETURNS TABLE (
  id UUID,
  name TEXT,
  smiles TEXT,
  tg FLOAT,
  toxicity FLOAT,
  permeability FLOAT,
  ice_inhibition FLOAT,
  score FLOAT
) AS $$
BEGIN
  RETURN QUERY
  WITH 
    tg_data AS (
      SELECT molecule_id, numeric_value AS tg_value
      FROM molecular_properties
      WHERE property_type_id = (SELECT id FROM property_types WHERE name = 'glass_transition_temperature')
      AND (min_tg IS NULL OR numeric_value >= min_tg)
    ),
    toxicity_data AS (
      SELECT molecule_id, numeric_value AS toxicity_value
      FROM molecular_properties
      WHERE property_type_id = (SELECT id FROM property_types WHERE name = 'cytotoxicity_lc50')
      AND (max_toxicity IS NULL OR numeric_value <= max_toxicity)
    ),
    permeability_data AS (
      SELECT molecule_id, numeric_value AS permeability_value
      FROM molecular_properties
      WHERE property_type_id = (SELECT id FROM property_types WHERE name = 'membrane_permeability')
      AND (min_permeability IS NULL OR numeric_value >= min_permeability)
    ),
    ice_inhibition_data AS (
      SELECT molecule_id, numeric_value AS ice_inhibition_value
      FROM molecular_properties
      WHERE property_type_id = (SELECT id FROM property_types WHERE name = 'ice_inhibition_factor')
      AND (min_ice_inhibition IS NULL OR numeric_value >= min_ice_inhibition)
    )
  SELECT 
    m.id,
    m.name,
    m.smiles,
    tg.tg_value AS tg,
    tox.toxicity_value AS toxicity,
    perm.permeability_value AS permeability,
    ice.ice_inhibition_value AS ice_inhibition,
    -- Calculate weighted score for ranking
    (COALESCE(tg.tg_value / 100, 0) * 0.3 + 
     COALESCE((1000 - tox.toxicity_value) / 1000, 0) * 0.3 + 
     COALESCE(perm.permeability_value / 10, 0) * 0.2 + 
     COALESCE(ice.ice_inhibition_value / 10, 0) * 0.2) AS score
  FROM 
    molecules m
  LEFT JOIN tg_data tg ON m.id = tg.molecule_id
  LEFT JOIN toxicity_data tox ON m.id = tox.molecule_id
  LEFT JOIN permeability_data perm ON m.id = perm.molecule_id
  LEFT JOIN ice_inhibition_data ice ON m.id = ice.molecule_id
  WHERE
    (min_tg IS NULL OR tg.molecule_id IS NOT NULL) AND
    (max_toxicity IS NULL OR tox.molecule_id IS NOT NULL) AND
    (min_permeability IS NULL OR perm.molecule_id IS NOT NULL) AND
    (min_ice_inhibition IS NULL OR ice.molecule_id IS NOT NULL)
  ORDER BY 
    score DESC;
END;
$$ LANGUAGE plpgsql;
```

## Combined Search Strategies

### Multi-criteria Search Implementation

CryoProtect implements combined property and structure-based searches for maximum scientific relevance:

```json
{
  "search_criteria": {
    "structural_similarity": {
      "reference_molecule": "DMSO",
      "reference_smiles": "CS(C)=O",
      "similarity_method": "morgan",
      "similarity_threshold": 0.6,
      "max_results": 100
    },
    "property_filters": [
      {
        "property": "molecular_weight",
        "operator": "between",
        "min_value": 50.0,
        "max_value": 250.0,
        "unit": "Da"
      },
      {
        "property": "glass_transition_temperature",
        "operator": "greater_than",
        "value": -45.0,
        "unit": "°C"
      },
      {
        "property": "cytotoxicity_lc50",
        "operator": "greater_than",
        "value": 500.0,
        "unit": "μM"
      }
    ],
    "scoring_weights": {
      "structural_similarity": 0.4,
      "glass_transition_temperature": 0.3,
      "cytotoxicity": 0.2,
      "membrane_permeability": 0.1
    }
  }
}
```

### Data Pipeline for Cryoprotectant Screening

The complete data pipeline combines multiple filtering and search strategies:

1. **Initial Broad Screening**:
   - Property-based filters applied to the complete molecular database
   - Exclusion of obviously unsuitable compounds (reactive, toxic groups, etc.)
   - Output: Subset of molecules meeting basic criteria (~10-20% of database)

2. **Structural Similarity Enrichment**:
   - Similarity search using known effective cryoprotectants as references
   - Multiple fingerprint methods applied in parallel
   - Output: Ranked list of structurally relevant compounds

3. **Advanced Property Filtering**:
   - Application of more stringent experimental property filters
   - Integration of toxicity and cytotoxicity data when available
   - Output: Highly filtered candidate list (~1-2% of database)

4. **Composite Scoring**:
   - Weighted scoring combining structural similarity and property data
   - Customizable weighting schemes for different preservation protocols
   - Output: Final ranked list with composite scores (0-100)

5. **Experimental Validation Pipeline**:
   - Recommendations for laboratory testing based on score
   - Tracking of experimental results to refine future screening
   - Output: Validation data fed back into the system

### Example Data Pipeline Implementation

```python
def cryoprotectant_screening_pipeline(db_connection, config):
    """
    Complete pipeline for cryoprotectant screening and selection.
    
    Parameters:
    -----------
    db_connection : psycopg2.connection
        Database connection object
    config : dict
        Pipeline configuration
        
    Returns:
    --------
    dict
        Results of screening with candidates and scores
    """
    # Step 1: Initial property-based filtering
    with db_connection.cursor() as cursor:
        cursor.execute("""
            SELECT id, name, smiles, molecular_weight 
            FROM molecules
            WHERE 
                molecular_weight BETWEEN %(min_mw)s AND %(max_mw)s
                AND EXISTS (
                    SELECT 1 FROM molecular_properties 
                    WHERE molecule_id = molecules.id 
                    AND property_type_id = (SELECT id FROM property_types WHERE name = 'logP')
                    AND numeric_value BETWEEN %(min_logp)s AND %(max_logp)s
                )
        """, {
            'min_mw': config['property_filters']['molecular_weight']['min'],
            'max_mw': config['property_filters']['molecular_weight']['max'],
            'min_logp': config['property_filters']['logP']['min'],
            'max_logp': config['property_filters']['logP']['max'],
        })
        initial_candidates = cursor.fetchall()
    
    # Step 2: Structural filtering using SMARTS patterns
    from rdkit import Chem
    from rdkit.Chem import AllChem
    
    filtered_candidates = []
    required_patterns = [Chem.MolFromSmarts(pattern) for pattern in config['structural_filters']['required']]
    exclusion_patterns = [Chem.MolFromSmarts(pattern) for pattern in config['structural_filters']['exclusion']]
    
    for candidate in initial_candidates:
        mol = Chem.MolFromSmiles(candidate[2])  # SMILES is at index 2
        if mol is None:
            continue
            
        # Check required patterns
        if all(mol.HasSubstructMatch(pattern) for pattern in required_patterns):
            # Check exclusion patterns
            if not any(mol.HasSubstructMatch(pattern) for pattern in exclusion_patterns):
                filtered_candidates.append(candidate)
    
    # Step 3: Similarity search enrichment
    reference_mols = []
    for ref in config['reference_molecules']:
        ref_mol = Chem.MolFromSmiles(ref['smiles'])
        if ref_mol:
            reference_mols.append((ref_mol, ref['weight']))
    
    similarity_scores = {}
    for candidate in filtered_candidates:
        mol = Chem.MolFromSmiles(candidate[2])
        if mol is None:
            continue
            
        # Calculate Morgan fingerprint for candidate
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        
        # Calculate weighted similarity to reference molecules
        total_similarity = 0
        total_weight = 0
        
        for ref_mol, weight in reference_mols:
            ref_fp = AllChem.GetMorganFingerprintAsBitVect(ref_mol, 2, nBits=2048)
            similarity = DataStructs.TanimotoSimilarity(fp, ref_fp)
            total_similarity += similarity * weight
            total_weight += weight
        
        # Store weighted average similarity
        if total_weight > 0:
            similarity_scores[candidate[0]] = total_similarity / total_weight
    
    # Step 4: Advanced property scoring
    final_scores = {}
    for candidate in filtered_candidates:
        # Skip if didn't pass similarity threshold
        if candidate[0] not in similarity_scores or similarity_scores[candidate[0]] < config['min_similarity']:
            continue
            
        # Query advanced properties
        with db_connection.cursor() as cursor:
            cursor.execute("""
                SELECT 
                    pt.name,
                    mp.numeric_value,
                    pt.units
                FROM 
                    molecular_properties mp
                JOIN 
                    property_types pt ON mp.property_type_id = pt.id
                WHERE 
                    mp.molecule_id = %(molecule_id)s
                    AND pt.name IN %(property_names)s
            """, {
                'molecule_id': candidate[0],
                'property_names': tuple(config['advanced_properties'])
            })
            properties = {row[0]: row[1] for row in cursor.fetchall()}
        
        # Calculate composite score
        score = similarity_scores[candidate[0]] * config['weights']['similarity']
        
        # Add property contributions to score
        for prop_name, weight in config['weights'].items():
            if prop_name == 'similarity':
                continue
                
            if prop_name in properties:
                # Normalize property value to 0-1 range based on min/max in config
                prop_min = config['property_ranges'][prop_name]['min']
                prop_max = config['property_ranges'][prop_name]['max']
                prop_val = properties[prop_name]
                
                # For some properties, higher is better; for others, lower is better
                if config['property_ranges'][prop_name].get('higher_is_better', True):
                    normalized = (prop_val - prop_min) / (prop_max - prop_min)
                else:
                    normalized = (prop_max - prop_val) / (prop_max - prop_min)
                
                normalized = max(0, min(1, normalized))  # Clamp to 0-1
                score += normalized * weight
        
        # Store final score
        final_scores[candidate[0]] = score
    
    # Step 5: Sort and return top candidates
    top_candidates = sorted(
        [(cid, final_scores[cid]) for cid in final_scores],
        key=lambda x: x[1],
        reverse=True
    )[:config['max_results']]
    
    # Retrieve full candidate information
    results = []
    for cid, score in top_candidates:
        with db_connection.cursor() as cursor:
            cursor.execute("""
                SELECT id, name, smiles, molecular_weight, formula
                FROM molecules
                WHERE id = %s
            """, (cid,))
            molecule = cursor.fetchone()
            
            if molecule:
                results.append({
                    'id': molecule[0],
                    'name': molecule[1],
                    'smiles': molecule[2],
                    'molecular_weight': molecule[3],
                    'formula': molecule[4],
                    'similarity_score': similarity_scores[cid],
                    'final_score': score
                })
    
    return {
        'total_initial_candidates': len(initial_candidates),
        'filtered_candidates': len(filtered_candidates),
        'scored_candidates': len(final_scores),
        'top_candidates': results
    }
```

## Future Enhancements for Molecular Filtering

### Machine Learning Integration

1. **Supervised Classification Models**:
   - Models trained on known cryoprotectant efficacy data
   - Features derived from molecular descriptors and experimental properties
   - Implementation: Scikit-learn models with PostgreSQL integration
   - Example:
   
```python
# Example model training code (implemented in separate Python module)
from sklearn.ensemble import RandomForestClassifier
import numpy as np

def train_cryoprotectant_model(training_data):
    # Extract features and labels
    X = np.array([data['features'] for data in training_data])
    y = np.array([data['is_effective_cryoprotectant'] for data in training_data])
    
    # Train model
    model = RandomForestClassifier(n_estimators=100, random_state=42)
    model.fit(X, y)
    
    return model
```

2. **Self-Supervised Structure-Property Models**:
   - Predicts cryoprotectant-relevant properties from structure
   - Implementation: Graph Neural Networks for molecular property prediction
   - Database integration: Stores predictions in molecular_properties table

### Advanced Structural Similarity Methods

1. **Pharmacophore Fingerprints**:
   - Custom fingerprints based on cryoprotectant pharmacophore models
   - Implementation: RDKit's 2D and 3D pharmacophore fingerprints
   - Database function: `pharmacophore_similarity(smiles1, smiles2)`

2. **Maximum Common Substructure Analysis**:
   - Identifies shared substructures between known cryoprotectants
   - Implementation: RDKit's `FindMCS` function
   - Application: Enrichment of screening results with known active scaffolds

### Real-time Feedback Loop

1. **Experimental Results Integration**:
   - Automatic update of molecule scores based on experimental results
   - Implementation: Feedback mechanism from experiment tables to search weights
   - Database trigger: Updates search relevance based on experimental confirmation

2. **Active Learning Pipeline**:
   - Intelligently selects molecules for experimental testing
   - Maximizes information gain with each experiment
   - Implementation: Uncertainty sampling and exploration strategies

## API Examples for Molecular Filtering

### Comprehensive Molecular Search API

The following endpoint provides comprehensive searching and filtering of molecules:

**Endpoint**: `/api/v1/molecules/search`

**Method**: POST

**Request Body**:
```json
{
  "property_filters": [
    {
      "property": "molecular_weight",
      "operator": "between",
      "min_value": 60.0,
      "max_value": 200.0
    },
    {
      "property": "logP",
      "operator": "less_than",
      "value": 1.0
    }
  ],
  "structural_filters": {
    "required_patterns": ["[OX2H]"],
    "excluded_patterns": ["[SX2]"]
  },
  "similarity_search": {
    "reference_smiles": "C(CO)O",
    "method": "morgan",
    "threshold": 0.6
  },
  "pagination": {
    "page": 1,
    "page_size": 20
  },
  "sort_by": "similarity_score",
  "sort_order": "desc"
}
```

**Response**:
```json
{
  "page": 1,
  "page_size": 20,
  "total_count": 156,
  "total_pages": 8,
  "results": [
    {
      "id": "550e8400-e29b-41d4-a716-446655440000",
      "name": "Propylene glycol",
      "smiles": "CC(O)CO",
      "formula": "C3H8O2",
      "molecular_weight": 76.09,
      "properties": {
        "logP": -0.92,
        "hydrogen_bond_donors": 2,
        "hydrogen_bond_acceptors": 2,
        "glass_transition_temperature": -60.0,
        "cytotoxicity_lc50": 1200.0
      },
      "similarity_score": 0.85,
      "relevance_score": 0.78
    },
    {
      "id": "6ba7b810-9dad-11d1-80b4-00c04fd430c8",
      "name": "Ethylene glycol",
      "smiles": "C(CO)O",
      "formula": "C2H6O2",
      "molecular_weight": 62.07,
      "properties": {
        "logP": -1.36,
        "hydrogen_bond_donors": 2,
        "hydrogen_bond_acceptors": 2,
        "glass_transition_temperature": -55.0,
        "cytotoxicity_lc50": 720.0
      },
      "similarity_score": 1.0,
      "relevance_score": 0.75
    }
  ]
}
```