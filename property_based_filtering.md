# Property-Based Molecule Filtering for Cryoprotectants

## Scientific Rationale for Filtering Criteria

Cryoprotectant molecules must possess specific physicochemical properties to effectively protect biological samples during freezing. The filtering criteria used in CryoProtect are based on established scientific principles of cryobiology and molecular biophysics.

### Core Physicochemical Properties

#### 1. Molecular Weight and Size Parameters

```json
{
  "molecular_weight_filter": {
    "description": "Ideal size range for balancing permeability and interaction strength",
    "critical_range": {
      "min": 62.0,
      "max": 250.0,
      "unit": "Da",
      "preferred_range": {
        "min": 76.0,
        "max": 150.0,
        "unit": "Da"
      }
    },
    "extended_range": {
      "min": 32.0,
      "max": 500.0,
      "unit": "Da"
    },
    "scientific_rationale": [
      "Small molecules (MW < 62 Da) generally have insufficient hydrogen bonding capacity",
      "Medium-sized molecules (76-150 Da) like glycerol have optimal balance of permeability and cryoprotection",
      "Larger molecules (> 250 Da) may have limited cell permeability but can still function as extracellular cryoprotectants"
    ],
    "example_compounds": {
      "optimal": ["glycerol (92.09 Da)", "propylene glycol (76.09 Da)", "DMSO (78.13 Da)"],
      "boundary": ["ethylene glycol (62.07 Da)", "trehalose (342.3 Da)"]
    },
    "implementation_sql": "molecular_weight BETWEEN 62.0 AND 250.0"
  }
}
```

#### 2. Hydrogen Bonding Capacity

```json
{
  "hydrogen_bonding_filter": {
    "description": "Hydrogen bonding parameters for water interaction and ice disruption",
    "hydrogen_bond_donors": {
      "min": 1,
      "max": 6,
      "preferred_range": {
        "min": 2,
        "max": 4
      },
      "scientific_rationale": [
        "Minimum of 1 donor required for interaction with water molecules",
        "Multiple donors (2-4) provide optimal ice crystal disruption",
        "Excessive donors (>6) may result in high viscosity and poor permeability"
      ]
    },
    "hydrogen_bond_acceptors": {
      "min": 2,
      "max": 10,
      "preferred_range": {
        "min": 2,
        "max": 6
      },
      "scientific_rationale": [
        "Minimum of 2 acceptors needed for effective water structuring",
        "Optimal range (2-6) balances water interaction with other properties",
        "High acceptor count (>10) often correlates with poor membrane penetration"
      ]
    },
    "hydrogen_bonding_ratio": {
      "description": "Ratio of acceptors to donors",
      "preferred_ratio": {
        "min": 0.75,
        "max": 3.0
      },
      "scientific_rationale": "Balanced ratio ensures proper interaction with water molecules"
    },
    "implementation_rdkit": "Lipinski.NumHDonors(mol) BETWEEN 1 AND 6 AND Lipinski.NumHAcceptors(mol) BETWEEN 2 AND 10"
  }
}
```

#### 3. Lipophilicity and Membrane Interaction

```json
{
  "lipophilicity_filter": {
    "description": "Properties affecting membrane penetration and water solubility",
    "logP": {
      "min": -3.0,
      "max": 1.5,
      "preferred_range": {
        "min": -2.0,
        "max": 0.5
      },
      "unit": "dimensionless",
      "scientific_rationale": [
        "Negative logP ensures water solubility",
        "Slightly negative to slightly positive logP (-2.0 to 0.5) provides optimal membrane permeability",
        "Values above 1.5 generally indicate poor water solubility unsuitable for cryoprotection"
      ]
    },
    "tPSA": {
      "min": 20.0,
      "max": 140.0,
      "preferred_range": {
        "min": 40.0,
        "max": 90.0
      },
      "unit": "Å²",
      "scientific_rationale": [
        "Topological Polar Surface Area correlates with membrane permeability",
        "Minimum of 20 Å² needed for water interaction",
        "Range of 40-90 Å² provides optimal balance between water binding and membrane penetration"
      ]
    },
    "example_compounds": {
      "optimal": ["DMSO (logP: -1.35, tPSA: 17.8)", "glycerol (logP: -1.76, tPSA: 60.7)"],
      "boundary": ["ethanol (logP: -0.31, tPSA: 20.2)", "trehalose (logP: -5.0, tPSA: 189.5)"]
    },
    "implementation_sql": "logP BETWEEN -3.0 AND 1.5 AND polar_surface_area BETWEEN 20.0 AND 140.0"
  }
}
```

### Cryoprotectant-Specific Properties

#### 1. Glass Transition Temperature (Tg)

```json
{
  "glass_transition_filter": {
    "description": "Filter based on glass transition temperature (Tg) of aqueous solutions",
    "known_compounds": {
      "min": -135.0,
      "max": -32.0,
      "unit": "°C",
      "optimal_range": {
        "min": -60.0,
        "max": -40.0,
        "unit": "°C"
      }
    },
    "predicted_tg": {
      "min": -100.0,
      "unit": "°C",
      "priority": "high"
    },
    "scientific_rationale": [
      "Higher Tg values indicate better glass-forming ability",
      "Tg values above -60°C are preferred for practical vitrification",
      "Compounds with Tg < -100°C typically form unstable glasses"
    ],
    "implementation_sql": "EXISTS (SELECT 1 FROM molecular_properties WHERE molecule_id = molecules.id AND property_type_id = (SELECT id FROM property_types WHERE name = 'glass_transition_temperature') AND numeric_value > -100.0)"
  }
}
```

#### 2. Cytotoxicity and Safety Parameters

```json
{
  "toxicity_filter": {
    "description": "Filter based on cytotoxicity and safety parameters",
    "ld50_oral": {
      "min": 1000.0,
      "unit": "mg/kg",
      "preferred_min": 5000.0,
      "priority": "critical"
    },
    "cell_viability": {
      "min": 70.0,
      "unit": "percent at 100mM",
      "preferred_min": 85.0
    },
    "excluded_toxicophores": [
      {
        "name": "epoxide",
        "smarts": "C1OC1",
        "rationale": "Reactive with proteins and DNA"
      },
      {
        "name": "acyl_halide",
        "smarts": "C(=O)[F,Cl,Br,I]",
        "rationale": "Highly reactive with biomolecules"
      },
      {
        "name": "aldehyde",
        "smarts": "[CH;R0]=O",
        "rationale": "Can crosslink proteins"
      }
    ],
    "scientific_rationale": [
      "Low toxicity critical for use with biological samples",
      "Cell viability above 70% at working concentrations required",
      "Exclusion of reactive functional groups that may cause damage"
    ],
    "implementation_rdkit": "not mol.HasSubstructMatch(Chem.MolFromSmarts('C1OC1')) AND not mol.HasSubstructMatch(Chem.MolFromSmarts('C(=O)[F,Cl,Br,I]'))"
  }
}
```

#### 3. Ice Crystal Inhibition Properties

```json
{
  "ice_inhibition_filter": {
    "description": "Filter based on ice crystal growth inhibition properties",
    "ice_recrystallization_inhibition": {
      "min": 0.3,
      "unit": "relative to PVA",
      "preferred_min": 0.5,
      "description": "Relative ice recrystallization inhibition activity compared to polyvinyl alcohol"
    },
    "structural_features": [
      {
        "name": "hydroxyl_spacing",
        "description": "Optimally spaced hydroxyl groups for ice binding",
        "smarts": "[OX2H][CX4][CX4][OX2H]"
      },
      {
        "name": "amide_group",
        "description": "Amide functionality that interacts with water structure",
        "smarts": "[NX3][CX3]=[OX1]"
      }
    ],
    "scientific_rationale": [
      "Inhibition of ice crystal growth reduces mechanical damage during freezing/thawing",
      "Specific structural features correlate with ice crystal inhibition",
      "Compounds with multiple ice-binding motifs have synergistic effects"
    ],
    "implementation_hybrid": "CASE WHEN EXISTS (SELECT 1 FROM molecular_properties WHERE molecule_id = molecules.id AND property_type_id = (SELECT id FROM property_types WHERE name = 'ice_recrystallization_inhibition') AND numeric_value >= 0.3) THEN TRUE WHEN molecules.smiles IS NOT NULL AND (molecules_contain_substructure(molecules.id, '[OX2H][CX4][CX4][OX2H]') OR molecules_contain_substructure(molecules.id, '[NX3][CX3]=[OX1]')) THEN TRUE ELSE FALSE END"
  }
}
```

#### 4. Cell Membrane Permeability

```json
{
  "permeability_filter": {
    "description": "Filter based on cell membrane permeability parameters",
    "experimental_permeability": {
      "min": 1.0e-6,
      "unit": "cm/s",
      "preferred_min": 5.0e-6,
      "description": "Measured permeability coefficient"
    },
    "calculated_parameters": {
      "rule_of_five_violations": {
        "max": 1,
        "description": "Lipinski's Rule of Five violations"
      },
      "rule_of_three_compliant": {
        "value": true,
        "description": "Compliance with Rule of Three for fragment permeability"
      }
    },
    "scientific_rationale": [
      "Sufficient permeability allows intracellular cryoprotection",
      "Multiple violations of permeability rules indicate poor cell penetration",
      "Balance required between size, lipophilicity, and hydrogen bonding"
    ],
    "implementation_hybrid": "(CASE WHEN EXISTS (SELECT 1 FROM molecular_properties WHERE molecule_id = molecules.id AND property_type_id = (SELECT id FROM property_types WHERE name = 'membrane_permeability') AND numeric_value >= 1.0e-6) THEN TRUE ELSE (lipinski_violations(molecules.id) <= 1) END)"
  }
}
```

### Functional Group Filters

```json
{
  "functional_group_filters": {
    "description": "Filters based on specific functional groups relevant to cryoprotection",
    "required_groups": [
      {
        "name": "hydroxyl",
        "smarts": "[OX2H]",
        "min_count": 1,
        "preferred_count": 3,
        "max_count": 8,
        "description": "Hydroxyl groups essential for water interaction",
        "scientific_rationale": "Forms hydrogen bonds with water, disrupting ice crystal formation"
      }
    ],
    "favorable_groups": [
      {
        "name": "ether",
        "smarts": "[OX2]([CX4])[CX4]",
        "score_bonus": 0.5,
        "description": "Ether linkages provide flexibility with hydrogen bond accepting capacity",
        "scientific_rationale": "Common in many effective cryoprotectants like PEG derivatives"
      },
      {
        "name": "amide",
        "smarts": "[NX3][CX3]=[OX1]",
        "score_bonus": 0.7,
        "description": "Amide groups interact strongly with water",
        "scientific_rationale": "Found in effective cryoprotectants like formamide derivatives"
      },
      {
        "name": "sulfoxide",
        "smarts": "[#16X3](=[OX1])([#6])[#6]",
        "score_bonus": 1.0,
        "description": "Sulfoxide group as in DMSO",
        "scientific_rationale": "Powerful hydrogen bond acceptor that disrupts water structure"
      }
    ],
    "unfavorable_groups": [
      {
        "name": "carboxylic_acid",
        "smarts": "[CX3](=O)[OX2H]",
        "score_penalty": 0.3,
        "description": "Carboxylic acids may increase toxicity at higher concentrations",
        "rationale": "pH effects and potential reactivity with biomolecules",
        "exceptions": ["Small alpha-hydroxy acids may be acceptable"]
      }
    ],
    "implementation_sql": "EXISTS (SELECT 1 FROM molecular_substructures WHERE molecule_id = molecules.id AND substructure_type = 'hydroxyl') AND NOT EXISTS (SELECT 1 FROM molecular_substructures WHERE molecule_id = molecules.id AND substructure_type IN ('epoxide', 'acyl_halide', 'isocyanate'))"
  }
}
```

## Implementation in CryoProtect

### Database Query Implementation

The property-based filtering criteria are implemented as configurable SQL queries in the database layer:

```sql
-- Example: Function to find cryoprotectant candidates based on physicochemical properties
CREATE OR REPLACE FUNCTION find_cryoprotectant_candidates(
    min_mw FLOAT DEFAULT 32.0,
    max_mw FLOAT DEFAULT 500.0,
    min_logp FLOAT DEFAULT -3.0,
    max_logp FLOAT DEFAULT 1.5,
    min_hbd INTEGER DEFAULT 1,
    max_hbd INTEGER DEFAULT 6,
    min_hba INTEGER DEFAULT 2,
    max_hba INTEGER DEFAULT 10,
    require_hydroxyl BOOLEAN DEFAULT TRUE
)
RETURNS TABLE (
    id UUID,
    name TEXT,
    smiles TEXT,
    molecular_weight FLOAT,
    logp FLOAT,
    hydrogen_bond_donors INTEGER,
    hydrogen_bond_acceptors INTEGER,
    cryoprotectant_score FLOAT
) AS $$
BEGIN
    RETURN QUERY
    WITH property_filters AS (
        SELECT m.id, m.name, m.smiles
        FROM molecules m
        WHERE 
            m.molecular_weight BETWEEN min_mw AND max_mw
            AND (
                -- Check logP from properties table
                EXISTS (
                    SELECT 1
                    FROM molecular_properties mp
                    JOIN property_types pt ON mp.property_type_id = pt.id
                    WHERE mp.molecule_id = m.id
                      AND pt.name = 'logP'
                      AND mp.numeric_value BETWEEN min_logp AND max_logp
                )
                OR
                -- Fallback to cached properties in JSONB
                (m.properties ? 'logP' AND 
                 (m.properties->>'logP')::float BETWEEN min_logp AND max_logp)
            )
            AND (
                -- Check hydrogen bond donors
                EXISTS (
                    SELECT 1
                    FROM molecular_properties mp
                    JOIN property_types pt ON mp.property_type_id = pt.id
                    WHERE mp.molecule_id = m.id
                      AND pt.name = 'hydrogen_bond_donors'
                      AND mp.numeric_value BETWEEN min_hbd AND max_hbd
                )
                OR
                -- Fallback to cached properties
                (m.properties ? 'hydrogen_bond_donors' AND 
                 (m.properties->>'hydrogen_bond_donors')::int BETWEEN min_hbd AND max_hbd)
            )
            AND (
                -- Check hydrogen bond acceptors
                EXISTS (
                    SELECT 1
                    FROM molecular_properties mp
                    JOIN property_types pt ON mp.property_type_id = pt.id
                    WHERE mp.molecule_id = m.id
                      AND pt.name = 'hydrogen_bond_acceptors'
                      AND mp.numeric_value BETWEEN min_hba AND max_hba
                )
                OR
                -- Fallback to cached properties
                (m.properties ? 'hydrogen_bond_acceptors' AND 
                 (m.properties->>'hydrogen_bond_acceptors')::int BETWEEN min_hba AND max_hba)
            )
            AND (
                -- Hydroxyl group requirement
                NOT require_hydroxyl 
                OR
                EXISTS (
                    SELECT 1
                    FROM molecular_substructures ms
                    WHERE ms.molecule_id = m.id
                      AND ms.substructure_type = 'hydroxyl'
                )
            )
    ),
    -- Get property values for scoring and output
    property_values AS (
        SELECT 
            pf.id,
            pf.name,
            pf.smiles,
            m.molecular_weight,
            COALESCE(
                (SELECT mp.numeric_value 
                 FROM molecular_properties mp 
                 JOIN property_types pt ON mp.property_type_id = pt.id 
                 WHERE mp.molecule_id = pf.id AND pt.name = 'logP'),
                (pf.properties->>'logP')::float
            ) AS logp,
            COALESCE(
                (SELECT mp.numeric_value 
                 FROM molecular_properties mp 
                 JOIN property_types pt ON mp.property_type_id = pt.id 
                 WHERE mp.molecule_id = pf.id AND pt.name = 'hydrogen_bond_donors'),
                (pf.properties->>'hydrogen_bond_donors')::int
            ) AS hydrogen_bond_donors,
            COALESCE(
                (SELECT mp.numeric_value 
                 FROM molecular_properties mp 
                 JOIN property_types pt ON mp.property_type_id = pt.id 
                 WHERE mp.molecule_id = pf.id AND pt.name = 'hydrogen_bond_acceptors'),
                (pf.properties->>'hydrogen_bond_acceptors')::int
            ) AS hydrogen_bond_acceptors,
            -- Calculate cryoprotectant score
            -- Normalized score from 0-100 based on optimal property ranges
            (
                -- MW component: optimal around 100-150 Da
                (CASE 
                    WHEN m.molecular_weight BETWEEN 100 AND 150 THEN 1.0
                    WHEN m.molecular_weight < 100 THEN 1.0 - (100 - m.molecular_weight) / 100
                    ELSE 1.0 - (m.molecular_weight - 150) / 350
                END) * 20 +
                
                -- LogP component: optimal around -1.5 to 0
                (CASE
                    WHEN logp BETWEEN -1.5 AND 0 THEN 1.0
                    WHEN logp < -1.5 THEN 1.0 - (ABS(logp + 1.5) / 1.5)
                    ELSE 1.0 - (logp / 1.5)
                END) * 25 +
                
                -- Hydrogen bonding component
                (CASE
                    WHEN hydrogen_bond_donors BETWEEN 2 AND 4 THEN 1.0
                    ELSE 0.7
                END) * 15 +
                
                (CASE
                    WHEN hydrogen_bond_acceptors BETWEEN 2 AND 6 THEN 1.0
                    ELSE 0.7
                END) * 15 +
                
                -- Bonus for known effective groups
                (CASE
                    WHEN EXISTS (
                        SELECT 1 FROM molecular_substructures 
                        WHERE molecule_id = pf.id AND substructure_type = 'sulfoxide'
                    ) THEN 15
                    WHEN EXISTS (
                        SELECT 1 FROM molecular_substructures 
                        WHERE molecule_id = pf.id AND substructure_type = 'amide'
                    ) THEN 10
                    WHEN EXISTS (
                        SELECT 1 FROM molecular_substructures 
                        WHERE molecule_id = pf.id AND substructure_type = 'ether'
                    ) THEN 5
                    ELSE 0
                END) +
                
                -- Bonus for favorable Tg if known
                (CASE
                    WHEN EXISTS (
                        SELECT 1 FROM molecular_properties mp
                        JOIN property_types pt ON mp.property_type_id = pt.id
                        WHERE mp.molecule_id = pf.id 
                          AND pt.name = 'glass_transition_temperature'
                          AND mp.numeric_value > -60
                    ) THEN 20
                    ELSE 0
                END)
            ) AS cryoprotectant_score
        FROM property_filters pf
        JOIN molecules m ON pf.id = m.id
    )
    
    -- Return final filtered and scored results
    SELECT
        id,
        name,
        smiles,
        molecular_weight,
        logp,
        hydrogen_bond_donors,
        hydrogen_bond_acceptors,
        cryoprotectant_score
    FROM
        property_values
    ORDER BY
        cryoprotectant_score DESC;
END;
$$ LANGUAGE plpgsql;
```

### Application Layer Implementation

In the application layer, these filtering criteria are implemented through a configurable filtering system:

```python
class CryoprotectantFilterManager:
    """
    Manages and applies property-based filters for cryoprotectant molecules.
    
    This class loads filter configurations from JSON files, applies
    filters to molecules in the database, and provides scoring mechanisms
    to rank potential cryoprotectants.
    """
    
    def __init__(self, db_connection, config_path=None):
        """
        Initialize filter manager.
        
        Parameters:
        -----------
        db_connection : psycopg2.connection
            Database connection
        config_path : str, optional
            Path to configuration file
        """
        self.db_connection = db_connection
        self.filters = {}
        
        # Load default filters
        self._load_default_filters()
        
        # Load custom filters from config file
        if config_path:
            self._load_filters_from_file(config_path)
    
    def _load_default_filters(self):
        """Load default cryoprotectant filter profiles"""
        self.filters["standard"] = {
            "name": "Standard Cryoprotectant Filter",
            "description": "General-purpose filter for potential cryoprotectants",
            "properties": {
                "molecular_weight": {"min": 62.0, "max": 250.0, "weight": 0.15},
                "logP": {"min": -3.0, "max": 1.5, "weight": 0.20},
                "hydrogen_bond_donors": {"min": 1, "max": 6, "weight": 0.15},
                "hydrogen_bond_acceptors": {"min": 2, "max": 10, "weight": 0.15}
            },
            "required_substructures": [
                {"name": "hydroxyl", "smarts": "[OX2H]", "min_count": 1}
            ],
            "favorable_substructures": [
                {"name": "sulfoxide", "smarts": "[#16X3](=[OX1])([#6])[#6]", "score_bonus": 10},
                {"name": "amide", "smarts": "[NX3][CX3]=[OX1]", "score_bonus": 5},
                {"name": "ether", "smarts": "[OX2]([CX4])[CX4]", "score_bonus": 3}
            ],
            "excluded_substructures": [
                {"name": "epoxide", "smarts": "C1OC1"},
                {"name": "acyl_halide", "smarts": "C(=O)[F,Cl,Br,I]"},
                {"name": "isocyanate", "smarts": "[NX2]=[CX2]=[OX1]"}
            ]
        }
        
        self.filters["intracellular"] = {
            "name": "Intracellular Cryoprotectant Filter",
            "description": "Filter for cryoprotectants with good cell penetration",
            "properties": {
                "molecular_weight": {"min": 32.0, "max": 180.0, "weight": 0.20},
                "logP": {"min": -2.0, "max": 1.0, "weight": 0.25},
                "polar_surface_area": {"max": 90.0, "weight": 0.15},
                "hydrogen_bond_donors": {"min": 1, "max": 4, "weight": 0.10},
                "hydrogen_bond_acceptors": {"min": 1, "max": 7, "weight": 0.10}
            },
            # Additional configuration similar to standard filter
        }
        
        self.filters["vitrification"] = {
            "name": "Vitrification Cryoprotectant Filter",
            "description": "Filter optimized for glass-forming cryoprotectants",
            "properties": {
                "molecular_weight": {"min": 60.0, "max": 400.0, "weight": 0.10},
                "logP": {"min": -3.0, "max": 0.5, "weight": 0.15},
                "hydrogen_bond_donors": {"min": 2, "max": 8, "weight": 0.20},
                "hydrogen_bond_acceptors": {"min": 3, "max": 12, "weight": 0.20}
            },
            # Additional configuration similar to standard filter
        }
    
    def _load_filters_from_file(self, config_path):
        """Load filter configurations from JSON file"""
        import json
        
        try:
            with open(config_path, 'r') as f:
                custom_filters = json.load(f)
                
            # Validate and add custom filters
            for filter_id, filter_config in custom_filters.items():
                if self._validate_filter_config(filter_config):
                    self.filters[filter_id] = filter_config
        except Exception as e:
            print(f"Error loading filter configuration: {e}")
    
    def _validate_filter_config(self, filter_config):
        """Validate filter configuration structure"""
        # Basic validation of required fields
        required_fields = ["name", "properties"]
        if not all(field in filter_config for field in required_fields):
            return False
            
        # Validate property configurations
        for prop, config in filter_config["properties"].items():
            if not isinstance(config, dict):
                return False
                
        return True
    
    def apply_filter(self, filter_id="standard", limit=100, offset=0):
        """
        Apply a specific filter profile to the database.
        
        Parameters:
        -----------
        filter_id : str
            ID of the filter profile to apply
        limit : int
            Maximum number of results to return
        offset : int
            Offset for pagination
            
        Returns:
        --------
        list
            List of filtered molecules with their properties and scores
        """
        if filter_id not in self.filters:
            raise ValueError(f"Unknown filter profile: {filter_id}")
            
        filter_config = self.filters[filter_id]
        
        # Build SQL query from filter configuration
        query, params = self._build_filter_query(filter_config, limit, offset)
        
        # Execute query
        with self.db_connection.cursor() as cursor:
            cursor.execute(query, params)
            results = cursor.fetchall()
            
            # Get column names
            columns = [desc[0] for desc in cursor.description]
            
            # Convert to dictionaries
            molecules = [dict(zip(columns, row)) for row in results]
            
        return molecules
    
    def _build_filter_query(self, filter_config, limit, offset):
        """Build SQL query from filter configuration"""
        # Extract property filters
        conditions = []
        params = {}
        
        # Process property filters
        for prop, config in filter_config["properties"].items():
            if "min" in config:
                conditions.append(f"({prop} >= %({prop}_min)s OR {prop} IS NULL)")
                params[f"{prop}_min"] = config["min"]
            if "max" in config:
                conditions.append(f"({prop} <= %({prop}_max)s OR {prop} IS NULL)")
                params[f"{prop}_max"] = config["max"]
        
        # Process required substructures
        if "required_substructures" in filter_config:
            for idx, substructure in enumerate(filter_config["required_substructures"]):
                conditions.append(f"EXISTS (SELECT 1 FROM molecular_substructures WHERE molecule_id = molecules.id AND substructure_type = %(req_substructure_{idx})s)")
                params[f"req_substructure_{idx}"] = substructure["name"]
        
        # Process excluded substructures
        if "excluded_substructures" in filter_config:
            for idx, substructure in enumerate(filter_config["excluded_substructures"]):
                conditions.append(f"NOT EXISTS (SELECT 1 FROM molecular_substructures WHERE molecule_id = molecules.id AND substructure_type = %(excl_substructure_{idx})s)")
                params[f"excl_substructure_{idx}"] = substructure["name"]
        
        # Build WHERE clause
        where_clause = " AND ".join(conditions) if conditions else "TRUE"
        
        # Build scoring expression
        scoring_terms = []
        
        # Property-based scoring
        for prop, config in filter_config["properties"].items():
            if "weight" in config:
                # Normalize property to 0-1 range and apply weight
                if "min" in config and "max" in config:
                    scoring_terms.append(f"""
                        CASE
                            WHEN {prop} IS NULL THEN 0
                            WHEN {prop} BETWEEN %({prop}_score_min)s AND %({prop}_score_max)s THEN %({prop}_weight)s
                            WHEN {prop} < %({prop}_score_min)s THEN 
                                %({prop}_weight)s * (1.0 - (%({prop}_score_min)s - {prop}) / %({prop}_score_min)s)
                            ELSE 
                                %({prop}_weight)s * (1.0 - ({prop} - %({prop}_score_max)s) / %({prop}_score_max)s)
                        END
                    """)
                    params[f"{prop}_score_min"] = config["min"]
                    params[f"{prop}_score_max"] = config["max"]
                    params[f"{prop}_weight"] = config["weight"]
                elif "min" in config:
                    scoring_terms.append(f"""
                        CASE
                            WHEN {prop} IS NULL THEN 0
                            WHEN {prop} >= %({prop}_score_min)s THEN %({prop}_weight)s
                            ELSE %({prop}_weight)s * ({prop} / %({prop}_score_min)s)
                        END
                    """)
                    params[f"{prop}_score_min"] = config["min"]
                    params[f"{prop}_weight"] = config["weight"]
                elif "max" in config:
                    scoring_terms.append(f"""
                        CASE
                            WHEN {prop} IS NULL THEN 0
                            WHEN {prop} <= %({prop}_score_max)s THEN %({prop}_weight)s
                            ELSE %({prop}_weight)s * (1.0 - ({prop} - %({prop}_score_max)s) / %({prop}_score_max)s)
                        END
                    """)
                    params[f"{prop}_score_max"] = config["max"]
                    params[f"{prop}_weight"] = config["weight"]
        
        # Substructure-based scoring
        if "favorable_substructures" in filter_config:
            for idx, substructure in enumerate(filter_config["favorable_substructures"]):
                scoring_terms.append(f"""
                    CASE
                        WHEN EXISTS (
                            SELECT 1 FROM molecular_substructures 
                            WHERE molecule_id = molecules.id 
                            AND substructure_type = %(fav_substructure_{idx})s
                        ) THEN %(fav_substructure_bonus_{idx})s
                        ELSE 0
                    END
                """)
                params[f"fav_substructure_{idx}"] = substructure["name"]
                params[f"fav_substructure_bonus_{idx}"] = substructure["score_bonus"]
        
        # Complete scoring expression
        score_expr = " + ".join(scoring_terms) if scoring_terms else "0"
        
        # Build complete query
        query = f"""
            SELECT 
                id,
                name,
                smiles,
                molecular_weight,
                logP,
                hydrogen_bond_donors,
                hydrogen_bond_acceptors,
                polar_surface_area,
                ({score_expr}) AS cryoprotectant_score
            FROM 
                molecules
            WHERE 
                {where_clause}
            ORDER BY 
                cryoprotectant_score DESC
            LIMIT %(limit)s
            OFFSET %(offset)s
        """
        
        params["limit"] = limit
        params["offset"] = offset
        
        return query, params
    
    def get_available_filters(self):
        """
        Get list of available filter profiles.
        
        Returns:
        --------
        list
            List of filter profile metadata
        """
        return [
            {
                "id": filter_id,
                "name": config["name"],
                "description": config.get("description", "")
            }
            for filter_id, config in self.filters.items()
        ]
    
    def create_custom_filter(self, filter_id, filter_config):
        """
        Create a new custom filter profile.
        
        Parameters:
        -----------
        filter_id : str
            ID for the new filter profile
        filter_config : dict
            Filter configuration
            
        Returns:
        --------
        bool
            True if successful, False otherwise
        """
        # Validate configuration
        if not self._validate_filter_config(filter_config):
            return False
            
        # Add to filters dictionary
        self.filters[filter_id] = filter_config
        return True
```

### API Implementation

The filtering system is exposed through a RESTful API:

```python
@app.route('/api/v1/molecules/filter', methods=['POST'])
def filter_molecules():
    """
    Filter molecules based on property and structural criteria.
    
    Request Body:
    {
        "filter_id": "standard",  # Use a predefined filter profile
        "custom_filter": {        # Or provide a custom filter
            "properties": {
                "molecular_weight": {"min": 60.0, "max": 200.0},
                "logP": {"max": 1.0}
            },
            "required_substructures": ["hydroxyl"],
            "excluded_substructures": ["epoxide"]
        },
        "pagination": {
            "page": 1,
            "page_size": 20
        },
        "sort": {
            "field": "cryoprotectant_score",
            "order": "desc"
        }
    }
    
    Response:
    {
        "molecules": [...],
        "pagination": {
            "page": 1,
            "page_size": 20,
            "total": 156,
            "pages": 8
        }
    }
    """
    # Get request data
    data = request.get_json()
    
    # Get pagination parameters
    page = data.get('pagination', {}).get('page', 1)
    page_size = data.get('pagination', {}).get('page_size', 20)
    offset = (page - 1) * page_size
    
    # Initialize filter manager
    filter_manager = CryoprotectantFilterManager(db_connection)
    
    # Apply filter
    if 'custom_filter' in data:
        # Create temporary custom filter
        custom_filter_id = f"temp_{uuid.uuid4()}"
        filter_manager.create_custom_filter(custom_filter_id, data['custom_filter'])
        molecules = filter_manager.apply_filter(custom_filter_id, limit=page_size, offset=offset)
    else:
        # Use predefined filter
        filter_id = data.get('filter_id', 'standard')
        molecules = filter_manager.apply_filter(filter_id, limit=page_size, offset=offset)
    
    # Get total count (for pagination)
    total_count = len(molecules)  # This should be optimized for actual counts
    
    # Format response
    response = {
        'molecules': molecules,
        'pagination': {
            'page': page,
            'page_size': page_size,
            'total': total_count,
            'pages': math.ceil(total_count / page_size)
        }
    }
    
    return jsonify(response)
```

## Extending the Filtering System

### User-Configurable Filters

Users can create and save custom filter configurations through the UI:

```json
{
  "user_saved_filters": [
    {
      "id": "user_intracellular_cryoprotectants",
      "name": "My Intracellular Cryoprotectants",
      "description": "Custom filter for small, cell-permeable cryoprotectants",
      "created_by": "user123",
      "created_at": "2025-05-01T10:32:16Z",
      "is_public": false,
      "properties": {
        "molecular_weight": {"min": 40.0, "max": 140.0, "weight": 0.25},
        "logP": {"min": -2.0, "max": 0.5, "weight": 0.25}
      },
      "required_substructures": [
        {"name": "hydroxyl", "smarts": "[OX2H]", "min_count": 2}
      ],
      "favorable_substructures": [
        {"name": "sulfoxide", "smarts": "[#16X3](=[OX1])([#6])[#6]", "score_bonus": 15}
      ]
    }
  ]
}
```

### Machine Learning Integration

The filtering system integrates with machine learning models for property prediction and cryoprotectant efficacy prediction:

```python
def apply_ml_filter(molecules, ml_model_path):
    """
    Apply machine learning model to further filter and score molecules.
    
    Parameters:
    -----------
    molecules : list
        List of molecule dictionaries from initial filtering
    ml_model_path : str
        Path to trained ML model
        
    Returns:
    --------
    list
        Filtered and rescored molecules
    """
    # Load ML model
    import joblib
    model = joblib.load(ml_model_path)
    
    # Convert SMILES to RDKit molecules and calculate features
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    
    features = []
    valid_indices = []
    
    for i, mol_dict in enumerate(molecules):
        mol = Chem.MolFromSmiles(mol_dict['smiles'])
        if mol:
            # Calculate descriptors for ML model
            descriptors = []
            for desc_name in model.feature_names:
                if hasattr(Descriptors, desc_name):
                    descriptors.append(getattr(Descriptors, desc_name)(mol))
                else:
                    descriptors.append(0)  # Default for unknown descriptors
            
            features.append(descriptors)
            valid_indices.append(i)
    
    # Apply model to predict efficacy scores
    efficacy_scores = model.predict_proba(features)[:, 1]  # Probability of positive class
    
    # Update molecules with ML scores
    for i, idx in enumerate(valid_indices):
        molecules[idx]['ml_efficacy_score'] = float(efficacy_scores[i])
        
        # Combine with rule-based score for final ranking
        molecules[idx]['combined_score'] = (
            0.6 * molecules[idx]['cryoprotectant_score'] + 
            0.4 * molecules[idx]['ml_efficacy_score'] * 100  # Scale to 0-100
        )
    
    # Remove molecules without ML scores
    filtered_molecules = [mol for i, mol in enumerate(molecules) if i in valid_indices]
    
    # Sort by combined score
    filtered_molecules.sort(key=lambda x: x['combined_score'], reverse=True)
    
    return filtered_molecules
```