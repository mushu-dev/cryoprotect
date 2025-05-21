# CryoProtect Data Pipeline Specifications

## Data Pipeline Overview

The CryoProtect data pipeline manages the flow of molecular data from external sources through processing, filtering, and into the database. This document specifies the input and output formats at each stage of the pipeline.

## 1. External Data Acquisition

### PubChem Data Import

**Input**: PubChem API requests for cryoprotectant-related compounds

```json
{
  "import_config": {
    "source": "pubchem",
    "method": "api",
    "query_type": "compound_search",
    "parameters": {
      "query": "cryoprotectant OR \"cryo protectant\" OR \"antifreeze compound\"",
      "max_compounds": 5000,
      "require_structures": true,
      "properties": [
        "MolecularFormula",
        "MolecularWeight",
        "CanonicalSMILES",
        "InChIKey",
        "IUPACName",
        "XLogP",
        "HydrogenBondDonorCount",
        "HydrogenBondAcceptorCount",
        "RotatableBondCount",
        "TPSA",
        "HeavyAtomCount",
        "ExactMass",
        "Complexity",
        "IsomericSMILES"
      ]
    },
    "rate_limit": {
      "requests_per_second": 5,
      "max_concurrent": 3
    }
  }
}
```

**Output**: Standardized molecule data objects

```json
{
  "molecules": [
    {
      "external_id": "PubChem:753",
      "name": "Glycerol",
      "synonyms": ["Glycerin", "1,2,3-Propanetriol", "1,2,3-Trihydroxypropane"],
      "smiles": "C(C(CO)O)O",
      "inchi": "InChI=1S/C3H8O3/c4-1-3(6)2-5/h3-6H,1-2H2",
      "inchikey": "PEDCQBHIVMGVHV-UHFFFAOYSA-N",
      "molecular_formula": "C3H8O3",
      "molecular_weight": 92.09,
      "properties": {
        "xlogp": -1.76,
        "hydrogen_bond_donors": 3,
        "hydrogen_bond_acceptors": 3,
        "rotatable_bonds": 2,
        "tpsa": 60.7,
        "heavy_atom_count": 6,
        "exact_mass": 92.0473,
        "complexity": 44.3
      },
      "source_data": {
        "cid": 753,
        "last_modified": "2024-01-15T09:23:47Z",
        "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/753"
      }
    }
  ],
  "metadata": {
    "query": "cryoprotectant",
    "timestamp": "2025-05-05T14:30:22Z",
    "total_found": 1243,
    "retrieved": 1000,
    "import_status": "partial"
  }
}
```

### ChEMBL Data Import

**Input**: ChEMBL API requests for bioactivity data

```json
{
  "import_config": {
    "source": "chembl",
    "method": "api",
    "parameters": {
      "target_keywords": ["cryoprotection", "freeze protection", "frost protection"],
      "activity_types": ["IC50", "EC50", "Potency", "Activity"],
      "min_confidence_score": 8,
      "max_compounds": 2000,
      "include_properties": true
    },
    "resilient_import": {
      "checkpoint_interval": 100,
      "retry_attempts": 3,
      "backoff_factor": 1.5
    }
  }
}
```

**Output**: Standardized bioactivity data linked to molecules

```json
{
  "bioactivities": [
    {
      "molecule": {
        "chembl_id": "CHEMBL2365062",
        "name": "Trehalose",
        "smiles": "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CO)O)O)O)O)O)O)O",
        "inchikey": "DBOWFGJUJNJVFF-SHIGYJNBSA-N",
        "molecular_weight": 342.3,
        "molecular_formula": "C12H22O11"
      },
      "activity": {
        "activity_id": "CHEMBL_ACT_25632654",
        "target_name": "Cryoprotection assay",
        "target_organism": "Saccharomyces cerevisiae",
        "activity_type": "EC50",
        "value": 2.3,
        "units": "mM",
        "relation": "=",
        "confidence_score": 9,
        "assay_description": "Protection against freeze-thaw damage in yeast cells",
        "reference": "J Biol Chem 2019;294(15):6086-6097",
        "pubmed_id": "30846547"
      }
    }
  ],
  "metadata": {
    "query_parameters": {
      "target_keywords": ["cryoprotection"],
      "activity_types": ["EC50"]
    },
    "timestamp": "2025-05-05T16:42:11Z",
    "total_found": 87,
    "retrieved": 87,
    "import_status": "complete"
  }
}
```

## 2. Data Normalization and Standardization

### Molecular Structure Standardization

**Input**: Raw molecular structures from various sources

```json
{
  "structures": [
    {
      "id": "compound1",
      "smiles": "C(C(CO)O)O", 
      "source": "PubChem"
    },
    {
      "id": "compound2",
      "smiles": "C(=O)NCCN.Cl", 
      "source": "ChEMBL"
    },
    {
      "id": "compound3",
      "smiles": "C[N+](C)(C)CCO.[Cl-]", 
      "source": "Manual"
    }
  ],
  "standardization_options": {
    "remove_salts": true,
    "neutralize_charges": true,
    "normalize_tautomers": true,
    "generate_canonical_smiles": true,
    "standardize_stereo": true
  }
}
```

**Output**: Standardized molecular structures

```json
{
  "standardized_structures": [
    {
      "id": "compound1",
      "original_smiles": "C(C(CO)O)O",
      "standardized_smiles": "C(C(CO)O)O",
      "inchikey": "PEDCQBHIVMGVHV-UHFFFAOYSA-N",
      "standardization_operations": []
    },
    {
      "id": "compound2",
      "original_smiles": "C(=O)NCCN.Cl",
      "standardized_smiles": "C(=O)NCCN",
      "inchikey": "LRQKBLKVPFOOAZ-UHFFFAOYSA-N",
      "standardization_operations": ["salt_removal"]
    },
    {
      "id": "compound3",
      "original_smiles": "C[N+](C)(C)CCO.[Cl-]",
      "standardized_smiles": "CN(C)(C)CCO",
      "inchikey": "MKUXAQIIEYXACX-UHFFFAOYSA-N",
      "standardization_operations": ["salt_removal", "charge_neutralization"]
    }
  ],
  "metadata": {
    "timestamp": "2025-05-05T17:10:45Z",
    "rdkit_version": "2024.03.4",
    "standardization_stats": {
      "total_compounds": 3,
      "unchanged": 1,
      "salt_removal": 2,
      "charge_neutralization": 1,
      "tautomer_standardization": 0,
      "stereo_standardization": 0,
      "failed": 0
    }
  }
}
```

### Property Calculation

**Input**: Standardized molecular structures

```json
{
  "molecules": [
    {
      "id": "mol1",
      "smiles": "C(C(CO)O)O",
      "name": "Glycerol"
    }
  ],
  "properties_to_calculate": [
    "molecular_weight",
    "logP",
    "hydrogen_bond_donors",
    "hydrogen_bond_acceptors",
    "rotatable_bonds",
    "polar_surface_area",
    "cryoprotectant_score"
  ],
  "custom_properties": [
    {
      "name": "hydroxyl_density",
      "description": "Ratio of hydroxyl groups to heavy atoms",
      "substructure": "[OX2H]",
      "calculation": "count/mol.GetNumHeavyAtoms()"
    }
  ]
}
```

**Output**: Molecules with calculated properties

```json
{
  "molecules_with_properties": [
    {
      "id": "mol1",
      "smiles": "C(C(CO)O)O",
      "name": "Glycerol",
      "properties": {
        "molecular_weight": 92.09,
        "logP": -1.76,
        "hydrogen_bond_donors": 3,
        "hydrogen_bond_acceptors": 3,
        "rotatable_bonds": 2,
        "polar_surface_area": 60.7,
        "hydroxyl_density": 0.5,
        "cryoprotectant_score": 87.4
      },
      "substructures": [
        {
          "type": "hydroxyl",
          "smarts": "[OX2H]",
          "count": 3,
          "matches": [[3], [5], [6]]
        }
      ]
    }
  ],
  "metadata": {
    "timestamp": "2025-05-05T17:25:33Z",
    "property_calculation_stats": {
      "total_molecules": 1,
      "successful": 1,
      "failed": 0
    }
  }
}
```

## 3. Filtering and Prioritization

### Initial Filtering

**Input**: Molecules with calculated properties

```json
{
  "filter_config": {
    "filter_id": "cryoprotectant_initial",
    "description": "Initial filtering of potential cryoprotectants",
    "property_filters": [
      {
        "property": "molecular_weight",
        "min": 32.0,
        "max": 500.0,
        "priority": "required"
      },
      {
        "property": "logP",
        "min": -3.0,
        "max": 1.5,
        "priority": "required"
      },
      {
        "property": "hydrogen_bond_donors",
        "min": 1,
        "priority": "required"
      }
    ],
    "structural_filters": {
      "required": [
        {
          "name": "hydroxyl",
          "smarts": "[OX2H]",
          "min_count": 1
        }
      ],
      "excluded": [
        {
          "name": "reactive_groups",
          "smarts": "[C,N]=[C,N,O][F,Cl,Br,I]"
        }
      ]
    }
  },
  "molecules": [
    "..." 
  ]
}
```

**Output**: Filtered molecules with filter status

```json
{
  "filtered_molecules": {
    "passed": [
      {
        "id": "mol1",
        "name": "Glycerol",
        "smiles": "C(C(CO)O)O",
        "filter_results": {
          "molecular_weight": {
            "value": 92.09,
            "passed": true,
            "reason": "92.09 is between 32.0 and 500.0"
          },
          "logP": {
            "value": -1.76,
            "passed": true,
            "reason": "-1.76 is between -3.0 and 1.5"
          },
          "hydrogen_bond_donors": {
            "value": 3,
            "passed": true,
            "reason": "3 is >= 1"
          },
          "structural_required": {
            "hydroxyl": {
              "passed": true,
              "matches": 3,
              "min_required": 1
            }
          }
        }
      }
    ],
    "failed": [
      {
        "id": "mol2",
        "name": "Benzene",
        "smiles": "c1ccccc1",
        "filter_results": {
          "molecular_weight": {
            "value": 78.11,
            "passed": true,
            "reason": "78.11 is between 32.0 and 500.0"
          },
          "logP": {
            "value": 2.13,
            "passed": false,
            "reason": "2.13 is greater than 1.5"
          },
          "hydrogen_bond_donors": {
            "value": 0,
            "passed": false,
            "reason": "0 is < 1"
          },
          "structural_required": {
            "hydroxyl": {
              "passed": false,
              "matches": 0,
              "min_required": 1
            }
          }
        }
      }
    ]
  },
  "filter_stats": {
    "total_molecules": 2,
    "passed": 1,
    "failed": 1,
    "failure_reasons": {
      "logP": 1,
      "hydrogen_bond_donors": 1,
      "structural_required.hydroxyl": 1
    }
  }
}
```

### Similarity-Based Enrichment

**Input**: Filtered molecules and reference cryoprotectants

```json
{
  "reference_molecules": [
    {
      "id": "ref1",
      "name": "Glycerol",
      "smiles": "C(C(CO)O)O",
      "importance_weight": 1.0
    },
    {
      "id": "ref2",
      "name": "DMSO",
      "smiles": "CS(C)=O",
      "importance_weight": 1.0
    },
    {
      "id": "ref3",
      "name": "Propylene Glycol",
      "smiles": "CC(O)CO",
      "importance_weight": 0.8
    }
  ],
  "candidate_molecules": [
    {
      "id": "cand1",
      "name": "1,3-Propanediol",
      "smiles": "C(CO)CO"
    },
    {
      "id": "cand2",
      "name": "1,2-Butanediol",
      "smiles": "CCC(O)CO"
    },
    {
      "id": "cand3",
      "name": "Dimethylacetamide",
      "smiles": "CC(=O)N(C)C"
    }
  ],
  "similarity_config": {
    "methods": [
      {
        "method": "morgan",
        "parameters": {
          "radius": 2,
          "nBits": 2048
        },
        "weight": 0.5
      },
      {
        "method": "maccs",
        "weight": 0.3
      },
      {
        "method": "rdkit",
        "weight": 0.2
      }
    ],
    "threshold": 0.6,
    "max_results": 100
  }
}
```

**Output**: Molecules with similarity scores

```json
{
  "similarity_results": [
    {
      "molecule": {
        "id": "cand1",
        "name": "1,3-Propanediol",
        "smiles": "C(CO)CO"
      },
      "similarities": [
        {
          "reference": "ref1",
          "reference_name": "Glycerol",
          "method": "morgan",
          "score": 0.85
        },
        {
          "reference": "ref2",
          "reference_name": "DMSO",
          "method": "morgan",
          "score": 0.32
        },
        {
          "reference": "ref3",
          "reference_name": "Propylene Glycol",
          "method": "morgan",
          "score": 0.76
        }
      ],
      "similarity_scores": {
        "morgan": 0.72,
        "maccs": 0.81,
        "rdkit": 0.68
      },
      "weighted_score": 0.75,
      "max_similarity": {
        "reference": "ref1",
        "reference_name": "Glycerol",
        "score": 0.85,
        "method": "morgan"
      }
    },
    {
      "molecule": {
        "id": "cand2",
        "name": "1,2-Butanediol",
        "smiles": "CCC(O)CO"
      },
      "similarities": [
        {"reference": "ref1", "reference_name": "Glycerol", "method": "morgan", "score": 0.72},
        {"reference": "ref2", "reference_name": "DMSO", "method": "morgan", "score": 0.28},
        {"reference": "ref3", "reference_name": "Propylene Glycol", "method": "morgan", "score": 0.89}
      ],
      "similarity_scores": {
        "morgan": 0.68,
        "maccs": 0.77,
        "rdkit": 0.64
      },
      "weighted_score": 0.70,
      "max_similarity": {
        "reference": "ref3",
        "reference_name": "Propylene Glycol",
        "score": 0.89,
        "method": "morgan"
      }
    },
    {
      "molecule": {
        "id": "cand3",
        "name": "Dimethylacetamide",
        "smiles": "CC(=O)N(C)C"
      },
      "similarities": [
        {"reference": "ref1", "reference_name": "Glycerol", "method": "morgan", "score": 0.21},
        {"reference": "ref2", "reference_name": "DMSO", "method": "morgan", "score": 0.67},
        {"reference": "ref3", "reference_name": "Propylene Glycol", "method": "morgan", "score": 0.18}
      ],
      "similarity_scores": {
        "morgan": 0.40,
        "maccs": 0.52,
        "rdkit": 0.35
      },
      "weighted_score": 0.43,
      "max_similarity": {
        "reference": "ref2",
        "reference_name": "DMSO",
        "score": 0.67,
        "method": "morgan"
      }
    }
  ],
  "metadata": {
    "timestamp": "2025-05-05T18:20:15Z",
    "similarity_stats": {
      "total_candidates": 3,
      "above_threshold": 2,
      "below_threshold": 1,
      "threshold": 0.6
    }
  }
}
```

### Advanced Scoring and Ranking

**Input**: Molecules with property and similarity data

```json
{
  "scoring_config": {
    "model_type": "weighted_combined",
    "property_weights": {
      "molecular_weight": { 
        "type": "gaussian",
        "target": 125.0,
        "width": 50.0,
        "weight": 0.1
      },
      "logP": {
        "type": "gaussian",
        "target": -1.0,
        "width": 1.0,
        "weight": 0.15
      },
      "hydrogen_bond_donors": {
        "type": "optimal_range",
        "min": 2,
        "max": 4,
        "weight": 0.1
      },
      "hydrogen_bond_acceptors": {
        "type": "optimal_range",
        "min": 2,
        "max": 6,
        "weight": 0.1
      }
    },
    "similarity_weight": 0.25,
    "structure_bonuses": [
      {
        "name": "hydroxyl_groups",
        "smarts": "[OX2H]",
        "score_per_match": 5,
        "max_score": 15,
        "weight": 0.1
      },
      {
        "name": "sulfoxide_group",
        "smarts": "[#16X3](=[OX1])([#6])[#6]",
        "flat_bonus": 10,
        "weight": 0.1
      },
      {
        "name": "amide_group",
        "smarts": "[NX3][CX3]=[OX1]",
        "flat_bonus": 8,
        "weight": 0.1
      }
    ],
    "ml_model_weight": 0.0
  },
  "molecules": [
    "..."
  ]
}
```

**Output**: Molecules with comprehensive scores

```json
{
  "scored_molecules": [
    {
      "id": "mol1",
      "name": "Glycerol",
      "smiles": "C(C(CO)O)O",
      "component_scores": {
        "properties": {
          "molecular_weight": {
            "value": 92.09,
            "score": 0.82,
            "contribution": 0.082
          },
          "logP": {
            "value": -1.76,
            "score": 0.71,
            "contribution": 0.107
          },
          "hydrogen_bond_donors": {
            "value": 3,
            "score": 1.0,
            "contribution": 0.1
          },
          "hydrogen_bond_acceptors": {
            "value": 3,
            "score": 1.0,
            "contribution": 0.1
          }
        },
        "similarity": {
          "score": 0.95,
          "contribution": 0.238,
          "best_match": "Glycerol (self)"
        },
        "structure_bonuses": {
          "hydroxyl_groups": {
            "matches": 3,
            "score": 15,
            "contribution": 0.1
          },
          "sulfoxide_group": {
            "matches": 0,
            "score": 0,
            "contribution": 0
          },
          "amide_group": {
            "matches": 0,
            "score": 0,
            "contribution": 0
          }
        },
        "ml_model": {
          "score": 0.0,
          "contribution": 0.0
        }
      },
      "total_score": 0.727,
      "normalized_score": 89.4,
      "rank": 1
    },
    {
      "id": "mol2",
      "name": "DMSO",
      "smiles": "CS(C)=O",
      "component_scores": "...",
      "total_score": 0.712,
      "normalized_score": 87.4,
      "rank": 2
    }
  ],
  "metadata": {
    "timestamp": "2025-05-05T18:45:22Z",
    "scoring_stats": {
      "molecules_scored": 50,
      "score_range": {
        "min": 0.123,
        "max": 0.727,
        "mean": 0.485,
        "median": 0.521
      }
    }
  }
}
```

## 4. Database Integration

### Database Insertion

**Input**: Processed molecules ready for database insertion

```json
{
  "molecules": [
    {
      "id": "temp_mol_id_12345",
      "name": "Glycerol",
      "smiles": "C(C(CO)O)O",
      "inchi": "InChI=1S/C3H8O3/c4-1-3(6)2-5/h3-6H,1-2H2",
      "inchikey": "PEDCQBHIVMGVHV-UHFFFAOYSA-N",
      "formula": "C3H8O3",
      "molecular_weight": 92.09,
      "data_source": "PubChem",
      "external_ids": {
        "pubchem_cid": 753,
        "cas": "56-81-5"
      },
      "properties": {
        "logP": -1.76,
        "hydrogen_bond_donors": 3,
        "hydrogen_bond_acceptors": 3,
        "rotatable_bonds": 2,
        "polar_surface_area": 60.7
      },
      "is_public": true,
      "cryoprotectant_score": 89.4,
      "substructures": [
        {"type": "hydroxyl", "smarts": "[OX2H]", "count": 3}
      ]
    }
  ],
  "options": {
    "update_existing": true,
    "batch_size": 100,
    "property_storage": "normalized",
    "return_ids": true
  }
}
```

**Output**: Database insertion results

```json
{
  "insertion_results": {
    "molecules": [
      {
        "temp_id": "temp_mol_id_12345",
        "db_id": "550e8400-e29b-41d4-a716-446655440000",
        "name": "Glycerol",
        "inchikey": "PEDCQBHIVMGVHV-UHFFFAOYSA-N",
        "status": "inserted"
      }
    ],
    "property_insertions": {
      "total": 5,
      "successful": 5,
      "failed": 0
    },
    "substructure_insertions": {
      "total": 1,
      "successful": 1,
      "failed": 0
    }
  },
  "metadata": {
    "timestamp": "2025-05-05T19:10:43Z",
    "database": "cryoprotect_production",
    "batch_id": "import_20250505_191043",
    "insertion_stats": {
      "new_molecules": 1,
      "updated_molecules": 0,
      "skipped_molecules": 0,
      "total_properties": 5
    }
  }
}
```

### Database Query and Retrieval

**Input**: Query parameters for molecule retrieval

```json
{
  "query": {
    "type": "advanced",
    "filters": {
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
          "value": 0.0
        }
      ],
      "substructure_filters": {
        "must_have": ["[OX2H]", "[OX2H][CX4][CX4][OX2H]"],
        "must_not_have": ["c1ccccc1"]
      },
      "cryoprotectant_score": {
        "min": 70.0
      }
    },
    "sort": {
      "field": "cryoprotectant_score",
      "direction": "desc"
    },
    "pagination": {
      "page": 1,
      "page_size": 20
    }
  }
}
```

**Output**: Query results with molecules and metadata

```json
{
  "results": {
    "molecules": [
      {
        "id": "550e8400-e29b-41d4-a716-446655440000",
        "name": "Glycerol",
        "smiles": "C(C(CO)O)O",
        "inchikey": "PEDCQBHIVMGVHV-UHFFFAOYSA-N",
        "formula": "C3H8O3",
        "molecular_weight": 92.09,
        "properties": {
          "logP": -1.76,
          "hydrogen_bond_donors": 3,
          "hydrogen_bond_acceptors": 3,
          "rotatable_bonds": 2,
          "polar_surface_area": 60.7,
          "glass_transition_temperature": -93.0
        },
        "cryoprotectant_score": 89.4,
        "is_public": true,
        "created_at": "2025-04-10T16:32:12Z",
        "updated_at": "2025-05-05T19:10:43Z"
      },
      {
        "id": "6ba7b810-9dad-11d1-80b4-00c04fd430c8",
        "name": "Propylene glycol",
        "smiles": "CC(O)CO",
        "inchikey": "DNIAPMSPPWPWGF-UHFFFAOYSA-N",
        "formula": "C3H8O2",
        "molecular_weight": 76.09,
        "properties": {
          "logP": -0.92,
          "hydrogen_bond_donors": 2,
          "hydrogen_bond_acceptors": 2,
          "rotatable_bonds": 1,
          "polar_surface_area": 40.46,
          "glass_transition_temperature": -108.0
        },
        "cryoprotectant_score": 84.7,
        "is_public": true,
        "created_at": "2025-04-10T16:32:15Z",
        "updated_at": "2025-05-05T19:10:44Z"
      }
    ],
    "pagination": {
      "page": 1,
      "page_size": 20,
      "total_results": 12,
      "total_pages": 1
    }
  },
  "metadata": {
    "timestamp": "2025-05-05T19:30:12Z",
    "query_execution_time_ms": 78,
    "query_params": {
      "molecular_weight": "60.0-200.0",
      "logP": "<0.0",
      "substructures": "[OX2H]"
    }
  }
}
```

## 5. Analytics and Feedback Loop

### Usage Analytics

**Input**: N/A (system-generated)

**Output**: Analytics on molecule filtering usage

```json
{
  "filtering_analytics": {
    "time_period": {
      "start": "2025-05-01T00:00:00Z",
      "end": "2025-05-07T23:59:59Z"
    },
    "query_statistics": {
      "total_queries": 853,
      "unique_users": 42,
      "average_results_returned": 23.5,
      "most_common_filters": [
        {
          "property": "molecular_weight",
          "count": 782,
          "percentage": 91.7
        },
        {
          "property": "logP",
          "count": 705,
          "percentage": 82.6
        },
        {
          "substructure": "hydroxyl",
          "count": 624,
          "percentage": 73.2
        }
      ],
      "most_viewed_molecules": [
        {
          "id": "550e8400-e29b-41d4-a716-446655440000",
          "name": "Glycerol",
          "views": 142
        },
        {
          "id": "6ba7b810-9dad-11d1-80b4-00c04fd430c8",
          "name": "Propylene glycol",
          "views": 118
        }
      ]
    },
    "filter_effectiveness": {
      "average_cryoprotectant_score": 67.3,
      "property_distribution": {
        "molecular_weight": {
          "min": 32.04,
          "max": 487.52,
          "mean": 153.84,
          "median": 121.6,
          "p10": 62.07,
          "p90": 302.4
        }
      },
      "experimental_feedback": {
        "total_experimental_results": 28,
        "successful_experiments": 21,
        "failed_experiments": 7,
        "success_rate": 0.75,
        "score_correlation": 0.68
      }
    }
  },
  "metadata": {
    "report_generation_time": "2025-05-08T01:15:33Z",
    "data_sources": ["query_logs", "user_events", "experimental_results"]
  }
}
```

### Experimental Feedback

**Input**: Experimental results for predicted cryoprotectants

```json
{
  "experimental_results": [
    {
      "molecule_id": "550e8400-e29b-41d4-a716-446655440000",
      "name": "Glycerol",
      "experiment_type": "cell_viability",
      "experimental_conditions": {
        "cell_type": "CHO-K1",
        "concentration": 10.0,
        "concentration_unit": "%w/v",
        "freezing_protocol": "slow_cooling_1C_per_min",
        "temperature": -80.0,
        "storage_duration": 7,
        "storage_duration_unit": "days"
      },
      "results": {
        "viability_post_thaw": 87.2,
        "viability_unit": "percent",
        "control_viability": 12.4,
        "protocol_effectiveness": "high",
        "notes": "Standard cryoprotectant performed as expected"
      },
      "experiment_metadata": {
        "experimenter": "user123",
        "experiment_date": "2025-05-03T14:30:00Z",
        "experimental_id": "EXP-2025-0503-01"
      }
    },
    {
      "molecule_id": "6ba7b810-9dad-11d1-80b4-00c04fd430c8",
      "name": "Propylene glycol",
      "experiment_type": "cell_viability",
      "experimental_conditions": "...",
      "results": "...",
      "experiment_metadata": "..."
    }
  ]
}
```

**Output**: Updated model and scoring parameters

```json
{
  "model_updates": {
    "updated_weights": {
      "molecular_weight": 0.12,
      "logP": 0.18,
      "hydrogen_bond_donors": 0.15,
      "hydrogen_bond_acceptors": 0.12,
      "similarity": 0.23,
      "structure_bonuses": 0.20
    },
    "property_target_adjustments": {
      "molecular_weight": {
        "target": 115.0,
        "width": 45.0
      },
      "logP": {
        "target": -1.2,
        "width": 0.9
      }
    },
    "effectiveness_prediction": {
      "model_performance": {
        "train_accuracy": 0.82,
        "validation_accuracy": 0.78,
        "f1_score": 0.81,
        "auc_roc": 0.84
      },
      "feature_importance": [
        {"feature": "hydrogen_bond_donors", "importance": 0.28},
        {"feature": "molecular_weight", "importance": 0.21},
        {"feature": "logP", "importance": 0.19},
        {"feature": "similarity_to_glycerol", "importance": 0.17},
        {"feature": "hydroxyl_count", "importance": 0.15}
      ]
    }
  },
  "metadata": {
    "update_timestamp": "2025-05-08T02:30:45Z",
    "experimental_data_points": 42,
    "model_version": "2.3.1"
  }
}
```

## Debugging and Diagnostics

### Error Handling

**Input**: Database insertion errors

```json
{
  "error_context": {
    "operation": "database_insertion",
    "batch_id": "import_20250505_191043",
    "error_count": 3,
    "errors": [
      {
        "molecule_id": "temp_mol_id_78901",
        "error_type": "duplicate_inchikey",
        "error_message": "Molecule with InChIKey MYMOFIZGZYHOMD-UHFFFAOYSA-N already exists with ID 7a8b9c10-11d1-e2f3-g4h5-i6j7k8l9m0n1",
        "attempted_operation": "insert",
        "timestamp": "2025-05-05T19:09:21Z"
      },
      {
        "molecule_id": "temp_mol_id_78902",
        "error_type": "invalid_smiles",
        "error_message": "Invalid SMILES string: C1CC(C)C1)",
        "attempted_operation": "insert",
        "timestamp": "2025-05-05T19:09:22Z"
      },
      {
        "molecule_id": "temp_mol_id_78903",
        "error_type": "property_validation",
        "error_message": "Value '-12.5' for property 'logP' is outside valid range [-3.0, 10.0]",
        "attempted_operation": "insert_property",
        "timestamp": "2025-05-05T19:09:23Z"
      }
    ]
  }
}
```

**Output**: Diagnostic results and recommendations

```json
{
  "diagnostic_results": {
    "error_analysis": {
      "duplicate_inchikey": {
        "count": 1,
        "recommendation": "Use upsert operation with update_existing=true to update existing molecules instead of inserting duplicates",
        "affected_molecules": ["temp_mol_id_78901"]
      },
      "invalid_smiles": {
        "count": 1,
        "recommendation": "Fix unbalanced parentheses in SMILES string: C1CC(C)C1)",
        "affected_molecules": ["temp_mol_id_78902"]
      },
      "property_validation": {
        "count": 1,
        "recommendation": "Check property value range or adjust validation constraints for logP",
        "affected_molecules": ["temp_mol_id_78903"]
      }
    },
    "corrective_actions": [
      {
        "action_type": "update_existing",
        "description": "Switch insertion mode to update existing molecules",
        "affected_count": 1
      },
      {
        "action_type": "fix_smiles",
        "description": "Fix SMILES syntax",
        "affected_count": 1
      },
      {
        "action_type": "validate_property_range",
        "description": "Verify logP calculation or adjust valid range",
        "affected_count": 1
      }
    ]
  },
  "resumption_strategy": {
    "recommendation": "Resume import with corrected data and update_existing=true",
    "molecules_to_retry": ["temp_mol_id_78901", "temp_mol_id_78902", "temp_mol_id_78903"],
    "suggested_batch_size": 50
  },
  "metadata": {
    "timestamp": "2025-05-05T19:15:45Z",
    "diagnostic_id": "DIAG-20250505-191545"
  }
}
```

## Machine Learning Integration

### Molecular Property Prediction

**Input**: Molecules for property prediction

```json
{
  "prediction_request": {
    "molecules": [
      {
        "id": "pred_mol_1",
        "smiles": "CCC(O)CO",
        "name": "1,2-Butanediol"
      },
      {
        "id": "pred_mol_2",
        "smiles": "CC(O)C(O)C",
        "name": "2,3-Butanediol"
      }
    ],
    "properties_to_predict": [
      "glass_transition_temperature",
      "cell_permeability",
      "ice_recrystallization_inhibition",
      "cytotoxicity_ic50"
    ],
    "model_version": "2.3.1",
    "include_uncertainty": true
  }
}
```

**Output**: Predicted properties with confidence intervals

```json
{
  "predictions": [
    {
      "molecule": {
        "id": "pred_mol_1",
        "smiles": "CCC(O)CO",
        "name": "1,2-Butanediol"
      },
      "predicted_properties": {
        "glass_transition_temperature": {
          "value": -115.3,
          "unit": "°C",
          "confidence_interval": [-120.6, -110.0],
          "confidence_level": 0.95,
          "uncertainty": 5.3
        },
        "cell_permeability": {
          "value": 3.8e-6,
          "unit": "cm/s",
          "confidence_interval": [2.5e-6, 5.1e-6],
          "confidence_level": 0.95,
          "uncertainty": 1.3e-6
        },
        "ice_recrystallization_inhibition": {
          "value": 0.42,
          "unit": "relative to PVA",
          "confidence_interval": [0.31, 0.53],
          "confidence_level": 0.95,
          "uncertainty": 0.11
        },
        "cytotoxicity_ic50": {
          "value": 850.2,
          "unit": "μM",
          "confidence_interval": [720.5, 980.0],
          "confidence_level": 0.95,
          "uncertainty": 129.8
        }
      },
      "prediction_quality": {
        "applicability_domain_distance": 0.12,
        "similar_training_compounds": ["Ethylene glycol", "Propylene glycol", "Glycerol"],
        "overall_confidence": "high"
      }
    },
    {
      "molecule": {
        "id": "pred_mol_2",
        "smiles": "CC(O)C(O)C",
        "name": "2,3-Butanediol"
      },
      "predicted_properties": "..."
    }
  ],
  "metadata": {
    "model_version": "2.3.1",
    "prediction_timestamp": "2025-05-06T10:15:33Z",
    "model_performance": {
      "glass_transition_temperature": {
        "rmse": 8.2,
        "r2": 0.86
      },
      "cell_permeability": {
        "rmse": 0.8e-6,
        "r2": 0.79
      },
      "ice_recrystallization_inhibition": {
        "rmse": 0.15,
        "r2": 0.82
      },
      "cytotoxicity_ic50": {
        "rmse": 150.2,
        "r2": 0.75
      }
    }
  }
}
```

## API Integration

### External API Response Format

**Input**: API request for cryoprotectant candidates

```http
GET /api/v1/molecules/cryoprotectants?min_mw=60&max_mw=200&min_score=70
```

**Output**: API response with molecule data for client consumption

```json
{
  "data": {
    "molecules": [
      {
        "id": "550e8400-e29b-41d4-a716-446655440000",
        "name": "Glycerol",
        "synonyms": ["Glycerin", "1,2,3-Propanetriol"],
        "smiles": "C(C(CO)O)O",
        "inchikey": "PEDCQBHIVMGVHV-UHFFFAOYSA-N",
        "formula": "C3H8O3",
        "molecular_weight": 92.09,
        "properties": {
          "logP": -1.76,
          "hydrogen_bond_donors": 3,
          "hydrogen_bond_acceptors": 3,
          "rotatable_bonds": 2,
          "polar_surface_area": 60.7,
          "glass_transition_temperature": -93.0
        },
        "cryoprotectant_score": 89.4,
        "experimental_data_available": true,
        "structure_visualization_url": "https://cryoprotect.example.com/api/v1/molecules/550e8400-e29b-41d4-a716-446655440000/image",
        "external_links": {
          "pubchem": "https://pubchem.ncbi.nlm.nih.gov/compound/753"
        }
      }
    ],
    "pagination": {
      "current_page": 1,
      "total_pages": 5,
      "total_results": 98,
      "page_size": 20
    }
  },
  "links": {
    "self": "/api/v1/molecules/cryoprotectants?min_mw=60&max_mw=200&min_score=70&page=1",
    "next": "/api/v1/molecules/cryoprotectants?min_mw=60&max_mw=200&min_score=70&page=2",
    "first": "/api/v1/molecules/cryoprotectants?min_mw=60&max_mw=200&min_score=70&page=1",
    "last": "/api/v1/molecules/cryoprotectants?min_mw=60&max_mw=200&min_score=70&page=5"
  },
  "meta": {
    "timestamp": "2025-05-06T13:25:44Z",
    "api_version": "1.0",
    "filters_applied": {
      "molecular_weight": "60-200",
      "cryoprotectant_score": "≥70"
    }
  }
}
```

## Logging and Monitoring

### Pipeline Execution Logs

**Output**: Detailed execution logs for data pipeline

```json
{
  "pipeline_execution_log": {
    "execution_id": "pipeline-20250506-120000",
    "start_time": "2025-05-06T12:00:00Z",
    "end_time": "2025-05-06T12:15:22Z",
    "total_duration_seconds": 922,
    "status": "completed",
    "stages": [
      {
        "stage_name": "data_acquisition",
        "start_time": "2025-05-06T12:00:00Z",
        "end_time": "2025-05-06T12:03:45Z",
        "duration_seconds": 225,
        "status": "completed",
        "details": {
          "pubchem_molecules_retrieved": 850,
          "chembl_molecules_retrieved": 125,
          "api_calls_made": 48,
          "rate_limit_pauses": 2
        }
      },
      {
        "stage_name": "structure_standardization",
        "start_time": "2025-05-06T12:03:45Z",
        "end_time": "2025-05-06T12:05:12Z",
        "duration_seconds": 87,
        "status": "completed",
        "details": {
          "molecules_processed": 975,
          "standardization_operations": {
            "salt_removal": 53,
            "charge_neutralization": 41,
            "tautomer_standardization": 102
          },
          "duplicate_inchikeys_found": 87
        }
      },
      {
        "stage_name": "property_calculation",
        "start_time": "2025-05-06T12:05:12Z",
        "end_time": "2025-05-06T12:08:33Z",
        "duration_seconds": 201,
        "status": "completed",
        "details": {
          "molecules_processed": 888,
          "properties_calculated": 8436,
          "calculation_errors": 12,
          "average_calculation_time_ms": 226
        }
      },
      {
        "stage_name": "filtering",
        "start_time": "2025-05-06T12:08:33Z",
        "end_time": "2025-05-06T12:09:15Z",
        "duration_seconds": 42,
        "status": "completed",
        "details": {
          "molecules_evaluated": 888,
          "molecules_passed": 342,
          "molecules_failed": 546,
          "failure_reasons": {
            "molecular_weight": 87,
            "logP": 156,
            "hydrogen_bond_donors": 201,
            "structural_filters": 271
          }
        }
      },
      {
        "stage_name": "similarity_calculation",
        "start_time": "2025-05-06T12:09:15Z",
        "end_time": "2025-05-06T12:11:52Z",
        "duration_seconds": 157,
        "status": "completed",
        "details": {
          "molecules_processed": 342,
          "similarity_methods": ["morgan", "maccs", "rdkit"],
          "reference_molecules": 5,
          "average_similarity_calculation_time_ms": 92
        }
      },
      {
        "stage_name": "scoring_and_ranking",
        "start_time": "2025-05-06T12:11:52Z",
        "end_time": "2025-05-06T12:12:41Z",
        "duration_seconds": 49,
        "status": "completed",
        "details": {
          "molecules_scored": 342,
          "score_distribution": {
            "min": 12.6,
            "max": 89.4,
            "mean": 61.3,
            "p25": 42.8,
            "p50": 65.7,
            "p75": 79.1
          }
        }
      },
      {
        "stage_name": "database_insertion",
        "start_time": "2025-05-06T12:12:41Z",
        "end_time": "2025-05-06T12:15:22Z",
        "duration_seconds": 161,
        "status": "completed",
        "details": {
          "molecules_to_insert": 342,
          "new_molecules_inserted": 255,
          "existing_molecules_updated": 87,
          "insertion_batches": 4,
          "properties_inserted": 3078,
          "database_transactions": 8
        }
      }
    ],
    "warnings": [
      {
        "stage": "property_calculation",
        "message": "12 molecules failed property calculation due to unusual structures",
        "timestamp": "2025-05-06T12:08:30Z"
      },
      {
        "stage": "database_insertion",
        "message": "Database connection pool nearing capacity (85%)",
        "timestamp": "2025-05-06T12:14:15Z"
      }
    ],
    "resources": {
      "peak_memory_usage_mb": 1248,
      "average_cpu_utilization": 68.4,
      "peak_database_connections": 17
    }
  },
  "metadata": {
    "pipeline_version": "3.2.1",
    "environment": "production",
    "trigger": "scheduled",
    "log_id": "LOG-20250506-120000"
  }
}
```