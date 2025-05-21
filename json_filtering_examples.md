# JSON Examples for Molecule Filtering Operations

This document provides detailed JSON examples for common molecular filtering operations in the CryoProtect system. These examples serve as a reference for API consumers and system developers.

## Basic Filtering Examples

### Simple Property Filter

**Request:**

```json
{
  "filter": {
    "property_filters": [
      {
        "property": "molecular_weight",
        "operator": "less_than",
        "value": 100.0
      }
    ]
  }
}
```

**Response:**

```json
{
  "count": 247,
  "results": [
    {
      "id": "550e8400-e29b-41d4-a716-446655440000",
      "name": "Ethylene glycol",
      "smiles": "C(CO)O",
      "molecular_weight": 62.07,
      "properties": {
        "logP": -1.36,
        "hydrogen_bond_donors": 2,
        "hydrogen_bond_acceptors": 2
      }
    },
    {
      "id": "6ba7b810-9dad-11d1-80b4-00c04fd430c8",
      "name": "Methanol",
      "smiles": "CO",
      "molecular_weight": 32.04,
      "properties": {
        "logP": -0.77,
        "hydrogen_bond_donors": 1,
        "hydrogen_bond_acceptors": 1
      }
    }
  ]
}
```

### Multiple Property Filters

**Request:**

```json
{
  "filter": {
    "property_filters": [
      {
        "property": "molecular_weight",
        "operator": "between",
        "min_value": 60.0,
        "max_value": 150.0
      },
      {
        "property": "logP",
        "operator": "less_than",
        "value": 0.0
      },
      {
        "property": "hydrogen_bond_donors",
        "operator": "greater_than_or_equal",
        "value": 2
      }
    ],
    "combine_with": "AND"
  }
}
```

**Response:**

```json
{
  "count": 58,
  "results": [
    {
      "id": "550e8400-e29b-41d4-a716-446655440000",
      "name": "Glycerol",
      "smiles": "C(C(CO)O)O",
      "molecular_weight": 92.09,
      "properties": {
        "logP": -1.76,
        "hydrogen_bond_donors": 3,
        "hydrogen_bond_acceptors": 3
      }
    },
    {
      "id": "6ba7b810-9dad-11d1-80b4-00c04fd430c8",
      "name": "Propylene glycol",
      "smiles": "CC(O)CO",
      "molecular_weight": 76.09,
      "properties": {
        "logP": -0.92,
        "hydrogen_bond_donors": 2,
        "hydrogen_bond_acceptors": 2
      }
    }
  ]
}
```

### Combined Property and Substructure Filters

**Request:**

```json
{
  "filter": {
    "property_filters": [
      {
        "property": "molecular_weight",
        "operator": "less_than",
        "value": 200.0
      }
    ],
    "substructure_filters": {
      "required": [
        {
          "smarts": "[OX2H]",
          "description": "hydroxyl group",
          "min_count": 2
        }
      ],
      "excluded": [
        {
          "smarts": "c1ccccc1",
          "description": "benzene ring"
        }
      ]
    }
  }
}
```

**Response:**

```json
{
  "count": 183,
  "results": [
    {
      "id": "550e8400-e29b-41d4-a716-446655440000",
      "name": "Glycerol",
      "smiles": "C(C(CO)O)O",
      "molecular_weight": 92.09,
      "substructure_matches": {
        "hydroxyl group": {
          "count": 3,
          "atom_indices": [[4], [5], [6]]
        }
      }
    }
  ]
}
```

## Advanced Filtering Examples

### Similarity-Based Search

**Request:**

```json
{
  "filter": {
    "similarity_search": {
      "reference_molecule": "C(C(CO)O)O",
      "reference_format": "smiles",
      "similarity_cutoff": 0.7,
      "similarity_metric": "tanimoto",
      "fingerprint_type": "morgan",
      "max_results": 100
    }
  }
}
```

**Response:**

```json
{
  "count": 24,
  "results": [
    {
      "id": "550e8400-e29b-41d4-a716-446655440000",
      "name": "Glycerol",
      "smiles": "C(C(CO)O)O",
      "molecular_weight": 92.09,
      "similarity": 1.0
    },
    {
      "id": "6ba7b810-9dad-11d1-80b4-00c04fd430c8",
      "name": "Propylene glycol",
      "smiles": "CC(O)CO",
      "molecular_weight": 76.09,
      "similarity": 0.86
    },
    {
      "id": "7cc8d910-a1a2-b3b4-c5c6-d7d8e9f0a1b2",
      "name": "Ethylene glycol",
      "smiles": "C(CO)O",
      "molecular_weight": 62.07,
      "similarity": 0.82
    }
  ]
}
```

### Combined Property, Substructure, and Similarity Filters

**Request:**

```json
{
  "filter": {
    "property_filters": [
      {
        "property": "molecular_weight",
        "operator": "less_than",
        "value": 150.0
      },
      {
        "property": "logP",
        "operator": "between",
        "min_value": -2.0,
        "max_value": 1.0
      }
    ],
    "substructure_filters": {
      "required": [
        {
          "smarts": "[OX2H]",
          "description": "hydroxyl group"
        }
      ]
    },
    "similarity_search": {
      "reference_molecule": "C(C(CO)O)O",
      "similarity_cutoff": 0.6
    },
    "combine_mode": "all"
  }
}
```

**Response:**

```json
{
  "count": 12,
  "results": [
    {
      "id": "550e8400-e29b-41d4-a716-446655440000",
      "name": "Glycerol",
      "smiles": "C(C(CO)O)O",
      "molecular_weight": 92.09,
      "properties": {
        "logP": -1.76,
        "hydrogen_bond_donors": 3,
        "hydrogen_bond_acceptors": 3
      },
      "similarity": 1.0
    },
    {
      "id": "6ba7b810-9dad-11d1-80b4-00c04fd430c8",
      "name": "Propylene glycol",
      "smiles": "CC(O)CO",
      "molecular_weight": 76.09,
      "properties": {
        "logP": -0.92,
        "hydrogen_bond_donors": 2,
        "hydrogen_bond_acceptors": 2
      },
      "similarity": 0.86
    }
  ]
}
```

### Pharmacophore-Based Search

**Request:**

```json
{
  "filter": {
    "pharmacophore_search": {
      "features": [
        {
          "type": "hydrogen_donor",
          "smarts": "[OX2H,NX3;!H0]",
          "min_count": 2
        },
        {
          "type": "hydrogen_acceptor",
          "smarts": "[O,N;!H]",
          "min_count": 2
        },
        {
          "type": "hydrophobic",
          "smarts": "[#6]~[#6]~[#6]",
          "min_count": 1
        }
      ],
      "distances": [
        {
          "between": ["hydrogen_donor", "hydrogen_acceptor"],
          "min_distance": 2.5,
          "max_distance": 4.5,
          "unit": "angstrom"
        }
      ]
    }
  }
}
```

**Response:**

```json
{
  "count": 38,
  "results": [
    {
      "id": "550e8400-e29b-41d4-a716-446655440000",
      "name": "Glycerol",
      "smiles": "C(C(CO)O)O",
      "molecular_weight": 92.09,
      "pharmacophore_matches": {
        "hydrogen_donor": {
          "count": 3,
          "atom_indices": [[4], [5], [6]]
        },
        "hydrogen_acceptor": {
          "count": 3,
          "atom_indices": [[4], [5], [6]]
        },
        "hydrophobic": {
          "count": 1,
          "atom_indices": [[0, 1, 2]]
        }
      }
    }
  ]
}
```

## Composite Scoring and Ranking Examples

### Cryoprotectant Scoring

**Request:**

```json
{
  "filter": {
    "cryoprotectant_filter": {
      "profile": "standard",
      "min_score": 70.0,
      "max_results": 50
    }
  }
}
```

**Response:**

```json
{
  "count": 32,
  "results": [
    {
      "id": "550e8400-e29b-41d4-a716-446655440000",
      "name": "Glycerol",
      "smiles": "C(C(CO)O)O",
      "molecular_weight": 92.09,
      "cryoprotectant_score": 89.4,
      "score_components": {
        "molecular_weight": 8.5,
        "logP": 14.2,
        "hydrogen_bond_donors": 15.0,
        "hydrogen_bond_acceptors": 15.0,
        "hydroxyl_groups": 15.0,
        "sulfoxide_group": 0.0,
        "similarity": 21.7
      }
    },
    {
      "id": "6ba7b810-9dad-11d1-80b4-00c04fd430c8",
      "name": "DMSO",
      "smiles": "CS(C)=O",
      "molecular_weight": 78.13,
      "cryoprotectant_score": 87.2,
      "score_components": {
        "molecular_weight": 8.9,
        "logP": 13.5,
        "hydrogen_bond_donors": 0.0,
        "hydrogen_bond_acceptors": 15.0,
        "hydroxyl_groups": 0.0,
        "sulfoxide_group": 25.0,
        "similarity": 24.8
      }
    }
  ]
}
```

### Custom Weighted Scoring

**Request:**

```json
{
  "filter": {
    "custom_scoring": {
      "property_weights": {
        "molecular_weight": {
          "target": 100.0,
          "width": 40.0,
          "weight": 0.15
        },
        "logP": {
          "target": -1.0,
          "width": 1.0,
          "weight": 0.25
        },
        "hydrogen_bond_donors": {
          "min": 2,
          "max": 4,
          "weight": 0.15
        },
        "hydrogen_bond_acceptors": {
          "min": 2,
          "max": 6,
          "weight": 0.15
        }
      },
      "structural_bonuses": [
        {
          "feature": "hydroxyl",
          "smarts": "[OX2H]",
          "bonus_per_match": 5.0,
          "max_bonus": 15.0,
          "weight": 0.15
        },
        {
          "feature": "sulfoxide",
          "smarts": "[#16X3](=[OX1])([#6])[#6]",
          "flat_bonus": 15.0,
          "weight": 0.15
        }
      ],
      "min_score": 0.6,
      "normalize_scores": true
    },
    "sort": {
      "field": "custom_score",
      "direction": "desc"
    }
  }
}
```

**Response:**

```json
{
  "count": 68,
  "results": [
    {
      "id": "550e8400-e29b-41d4-a716-446655440000",
      "name": "Glycerol",
      "smiles": "C(C(CO)O)O",
      "molecular_weight": 92.09,
      "custom_score": 0.87,
      "normalized_score": 87.0,
      "score_breakdown": {
        "molecular_weight": {
          "value": 92.09,
          "score": 0.94,
          "contribution": 0.141
        },
        "logP": {
          "value": -1.76,
          "score": 0.71,
          "contribution": 0.178
        },
        "hydrogen_bond_donors": {
          "value": 3,
          "score": 1.0,
          "contribution": 0.15
        },
        "hydrogen_bond_acceptors": {
          "value": 3,
          "score": 1.0,
          "contribution": 0.15
        },
        "structural_bonuses": {
          "hydroxyl": {
            "matches": 3,
            "bonus": 15.0,
            "contribution": 0.15
          },
          "sulfoxide": {
            "matches": 0,
            "bonus": 0.0,
            "contribution": 0.0
          }
        }
      }
    }
  ]
}
```

## Machine Learning Integration Examples

### ML-Enhanced Property Prediction

**Request:**

```json
{
  "filter": {
    "property_filters": [
      {
        "property": "molecular_weight",
        "operator": "less_than",
        "value": 200.0
      }
    ],
    "ml_property_prediction": {
      "properties": ["glass_transition_temperature", "cell_permeability"],
      "filters": [
        {
          "property": "glass_transition_temperature",
          "operator": "greater_than",
          "value": -60.0
        }
      ],
      "include_prediction_details": true
    }
  }
}
```

**Response:**

```json
{
  "count": 12,
  "results": [
    {
      "id": "550e8400-e29b-41d4-a716-446655440000",
      "name": "Compound XYZ",
      "smiles": "CCC(O)(O)CC",
      "molecular_weight": 118.17,
      "predicted_properties": {
        "glass_transition_temperature": {
          "value": -52.3,
          "unit": "°C",
          "confidence_interval": [-58.1, -46.5],
          "confidence": 0.85,
          "method": "graph_neural_network"
        },
        "cell_permeability": {
          "value": 2.8e-6,
          "unit": "cm/s",
          "confidence_interval": [1.9e-6, 3.7e-6],
          "confidence": 0.78,
          "method": "graph_neural_network"
        }
      },
      "prediction_model_info": {
        "model_version": "2.3.1",
        "training_set_size": 1240,
        "feature_importance": [
          {"feature": "logP", "importance": 0.28},
          {"feature": "hydrogen_bond_count", "importance": 0.21}
        ]
      }
    }
  ]
}
```

### ML Cryoprotectant Classification

**Request:**

```json
{
  "filter": {
    "property_filters": [
      {
        "property": "molecular_weight",
        "operator": "between",
        "min_value": 60.0,
        "max_value": 250.0
      }
    ],
    "ml_classification": {
      "model": "cryoprotectant_efficacy",
      "min_probability": 0.7,
      "include_explanation": true
    }
  }
}
```

**Response:**

```json
{
  "count": 28,
  "results": [
    {
      "id": "550e8400-e29b-41d4-a716-446655440000",
      "name": "Glycerol",
      "smiles": "C(C(CO)O)O",
      "molecular_weight": 92.09,
      "ml_classification": {
        "classification": "effective_cryoprotectant",
        "probability": 0.92,
        "confidence": "high",
        "explanation": {
          "key_features": [
            {"feature": "hydroxyl_groups", "contribution": "positive", "importance": 0.35},
            {"feature": "logP", "contribution": "positive", "importance": 0.27},
            {"feature": "molecular_weight", "contribution": "positive", "importance": 0.18}
          ],
          "similar_training_examples": [
            "Ethylene glycol", "Propylene glycol", "Sorbitol"
          ]
        }
      }
    }
  ]
}
```

## Pagination and Sorting Examples

**Request:**

```json
{
  "filter": {
    "property_filters": [
      {
        "property": "molecular_weight",
        "operator": "less_than",
        "value": 300.0
      }
    ]
  },
  "pagination": {
    "page": 2,
    "page_size": 10
  },
  "sort": [
    {
      "field": "molecular_weight",
      "direction": "asc"
    },
    {
      "field": "logP",
      "direction": "desc"
    }
  ]
}
```

**Response:**

```json
{
  "count": 386,
  "results": [
    {
      "id": "550e8400-e29b-41d4-a716-446655440000",
      "name": "Ethylene glycol",
      "smiles": "C(CO)O",
      "molecular_weight": 62.07,
      "properties": {
        "logP": -1.36,
        "hydrogen_bond_donors": 2,
        "hydrogen_bond_acceptors": 2
      }
    },
    "..."
  ],
  "pagination": {
    "page": 2,
    "page_size": 10,
    "total_pages": 39,
    "total_results": 386,
    "next_page": 3,
    "previous_page": 1
  }
}
```

## Export Format Examples

### SDF Export

**Request:**

```json
{
  "filter": {
    "cryoprotectant_filter": {
      "profile": "standard",
      "min_score": 80.0
    }
  },
  "export": {
    "format": "sdf",
    "properties": ["molecular_weight", "logP", "cryoprotectant_score"],
    "include_3d_coordinates": true
  }
}
```

**Response:**

```
Molecule data in SDF format...
```

### CSV Export

**Request:**

```json
{
  "filter": {
    "property_filters": [
      {
        "property": "molecular_weight",
        "operator": "less_than",
        "value": 200.0
      }
    ]
  },
  "export": {
    "format": "csv",
    "properties": ["id", "name", "smiles", "molecular_weight", "logP", "hydrogen_bond_donors", "hydrogen_bond_acceptors"]
  }
}
```

**Response:**

```
id,name,smiles,molecular_weight,logP,hydrogen_bond_donors,hydrogen_bond_acceptors
550e8400-e29b-41d4-a716-446655440000,Glycerol,C(C(CO)O)O,92.09,-1.76,3,3
6ba7b810-9dad-11d1-80b4-00c04fd430c8,Propylene glycol,CC(O)CO,76.09,-0.92,2,2
...
```

## Batch Processing Examples

### Batch Similarity Search

**Request:**

```json
{
  "batch_similarity_search": {
    "query_molecules": [
      {
        "id": "query1",
        "smiles": "C(C(CO)O)O",
        "name": "Glycerol"
      },
      {
        "id": "query2",
        "smiles": "CS(C)=O",
        "name": "DMSO"
      }
    ],
    "similarity_cutoff": 0.7,
    "max_results_per_query": 5
  }
}
```

**Response:**

```json
{
  "batch_results": [
    {
      "query_id": "query1",
      "query_name": "Glycerol",
      "results": [
        {
          "id": "550e8400-e29b-41d4-a716-446655440000",
          "name": "Glycerol",
          "smiles": "C(C(CO)O)O",
          "similarity": 1.0
        },
        {
          "id": "6ba7b810-9dad-11d1-80b4-00c04fd430c8",
          "name": "Propylene glycol",
          "smiles": "CC(O)CO",
          "similarity": 0.86
        }
      ],
      "result_count": 2
    },
    {
      "query_id": "query2",
      "query_name": "DMSO",
      "results": [
        {
          "id": "7cc8d910-a1a2-b3b4-c5c6-d7d8e9f0a1b2",
          "name": "DMSO",
          "smiles": "CS(C)=O",
          "similarity": 1.0
        },
        {
          "id": "8dd9e920-b2b3-c4c5-d6d7-e8e9f0a1b2c3",
          "name": "Acetone",
          "smiles": "CC(=O)C",
          "similarity": 0.72
        }
      ],
      "result_count": 2
    }
  ],
  "metadata": {
    "total_queries": 2,
    "successful_queries": 2,
    "total_results": 4
  }
}
```

## Custom Filter Profile Examples

### Creating a Custom Filter Profile

**Request:**

```json
{
  "create_filter_profile": {
    "profile_id": "user_custom_profile",
    "name": "My Custom Cryoprotectant Filter",
    "description": "Custom filter for specialized cryoprotectant applications",
    "property_filters": [
      {
        "property": "molecular_weight",
        "min": 50.0,
        "max": 180.0,
        "weight": 0.2
      },
      {
        "property": "logP",
        "min": -2.5,
        "max": 0.5,
        "weight": 0.25
      },
      {
        "property": "hydrogen_bond_donors",
        "min": 2,
        "max": 5,
        "weight": 0.15
      },
      {
        "property": "hydrogen_bond_acceptors",
        "min": 2,
        "max": 7,
        "weight": 0.15
      }
    ],
    "structural_filters": {
      "required": [
        {
          "name": "hydroxyl",
          "smarts": "[OX2H]",
          "min_count": 2
        }
      ],
      "excluded": [
        {
          "name": "aldehyde",
          "smarts": "[CX3H][=OX1]"
        }
      ]
    },
    "save": true,
    "is_public": false
  }
}
```

**Response:**

```json
{
  "status": "success",
  "profile": {
    "profile_id": "user_custom_profile",
    "name": "My Custom Cryoprotectant Filter",
    "created_by": "user123",
    "created_at": "2025-05-07T10:15:33Z",
    "is_public": false
  },
  "message": "Custom filter profile created successfully."
}
```

### Using a Custom Filter Profile

**Request:**

```json
{
  "filter": {
    "profile_id": "user_custom_profile",
    "min_score": 75.0
  }
}
```

**Response:**

```json
{
  "count": 18,
  "results": [
    {
      "id": "550e8400-e29b-41d4-a716-446655440000",
      "name": "Glycerol",
      "smiles": "C(C(CO)O)O",
      "molecular_weight": 92.09,
      "custom_profile_score": 89.2
    },
    "..."
  ],
  "profile_info": {
    "profile_id": "user_custom_profile",
    "name": "My Custom Cryoprotectant Filter",
    "created_by": "user123",
    "created_at": "2025-05-07T10:15:33Z"
  }
}
```

## Integration with Experimental Data

### Filter with Experimental Results

**Request:**

```json
{
  "filter": {
    "property_filters": [
      {
        "property": "molecular_weight",
        "operator": "less_than",
        "value": 200.0
      }
    ],
    "experimental_data": {
      "required": true,
      "experiment_type": "cell_viability",
      "cell_type": "CHO-K1",
      "min_viability": 70.0
    }
  }
}
```

**Response:**

```json
{
  "count": 8,
  "results": [
    {
      "id": "550e8400-e29b-41d4-a716-446655440000",
      "name": "Glycerol",
      "smiles": "C(C(CO)O)O",
      "molecular_weight": 92.09,
      "experimental_data": {
        "cell_viability": {
          "value": 87.2,
          "unit": "percent",
          "conditions": {
            "cell_type": "CHO-K1",
            "concentration": "10% w/v",
            "freezing_protocol": "slow_cooling_1C_per_min"
          },
          "experiment_id": "EXP-2025-0503-01",
          "date": "2025-05-03"
        }
      }
    },
    "..."
  ]
}
```