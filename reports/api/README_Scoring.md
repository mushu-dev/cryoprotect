# CryoProtect Analyzer - Cryoprotection Effectiveness Scoring System

This document describes the cryoprotection effectiveness scoring system implemented in the CryoProtect Analyzer project. The scoring system evaluates how well a molecule or mixture might function as a cryoprotectant based on its molecular properties.

## Overview

The scoring system calculates a composite score (0-100) based on multiple molecular properties that are known to be important for cryoprotection. Each property contributes to the overall score with a specific weight based on its importance in cryoprotection effectiveness.

### Scored Properties

The following properties are evaluated:

1. **Hydrogen Bonding Capacity (25%)**: The ability to form hydrogen bonds with water and biomolecules is critical for cryoprotection. Molecules with optimal hydrogen bonding capacity can displace water molecules and prevent ice crystal formation.

2. **LogP (15%)**: The partition coefficient (LogP) measures hydrophobicity/hydrophilicity balance. Effective cryoprotectants typically have moderate hydrophilicity (negative LogP values).

3. **Molecular Size and Weight (15%)**: Smaller molecules can penetrate cells more easily, but extremely small molecules may not provide sufficient cryoprotection. The optimal range is around 60-180 Da.

4. **Topological Polar Surface Area (15%)**: TPSA is related to a molecule's ability to permeate cell membranes. Optimal cryoprotectants have moderate TPSA values (40-90 Å²).

5. **Functional Groups (20%)**: Certain functional groups, especially hydroxyl (-OH) groups, are particularly important for cryoprotection. They can form hydrogen bonds with water and biomolecules.

6. **Permeability (10%)**: The ability to penetrate cell membranes is important for intracellular cryoprotection. This is estimated based on various physicochemical properties.

### Scoring Algorithm

Each property is scored on a scale of 0-1 based on how close it is to the optimal range for cryoprotection. The overall score is a weighted sum of these individual scores, scaled to 0-100.

For mixtures, the system calculates a weighted average of the component scores based on their concentrations, with a small synergy bonus for mixtures with complementary properties.

## API Endpoints

The scoring system is accessible through the following API endpoints:

### 1. Score a Molecule by SMILES/MOL/SDF

**Endpoint:** `POST /api/v1/scoring/molecules`

**Request Body:**
```json
{
  "molecule_data": "C(C(CO)O)O",  // SMILES, MOL, or SDF data
  "input_format": "smiles",        // "smiles", "mol", or "sdf"
  "store_result": false            // Whether to store the result in the database
}
```

**Response:**
```json
{
  "overall_score": 85,
  "component_scores": {
    "hydrogen_bonding": 95,
    "logp": 80,
    "molecular_size": 85,
    "tpsa": 82,
    "functional_groups": 90,
    "permeability": 75
  },
  "properties": {
    // Detailed molecular properties used for scoring
  }
}
```

### 2. Score a Molecule by ID

**Endpoint:** `POST /api/v1/molecules/{molecule_id}/score`

**Request Body:**
```json
{
  "store_result": true  // Whether to store the result in the database
}
```

**Response:** Same as above

### 3. Score a Mixture by ID

**Endpoint:** `POST /api/v1/mixtures/{mixture_id}/score`

**Request Body:**
```json
{
  "store_result": true  // Whether to store the result in the database
}
```

**Response:**
```json
{
  "mixture_id": "123e4567-e89b-12d3-a456-426614174000",
  "name": "Glycerol-DMSO Mixture",
  "overall_score": 82,
  "component_scores": [
    {
      "molecule_id": "123e4567-e89b-12d3-a456-426614174001",
      "name": "Glycerol",
      "concentration": 30.0,
      "concentration_unit": "% w/v",
      "score": 85
    },
    {
      "molecule_id": "123e4567-e89b-12d3-a456-426614174002",
      "name": "DMSO",
      "concentration": 10.0,
      "concentration_unit": "% w/v",
      "score": 75
    }
  ]
}
```

## Scientific Basis

The scoring system is based on scientific literature on cryoprotectants and their properties. Key references include:

1. Fuller, B. J. (2004). Cryoprotectants: the essential antifreezes to protect life in the frozen state. CryoLetters, 25(6), 375-388.

2. Best, B. P. (2015). Cryoprotectant toxicity: facts, issues, and questions. Rejuvenation research, 18(5), 422-436.

3. Elliott, G. D., Wang, S., & Fuller, B. J. (2017). Cryoprotectants: A review of the actions and applications of cryoprotective solutes that modulate cell recovery from ultra-low temperatures. Cryobiology, 76, 74-91.

## Examples of Good Cryoprotectants

The system has been calibrated to give high scores to known effective cryoprotectants:

1. **Glycerol** (SMILES: `C(C(CO)O)O`) - Score: ~85
   - Multiple hydroxyl groups for hydrogen bonding
   - Good hydrophilicity
   - Optimal molecular size

2. **DMSO (Dimethyl sulfoxide)** (SMILES: `CS(=O)C`) - Score: ~75
   - Good hydrogen bond acceptor
   - Excellent cell permeability
   - Moderate size

3. **Ethylene Glycol** (SMILES: `C(CO)O`) - Score: ~80
   - Multiple hydroxyl groups
   - Small size for good permeability
   - Good hydrophilicity

4. **Propylene Glycol** (SMILES: `CC(O)CO`) - Score: ~78
   - Multiple hydroxyl groups
   - Good balance of properties

## Implementation Details

The scoring system is implemented in the `api/scoring.py` module, with the following key functions:

- `calculate_molecule_score`: Calculates the overall score for a molecule
- `calculate_mixture_score`: Calculates the overall score for a mixture
- `store_molecule_score`: Stores a molecule's score in the database
- `store_mixture_score`: Stores a mixture's score in the database

The API endpoints are implemented in the `api/scoring_resources.py` module.

## Testing

The scoring system has been tested with known cryoprotectants and non-cryoprotectants to ensure it correctly differentiates between them. Test cases can be found in `tests/test_scoring.py`.

## Future Improvements

Potential future improvements to the scoring system include:

1. Machine learning-based scoring using experimental data
2. More sophisticated synergy modeling for mixtures
3. Integration with molecular dynamics simulations
4. Consideration of toxicity and biocompatibility
5. Customizable scoring weights for different applications