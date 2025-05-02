# CryoProtect v2 Database Schema

This document provides a comprehensive overview of the CryoProtect v2 database schema, including tables, relationships, constraints, and a visual relationship diagram. It is intended for both developers and scientists to understand the structure and logic of the database.

---

## Table Summaries

### 1. molecules
- Stores basic information about cryoprotectant molecules (from PubChem)
- Fields: id (UUID, PK), pubchem_cid, name, formula, smiles, inchikey, etc.

### 2. property_types
- Defines types of properties that can be measured (e.g., Molecular Weight, LogP)
- Fields: id (UUID, PK), name, data_type

### 3. molecular_properties
- Stores properties of individual molecules
- Fields: id (UUID, PK), molecule_id (FK), property_type_id (FK), value (numeric/text/boolean)

### 4. mixtures
- Information about mixtures of molecules
- Fields: id (UUID, PK), name, description, created_by

### 5. mixture_components
- Which molecules are in which mixtures, with concentrations
- Fields: id (UUID, PK), mixture_id (FK), molecule_id (FK), concentration, concentration_unit

### 6. calculation_methods
- Methods used for property predictions
- Fields: id (UUID, PK), name, description, version

### 7. predictions
- Computational predictions of mixture or molecule properties
- Fields: id (UUID, PK), molecule_id (FK, optional), mixture_id (FK, optional), property_type_id (FK), calculation_method_id (FK), predicted_value, confidence

### 8. experiments
- Experimental data for mixtures
- Fields: id (UUID, PK), mixture_id (FK), name, date_performed, conditions

### 9. experiment_properties
- Properties measured in experiments
- Fields: id (UUID, PK), experiment_id (FK), property_type_id (FK), value, unit

### 10. projects
- Research projects grouping data
- Fields: id (UUID, PK), name, description, team_id (FK)

### 11. teams
- Research teams
- Fields: id (UUID, PK), name, description

### 12. user_profile
- User accounts and metadata
- Fields: id (UUID, PK), user_id, name, email, etc.

---

## Relationships and Constraints

- **molecular_properties.molecule_id** → molecules.id (FK)
- **molecular_properties.property_type_id** → property_types.id (FK)
- **mixture_components.mixture_id** → mixtures.id (FK)
- **mixture_components.molecule_id** → molecules.id (FK)
- **predictions.molecule_id** → molecules.id (FK, optional)
- **predictions.mixture_id** → mixtures.id (FK, optional)
- **predictions.property_type_id** → property_types.id (FK)
- **predictions.calculation_method_id** → calculation_methods.id (FK)
- **experiments.mixture_id** → mixtures.id (FK)
- **experiment_properties.experiment_id** → experiments.id (FK)
- **experiment_properties.property_type_id** → property_types.id (FK)
- **projects.team_id** → teams.id (FK)
- **user_profile.team_id** → teams.id (FK, optional)

**Constraints:**
- No duplicate molecules by inchikey
- Required fields: molecules (inchikey, name, cid), molecular_properties (molecule_id, property_type_id, value), mixture_components (mixture_id, molecule_id, concentration, concentration_unit)
- Logical: Mixture component concentrations (percent) sum to ~100%; no negative concentrations

---

## Entity-Relationship Diagram

```plaintext
+----------------+      +---------------------+      +------------------+
|   molecules    |<-----| molecular_properties|----->|  property_types  |
+----------------+      +---------------------+      +------------------+
       ^                        ^
       |                        |
       |                        |
       |                        |
+----------------+      +---------------------+      +------------------+
|mixture_components|<----|    mixtures        |      | calculation_methods|
+----------------+      +---------------------+      +------------------+
       ^                        ^                        ^
       |                        |                        |
       |                        |                        |
+----------------+      +---------------------+      +------------------+
|  predictions   |<-----|   experiments       |----->|experiment_properties|
+----------------+      +---------------------+      +------------------+
       ^                        ^
       |                        |
+----------------+      +---------------------+
|   projects     |<-----|      teams          |
+----------------+      +---------------------+
       ^
       |
+----------------+
| user_profile   |
+----------------+
```

**Legend:**
- Arrows indicate foreign key relationships (child → parent)
- Some relationships (e.g., predictions) may reference either molecules or mixtures

---

## Views and Functions

- **Views:**
  - `molecule_with_properties`: Molecules with all properties (JSON)
  - `mixture_with_components`: Mixtures with all components (JSON)
- **Functions:**
  - `import_molecule_from_pubchem(cid)`
  - `calculate_mixture_score(mixture_id)`
  - `compare_prediction_with_experiment(prediction_id, experiment_id)`

---

## Row-Level Security (RLS)

- All data viewable by everyone
- Only authenticated users can insert new data
- Only creators can update/delete their own data

---

## Notes

- When adding a new molecule, use the import function for PubChem
- When creating mixtures, add components with concentrations
- Use the views for easy querying
- Use the scoring and comparison functions for analysis