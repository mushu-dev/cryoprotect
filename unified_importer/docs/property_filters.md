# Property Filters for Molecular Data Sources

This document describes the property filtering capabilities of the unified molecular importer, with a focus on the PubChem data source.

## Overview

Property filters allow you to selectively import compounds that match specific criteria, such as:

- Compounds matching specific search terms
- Compounds with specific physicochemical properties
- Compounds that fit into predefined categories

Property filters can be configured in the unified importer configuration and are applied during the import process.

## Configuring Property Filters

Property filters are specified in the configuration JSON file under the `sources.pubchem.property_filters` section:

```json
{
  "sources": {
    "pubchem": {
      "property_filters": [
        {
          "name": "cryoprotectants", 
          "description": "Known cryoprotectants and related compounds",
          "terms": ["glycerol", "dmso", "ethylene glycol", "trehalose", "sucrose", "cryoprotect"]
        },
        {
          "name": "small_molecules",
          "description": "Small molecules suitable for cryoprotection",
          "molecular_weight_max": 200,
          "logp_max": 1.0,
          "rotatable_bonds_max": 5
        }
      ]
    }
  }
}
```

## Types of Filters

### Term-Based Filters

Term-based filters search for specific terms in compound names and synonyms:

```json
{
  "name": "filter_name",
  "description": "Filter description",
  "terms": ["term1", "term2", "term3"]
}
```

These filters will match compounds that have any of the specified terms in their name or synonyms.

### Property-Based Filters

Property-based filters use numerical constraints on molecular properties:

```json
{
  "name": "filter_name",
  "description": "Filter description",
  "molecular_weight_min": 100,
  "molecular_weight_max": 500,
  "logp_min": -2.0,
  "logp_max": 5.0
}
```

### Supported Property Constraints

The following property constraints are supported:

| Constraint | Description |
|------------|-------------|
| `molecular_weight_min` | Minimum molecular weight |
| `molecular_weight_max` | Maximum molecular weight |
| `logp_min` | Minimum calculated LogP value |
| `logp_max` | Maximum calculated LogP value |
| `hbond_donors_min` | Minimum number of hydrogen bond donors |
| `hbond_donors_max` | Maximum number of hydrogen bond donors |
| `hbond_acceptors_min` | Minimum number of hydrogen bond acceptors |
| `hbond_acceptors_max` | Maximum number of hydrogen bond acceptors |
| `rotatable_bonds_min` | Minimum number of rotatable bonds |
| `rotatable_bonds_max` | Maximum number of rotatable bonds |
| `tpsa_min` | Minimum topological polar surface area |
| `tpsa_max` | Maximum topological polar surface area |

### Combined Filters

You can combine term-based and property-based criteria in a single filter:

```json
{
  "name": "combined_filter",
  "description": "Filter with both terms and properties",
  "terms": ["term1", "term2"],
  "molecular_weight_max": 300,
  "logp_max": 3.0
}
```

Combined filters will match compounds that:
1. Match ANY of the search terms, AND
2. Match ALL of the property constraints

## Using Property Filters in the API

### Searching by Filter Name

You can search directly for compounds matching a specific filter:

```python
# Search using the 'cryoprotectants' filter
cryoprotectants = await pubchem.search_compounds_by_filter(
    'cryoprotectants', 
    max_results=10
)
```

### Applying Filters to Regular Searches

You can apply all configured filters to regular search queries:

```python
# Search for 'alcohol' with filters applied
# This will only return compounds that match 'alcohol' AND at least one filter
alcohols = await pubchem.search_compounds(
    'alcohol', 
    max_results=10, 
    apply_filters=True
)
```

### Streaming Compounds with Filters

Filters can also be applied when streaming compound identifiers:

```python
# Stream compounds matching 'sugar' with filters applied
async for batch in pubchem.stream_compound_identifiers(
    'sugar', 
    batch_size=10, 
    apply_filters=True
):
    # Process batch of compounds...
    pass
```

## Implementation Details

Property filters in the PubChemDataSource class are implemented using the following methods:

- `_initialize_property_filters()`: Loads and validates filters from configuration
- `_validate_property_filter()`: Validates filter structure and parameters
- `filter_compound_by_properties()`: Checks if a compound matches any configured filter
- `_matches_filter()`: Checks if a compound matches a specific filter
- `_get_properties_for_filter()`: Extracts or fetches properties needed for filtering

The filtering process depends on PubChem's API and may require additional API calls to fetch property data not included in the basic compound information.

## Performance Considerations

- Term-based filters are faster as they often require fewer API calls
- Property-based filters may require additional API calls to fetch property data
- Using property filters will increase the number of API calls and may slow down the import process
- Consider using batch processing and checkpoints for large imports with property filters

## Example

See the [property filters example script](../examples/pubchem_property_filters.py) for a complete demonstration of configuring and using property filters.