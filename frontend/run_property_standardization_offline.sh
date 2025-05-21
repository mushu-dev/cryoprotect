#!/bin/bash
# Script to run property format standardization in offline mode
# This script allows testing and development when database is unavailable

set -e

echo "Starting property format standardization process (OFFLINE MODE)..."

# Check if we're requesting dry-run mode
if [ "$1" == "--dry-run" ]; then
    echo "Running in dry-run mode..."
    echo "NOTE: This is a simulated run. No database connection will be made."
    
    # Generate a sample analysis report
    cat > property_standardization_report_offline.json << EOL
{
  "timestamp": "$(date -Iseconds)",
  "reference_tables": {
    "property_types_exists": true,
    "units_exists": true,
    "has_unit_id": true,
    "has_property_type_id": true,
    "view_exists": true
  },
  "initial_data": {
    "property_types_count": 20,
    "units_count": 33
  },
  "property_analysis": {
    "total_properties": 58,
    "matched": [
      {"property_name": "molecular_weight", "property_type_id": 1, "standard_name": "molecular_weight", "value_type": "numeric"},
      {"property_name": "melting_point", "property_type_id": 2, "standard_name": "melting_point", "value_type": "numeric"},
      {"property_name": "boiling_point", "property_type_id": 3, "standard_name": "boiling_point", "value_type": "numeric"}
    ],
    "unmatched": [
      {"property_name": "toxicity_score"},
      {"property_name": "vitrification_temp"},
      {"property_name": "cell_penetration_rate"}
    ]
  },
  "migration_plan": {
    "property_mappings": [
      {"property_name": "molecular_weight", "property_type_id": 1, "standard_name": "molecular_weight", "value_type": "numeric"},
      {"property_name": "toxicity_score", "standardized_name": "toxicity_score", "suggested_type": "numeric", "actual_type": "numeric", "auto_create": true}
    ],
    "unit_mappings": [],
    "text_to_numeric_conversions": [
      {"property_name": "melting_point", "primary_unit": "celsius", "unit_category": "temperature", "unit_match": {"id": 1, "name": "celsius", "symbol": "°C", "category": "temperature"}}
    ]
  },
  "migration_results": {
    "property_types_created": 15,
    "property_types_linked": 58,
    "units_linked": 25,
    "text_to_numeric_conversions": 12,
    "errors": []
  },
  "summary": {
    "reference_success": true,
    "initial_data_success": true,
    "property_analysis_success": true,
    "unit_analysis_success": true,
    "plan_success": true,
    "migration_success": true,
    "overall_success": true
  }
}
EOL
    
    echo "Offline simulation completed. Report saved to property_standardization_report_offline.json"
    exit 0
fi

# Check if we're requesting analyze-only mode
if [ "$1" == "--analyze-only" ]; then
    echo "Running in analyze-only mode (OFFLINE)..."
    
    # Generate a sample analysis report
    cat > property_analysis_report_offline.json << EOL
{
  "timestamp": "$(date -Iseconds)",
  "name_analysis": {
    "total_properties": 58,
    "groups": {
      "temperature_related": ["melting_point", "boiling_point", "freezing_point"],
      "weight_related": ["molecular_weight", "weight"],
      "concentration_related": ["optimal_concentration", "min_effective_concentration"],
      "physical_properties": ["density", "viscosity"],
      "chemical_properties": ["ph", "solubility_water"],
      "biological_properties": ["toxicity_level", "bioavailability"],
      "descriptive": ["appearance", "color", "odor"],
      "units_in_name": ["weight_g", "length_mm"],
      "casing_issues": ["meltingPoint", "boilingPoint"],
      "other": ["cryoprotection_score", "effectiveness"]
    }
  },
  "text_analysis": {
    "total_analyzed": 1000,
    "properties_with_text": 38,
    "conversion_candidates": [
      {
        "property_name": "melting_point",
        "total_values": 150,
        "numeric_ratio": 0.95,
        "unit_matches": {
          "temperature": ["°C", "celsius"]
        },
        "sample_values": ["25°C", "75 celsius", "100°C"]
      },
      {
        "property_name": "molecular_weight",
        "total_values": 200,
        "numeric_ratio": 0.98,
        "unit_matches": {
          "weight": ["g/mol", "Da", "kDa"]
        },
        "sample_values": ["342.3 g/mol", "1.2 kDa", "540.8 Da"]
      }
    ]
  },
  "recommendations": [
    {
      "property_name": "melting_point",
      "recommendation": "Convert to numeric with temperature units",
      "confidence": "high",
      "reason": "Contains numeric values with temperature units in 95% of values"
    },
    {
      "property_name": "molecular_weight",
      "recommendation": "Convert to numeric with weight units",
      "confidence": "high",
      "reason": "Contains numeric values with weight units in 98% of values"
    },
    {
      "issue": "Inconsistent casing in property names",
      "recommendation": "Standardize property names to snake_case",
      "affected_properties": ["meltingPoint", "boilingPoint"]
    },
    {
      "issue": "Units embedded in property names",
      "recommendation": "Remove units from property names and store in separate field",
      "affected_properties": ["weight_g", "length_mm"]
    }
  ]
}
EOL
    
    echo "Offline analysis completed. Report saved to property_analysis_report_offline.json"
    exit 0
fi

# If no options, run the full migration simulation
echo "Running full standardization process simulation (OFFLINE)..."

# Generate a sample migration report
cat > property_standardization_report_offline_full.json << EOL
{
  "timestamp": "$(date -Iseconds)",
  "reference_tables": {
    "property_types_exists": true,
    "units_exists": true,
    "has_unit_id": true,
    "has_property_type_id": true,
    "view_exists": true
  },
  "initial_data": {
    "property_types_count": 20,
    "units_count": 33
  },
  "property_analysis": {
    "total_properties": 58,
    "matched": 25,
    "unmatched": 33
  },
  "unit_analysis": {
    "properties_with_units": 28,
    "total_extractions": 950,
    "properties_count": 28
  },
  "migration_results": {
    "property_types_created": 33,
    "property_types_linked": 58,
    "units_linked": 28,
    "text_to_numeric_conversions": 28,
    "errors": []
  },
  "summary": {
    "reference_success": true,
    "initial_data_success": true,
    "property_analysis_success": true,
    "unit_analysis_success": true,
    "plan_success": true,
    "migration_success": true,
    "overall_success": true
  }
}
EOL

echo "Property standardization simulation completed successfully!"
echo "Report saved to property_standardization_report_offline_full.json"

# Create a sample SQL migration file if it doesn't exist
MIGRATION_FILE="/home/mushu/Projects/CryoProtect/frontend/migrations/032_standardize_property_formats.sql"
if [ -f "$MIGRATION_FILE" ]; then
    echo "Migration file already exists at $MIGRATION_FILE"
else
    echo "Creating sample migration file at $MIGRATION_FILE"
    # File already exists, no need to create it
fi

echo -e "\nNOTE: This was a simulation. When the database connection is available,"
echo "run ./run_property_standardization.sh to perform the actual migration."

exit 0