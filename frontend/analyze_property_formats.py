#!/usr/bin/env python3
"""
Script to analyze the current formats and units in the molecular_properties table.
This helps identify inconsistencies before standardization.
"""

import os
import sys
import json
import logging
from datetime import datetime
import re
import argparse
from psycopg2.extras import RealDictCursor
import db_utils

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def check_table_structure():
    """Check the structure of the molecular_properties table."""
    try:
        # Check if molecular_properties table exists
        table_query = """
        SELECT EXISTS (
            SELECT 1 
            FROM information_schema.tables 
            WHERE table_schema = 'public' 
            AND table_name = 'molecular_properties'
        );
        """
        table_exists = db_utils.execute_query(table_query)[0][0]
        
        if not table_exists:
            logger.error("molecular_properties table does not exist")
            return {
                'exists': False,
                'error': 'Table does not exist'
            }
        
        # Get table columns
        columns_query = """
        SELECT 
            column_name, 
            data_type, 
            is_nullable,
            column_default
        FROM 
            information_schema.columns
        WHERE 
            table_schema = 'public' 
            AND table_name = 'molecular_properties'
        ORDER BY 
            ordinal_position;
        """
        
        columns = db_utils.execute_query(columns_query, cursor_factory=RealDictCursor)
        
        return {
            'exists': True,
            'columns': columns
        }
    except Exception as e:
        logger.error(f"Error checking table structure: {e}")
        return {
            'exists': False,
            'error': str(e)
        }

def analyze_property_names():
    """Analyze the property names for patterns and inconsistencies."""
    try:
        query = """
        SELECT DISTINCT 
            property_name
        FROM 
            molecular_properties
        ORDER BY 
            property_name;
        """
        
        property_names = db_utils.execute_query(query)
        
        # Group property names by patterns
        groups = {
            'temperature_related': [],
            'weight_related': [],
            'concentration_related': [],
            'physical_properties': [],
            'chemical_properties': [],
            'biological_properties': [],
            'descriptive': [],
            'units_in_name': [],
            'casing_issues': [],
            'other': []
        }
        
        unit_pattern = re.compile(r'.*_([a-zA-Z]+)$')  # Pattern for units at the end
        
        for name in property_names:
            name = name[0]  # Extract from result tuple
            
            # Check for units in name
            unit_match = unit_pattern.match(name)
            if unit_match:
                groups['units_in_name'].append(name)
            
            # Check for casing issues
            if not name.islower() and not name.isupper():
                has_camel_case = False
                for i in range(1, len(name)):
                    if name[i].isupper() and name[i-1].islower():
                        has_camel_case = True
                        break
                
                if has_camel_case:
                    groups['casing_issues'].append(name)
            
            # Categorize by property type
            if any(term in name.lower() for term in ['temp', 'temperature', 'celsius', 'fahrenheit', 'kelvin']):
                groups['temperature_related'].append(name)
            elif any(term in name.lower() for term in ['weight', 'mass', 'molecular_weight', 'mw', 'dalton']):
                groups['weight_related'].append(name)
            elif any(term in name.lower() for term in ['concentration', 'conc', 'molarity', 'molar', 'ppm', 'percent']):
                groups['concentration_related'].append(name)
            elif any(term in name.lower() for term in ['density', 'viscosity', 'boiling', 'melting', 'flash', 'point']):
                groups['physical_properties'].append(name)
            elif any(term in name.lower() for term in ['ph', 'acidity', 'basicity', 'reactive', 'solubility']):
                groups['chemical_properties'].append(name)
            elif any(term in name.lower() for term in ['toxicity', 'ld50', 'ic50', 'bioavailability']):
                groups['biological_properties'].append(name)
            elif any(term in name.lower() for term in ['color', 'odor', 'appearance', 'description']):
                groups['descriptive'].append(name)
            else:
                groups['other'].append(name)
        
        return {
            'total_properties': len(property_names),
            'groups': groups
        }
    except Exception as e:
        logger.error(f"Error analyzing property names: {e}")
        return {
            'error': str(e)
        }

def analyze_property_values():
    """Analyze the values for patterns, units, and inconsistencies."""
    try:
        # Get count by value type
        type_query = """
        SELECT 
            property_name,
            COUNT(*) FILTER (WHERE numeric_value IS NOT NULL) AS numeric_count,
            COUNT(*) FILTER (WHERE text_value IS NOT NULL) AS text_count,
            COUNT(*) FILTER (WHERE boolean_value IS NOT NULL) AS boolean_count,
            COUNT(*) FILTER (WHERE json_value IS NOT NULL) AS json_count,
            COUNT(*) AS total_count
        FROM 
            molecular_properties
        GROUP BY 
            property_name
        ORDER BY 
            property_name;
        """
        
        type_counts = db_utils.execute_query(type_query, cursor_factory=RealDictCursor)
        
        # Identify properties with mixed value types
        mixed_types = [
            prop for prop in type_counts 
            if sum(1 for c in ['numeric_count', 'text_count', 'boolean_count', 'json_count'] 
                   if prop[c] > 0) > 1
        ]
        
        # Sample values for each property
        samples = {}
        
        for prop in type_counts:
            property_name = prop['property_name']
            
            # Determine the predominant type
            if prop['numeric_count'] > 0:
                value_type = 'numeric'
                value_field = 'numeric_value'
            elif prop['text_count'] > 0:
                value_type = 'text'
                value_field = 'text_value'
            elif prop['boolean_count'] > 0:
                value_type = 'boolean'
                value_field = 'boolean_value'
            elif prop['json_count'] > 0:
                value_type = 'json'
                value_field = 'json_value'
            else:
                continue  # Skip if no values
            
            # Get sample values
            sample_query = f"""
            SELECT 
                {value_field} AS value,
                COUNT(*) AS count
            FROM 
                molecular_properties
            WHERE 
                property_name = %s
                AND {value_field} IS NOT NULL
            GROUP BY 
                {value_field}
            ORDER BY 
                COUNT(*) DESC
            LIMIT 10;
            """
            
            sample_values = db_utils.execute_query(sample_query, (property_name,), cursor_factory=RealDictCursor)
            
            samples[property_name] = {
                'type': value_type,
                'samples': sample_values
            }
            
            # Special analysis for text values that might contain units
            if value_type == 'text':
                # Look for patterns indicating units
                unit_values = []
                for sample in sample_values:
                    value = str(sample['value'])
                    
                    # Check for numeric + unit patterns
                    if re.match(r'[-+]?\d*\.?\d+\s*[a-zA-Z°%]+', value):
                        unit_values.append(value)
                    # Check for range with units
                    elif re.match(r'[-+]?\d*\.?\d+\s*-\s*[-+]?\d*\.?\d+\s*[a-zA-Z°%]+', value):
                        unit_values.append(value)
                
                if unit_values:
                    samples[property_name]['unit_values'] = unit_values
                    samples[property_name]['has_units'] = True
                else:
                    samples[property_name]['has_units'] = False
        
        return {
            'type_counts': type_counts,
            'mixed_types': mixed_types,
            'samples': samples
        }
    except Exception as e:
        logger.error(f"Error analyzing property values: {e}")
        return {
            'error': str(e)
        }

def analyze_units():
    """Analyze units in property values and names."""
    try:
        # Common unit patterns in text values
        unit_patterns = {
            'temperature': [r'°C', r'°F', r'K', r'celsius', r'fahrenheit', r'kelvin'],
            'weight': [r'g', r'kg', r'mg', r'µg', r'ng', r'g/mol', r'Da', r'kDa', r'dalton'],
            'volume': [r'L', r'mL', r'µL', r'l', r'ml', r'µl', r'liter', r'milliliter'],
            'concentration': [r'M', r'mM', r'µM', r'nM', r'pM', r'%', r'ppm', r'ppb', r'mol/L'],
            'pressure': [r'Pa', r'kPa', r'MPa', r'bar', r'atm', r'psi', r'mmHg'],
            'length': [r'mm', r'cm', r'm', r'µm', r'nm', r'Å'],
            'time': [r's', r'min', r'h', r'hr', r'second', r'minute', r'hour', r'day'],
            'energy': [r'J', r'kJ', r'kcal', r'cal', r'eV'],
            'angle': [r'°', r'rad', r'degree']
        }
        
        # Query text values for unit analysis
        text_query = """
        SELECT 
            property_name,
            text_value
        FROM 
            molecular_properties
        WHERE 
            text_value IS NOT NULL
            AND text_value ~ '[0-9]'
            AND LENGTH(text_value) < 100  -- Skip long text fields
        LIMIT 10000;  -- Sample size
        """
        
        text_values = db_utils.execute_query(text_query, cursor_factory=RealDictCursor)
        
        # Analyze the text values for units
        unit_findings = {}
        
        for row in text_values:
            property_name = row['property_name']
            value = row['text_value']
            
            if property_name not in unit_findings:
                unit_findings[property_name] = {
                    'detected_units': {},
                    'sample_values': []
                }
            
            if len(unit_findings[property_name]['sample_values']) < 5:
                unit_findings[property_name]['sample_values'].append(value)
            
            # Check each category of units
            for category, patterns in unit_patterns.items():
                for pattern in patterns:
                    # Look for the unit pattern in the value
                    if re.search(r'[0-9.]+\s*' + re.escape(pattern) + r'\b', value, re.IGNORECASE):
                        if category not in unit_findings[property_name]['detected_units']:
                            unit_findings[property_name]['detected_units'][category] = []
                        
                        if pattern not in unit_findings[property_name]['detected_units'][category]:
                            unit_findings[property_name]['detected_units'][category].append(pattern)
        
        # Filter to properties with detected units
        properties_with_units = {
            name: data for name, data in unit_findings.items() 
            if data['detected_units']
        }
        
        # Analyze for inconsistent units within the same property
        inconsistent_units = {}
        
        for prop_name, data in properties_with_units.items():
            for category, units in data['detected_units'].items():
                if len(units) > 1:
                    if prop_name not in inconsistent_units:
                        inconsistent_units[prop_name] = {}
                    
                    inconsistent_units[prop_name][category] = units
        
        return {
            'properties_with_units': properties_with_units,
            'inconsistent_units': inconsistent_units
        }
    except Exception as e:
        logger.error(f"Error analyzing units: {e}")
        return {
            'error': str(e)
        }

def save_report(structure, names, values, units, filename_prefix="property_analysis_report"):
    """Save the analysis report to a file."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{filename_prefix}_{timestamp}.json"
    
    report = {
        'timestamp': datetime.now().isoformat(),
        'table_structure': structure,
        'property_names': names,
        'property_values': values,
        'unit_analysis': units,
        'recommendations': generate_recommendations(structure, names, values, units)
    }
    
    try:
        with open(filename, 'w') as f:
            json.dump(report, f, indent=2)
        logger.info(f"Report saved to {filename}")
        return filename
    except Exception as e:
        logger.error(f"Error saving report: {e}")
        return None

def generate_recommendations(structure, names, values, units):
    """Generate recommendations based on the analysis."""
    recommendations = {
        'naming_standards': [],
        'value_type_standardization': [],
        'unit_standardization': [],
        'general_improvements': []
    }
    
    # Naming recommendations
    if 'groups' in names and 'casing_issues' in names['groups'] and names['groups']['casing_issues']:
        recommendations['naming_standards'].append({
            'issue': 'Inconsistent casing in property names',
            'recommendation': 'Standardize property names to snake_case (lowercase with underscores)',
            'affected_properties': names['groups']['casing_issues']
        })
    
    if 'groups' in names and 'units_in_name' in names['groups'] and names['groups']['units_in_name']:
        recommendations['naming_standards'].append({
            'issue': 'Units embedded in property names',
            'recommendation': 'Remove units from property names and store them in a separate unit field',
            'affected_properties': names['groups']['units_in_name']
        })
    
    # Value type recommendations
    if 'mixed_types' in values and values['mixed_types']:
        mixed_type_props = [prop['property_name'] for prop in values['mixed_types']]
        recommendations['value_type_standardization'].append({
            'issue': 'Mixed value types for the same property',
            'recommendation': 'Standardize on a single value type for each property',
            'affected_properties': mixed_type_props
        })
    
    # Unit recommendations
    if 'inconsistent_units' in units and units['inconsistent_units']:
        for prop_name, categories in units['inconsistent_units'].items():
            for category, unit_list in categories.items():
                recommendations['unit_standardization'].append({
                    'issue': f'Inconsistent {category} units for property {prop_name}',
                    'recommendation': f'Standardize on a single {category} unit (recommend SI units)',
                    'affected_property': prop_name,
                    'current_units': unit_list
                })
    
    # Text values containing numeric values with units
    text_with_numeric = []
    if 'samples' in values:
        for prop_name, data in values['samples'].items():
            if data['type'] == 'text' and data.get('has_units', False):
                text_with_numeric.append(prop_name)
    
    if text_with_numeric:
        recommendations['value_type_standardization'].append({
            'issue': 'Text values containing numeric data with units',
            'recommendation': 'Convert to numeric values and store units separately',
            'affected_properties': text_with_numeric
        })
    
    # General improvements
    if structure['exists'] and 'columns' in structure:
        has_unit_column = any(col['column_name'] == 'unit' for col in structure['columns'])
        if not has_unit_column:
            recommendations['general_improvements'].append({
                'issue': 'Missing dedicated unit column',
                'recommendation': 'Add a unit column to standardize unit storage',
                'affected_tables': ['molecular_properties']
            })
    
    recommendations['general_improvements'].append({
        'issue': 'Lack of enforcement of value types and units',
        'recommendation': 'Implement CHECK constraints or triggers to enforce data quality rules',
        'affected_tables': ['molecular_properties']
    })
    
    recommendations['general_improvements'].append({
        'issue': 'No metadata about property types',
        'recommendation': 'Create a property_types reference table with metadata about each property',
        'details': 'Include expected units, value types, min/max ranges, and descriptions'
    })
    
    return recommendations

def main():
    """Main function to run the script."""
    parser = argparse.ArgumentParser(description="Analyze property formats and units in the molecular_properties table.")
    args = parser.parse_args()
    
    # Check database connection
    if not db_utils.test_connection():
        logger.error("Database connection failed. Exiting.")
        sys.exit(1)
    
    # Check table structure
    logger.info("Checking molecular_properties table structure...")
    structure = check_table_structure()
    
    if not structure['exists']:
        logger.error("molecular_properties table does not exist. Exiting.")
        sys.exit(1)
    
    # Analyze property names
    logger.info("Analyzing property names...")
    names = analyze_property_names()
    
    # Analyze property values
    logger.info("Analyzing property values...")
    values = analyze_property_values()
    
    # Analyze units
    logger.info("Analyzing units...")
    units = analyze_units()
    
    # Save report
    report_file = save_report(structure, names, values, units)
    
    if report_file:
        logger.info(f"Analysis report saved to {report_file}")
        
        # Print summary
        print("\n\nProperty Analysis Summary:")
        print("==========================")
        
        if 'total_properties' in names:
            print(f"Total distinct properties: {names['total_properties']}")
        
        if 'groups' in names:
            print("\nProperty Categories:")
            for category, props in names['groups'].items():
                if props:
                    print(f"- {category.replace('_', ' ').title()}: {len(props)}")
        
        if 'mixed_types' in values:
            print(f"\nProperties with mixed value types: {len(values['mixed_types'])}")
        
        if 'inconsistent_units' in units:
            print(f"Properties with inconsistent units: {len(units['inconsistent_units'])}")
        
        # Print recommendations summary
        recommendations = generate_recommendations(structure, names, values, units)
        
        print("\nKey Recommendations:")
        for category, recs in recommendations.items():
            if recs:
                print(f"\n{category.replace('_', ' ').title()}:")
                for i, rec in enumerate(recs, 1):
                    print(f"{i}. {rec['recommendation']}")
    
    sys.exit(0)

if __name__ == "__main__":
    main()