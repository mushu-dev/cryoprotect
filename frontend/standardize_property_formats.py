#!/usr/bin/env python3
"""
Script to standardize property data formats and units in the molecular_properties table.
This script applies the migration and helps migrate existing data to the new format.
"""

import os
import sys
import json
import logging
import re
from datetime import datetime
import argparse
from psycopg2.extras import RealDictCursor
import db_utils

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def execute_migration():
    """Execute the migration to standardize property formats."""
    sql_file = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), 
        'migrations', 
        '032_standardize_property_formats.sql'
    )
    
    try:
        if db_utils.execute_sql_file(sql_file):
            logger.info("Successfully executed migration to standardize property formats")
            return True
        else:
            logger.error("Failed to execute migration")
            return False
    except Exception as e:
        logger.error(f"Error executing migration: {e}")
        return False

def check_reference_tables():
    """Check if the reference tables were created successfully."""
    try:
        # Check if property_types table exists
        property_types_query = """
        SELECT EXISTS (
            SELECT 1 
            FROM information_schema.tables 
            WHERE table_schema = 'public' 
            AND table_name = 'property_types'
        );
        """
        property_types_exists = db_utils.execute_query(property_types_query)[0][0]
        
        # Check if units table exists
        units_query = """
        SELECT EXISTS (
            SELECT 1 
            FROM information_schema.tables 
            WHERE table_schema = 'public' 
            AND table_name = 'units'
        );
        """
        units_exists = db_utils.execute_query(units_query)[0][0]
        
        # Check if new columns were added to molecular_properties
        columns_query = """
        SELECT 
            column_name
        FROM 
            information_schema.columns
        WHERE 
            table_schema = 'public'
            AND table_name = 'molecular_properties'
            AND column_name IN ('unit_id', 'property_type_id');
        """
        
        new_columns = db_utils.execute_query(columns_query, cursor_factory=RealDictCursor)
        has_unit_id = any(col['column_name'] == 'unit_id' for col in new_columns)
        has_property_type_id = any(col['column_name'] == 'property_type_id' for col in new_columns)
        
        # Check if view was created
        view_query = """
        SELECT EXISTS (
            SELECT 1 
            FROM information_schema.views 
            WHERE table_schema = 'public' 
            AND table_name = 'standardized_properties'
        );
        """
        view_exists = db_utils.execute_query(view_query)[0][0]
        
        return {
            'property_types_exists': property_types_exists,
            'units_exists': units_exists,
            'has_unit_id': has_unit_id,
            'has_property_type_id': has_property_type_id,
            'view_exists': view_exists
        }
    except Exception as e:
        logger.error(f"Error checking reference tables: {e}")
        return {
            'error': str(e)
        }

def check_initial_data():
    """Check if the initial data was populated in the reference tables."""
    try:
        # Count property types
        property_types_query = "SELECT COUNT(*) FROM property_types;"
        property_types_count = db_utils.execute_query(property_types_query)[0][0]
        
        # Count units
        units_query = "SELECT COUNT(*) FROM units;"
        units_count = db_utils.execute_query(units_query)[0][0]
        
        return {
            'property_types_count': property_types_count,
            'units_count': units_count
        }
    except Exception as e:
        logger.error(f"Error checking initial data: {e}")
        return {
            'error': str(e)
        }

def analyze_existing_properties():
    """Analyze existing properties for migration."""
    try:
        query = """
        SELECT DISTINCT 
            property_name
        FROM 
            molecular_properties
        ORDER BY 
            property_name;
        """
        
        property_names = [row[0] for row in db_utils.execute_query(query)]
        logger.info(f"Found {len(property_names)} distinct property names")
        
        # Check which properties match existing property types
        match_query = """
        SELECT 
            mp.property_name,
            pt.id AS property_type_id,
            pt.name AS standard_name,
            pt.value_type
        FROM 
            (SELECT DISTINCT property_name FROM molecular_properties) mp
        LEFT JOIN 
            property_types pt ON LOWER(mp.property_name) = pt.name
            OR LOWER(mp.property_name) = REPLACE(LOWER(pt.display_name), ' ', '_')
            OR LOWER(REPLACE(mp.property_name, '_', ' ')) = LOWER(pt.display_name);
        """
        
        matches = db_utils.execute_query(match_query, cursor_factory=RealDictCursor)
        
        # Group into matched and unmatched
        matched = [m for m in matches if m['property_type_id'] is not None]
        unmatched = [m for m in matches if m['property_type_id'] is None]
        
        logger.info(f"Found {len(matched)} properties matching existing types")
        logger.info(f"Found {len(unmatched)} properties not matching existing types")
        
        return {
            'total_properties': len(property_names),
            'matched': matched,
            'unmatched': unmatched
        }
    except Exception as e:
        logger.error(f"Error analyzing existing properties: {e}")
        return {
            'error': str(e)
        }

def analyze_existing_units():
    """Analyze property values for unit patterns."""
    try:
        # Use the extract_unit_from_text function on text values
        query = """
        WITH text_properties AS (
            SELECT
                id,
                property_name,
                text_value
            FROM
                molecular_properties
            WHERE
                text_value IS NOT NULL
                AND text_value ~ '[0-9]'  -- Contains numeric
                AND LENGTH(text_value) < 100  -- Not too long
            LIMIT 1000  -- Sample size
        )
        SELECT
            tp.id,
            tp.property_name,
            tp.text_value,
            (public.extract_unit_from_text(tp.text_value)).*
        FROM
            text_properties tp;
        """
        
        extracted = db_utils.execute_query(query, cursor_factory=RealDictCursor)
        
        # Group by property name and unit
        by_property = {}
        
        for row in extracted:
            if row['extracted_unit']:
                prop_name = row['property_name']
                if prop_name not in by_property:
                    by_property[prop_name] = {
                        'property_name': prop_name,
                        'units': {},
                        'examples': []
                    }
                
                unit = row['extracted_unit']
                if unit not in by_property[prop_name]['units']:
                    by_property[prop_name]['units'][unit] = 0
                
                by_property[prop_name]['units'][unit] += 1
                
                if len(by_property[prop_name]['examples']) < 3:
                    by_property[prop_name]['examples'].append({
                        'text_value': row['text_value'],
                        'extracted_value': row['extracted_value'],
                        'extracted_unit': row['extracted_unit']
                    })
        
        properties_with_units = list(by_property.values())
        
        return {
            'properties_with_units': properties_with_units,
            'total_extractions': len(extracted),
            'properties_count': len(properties_with_units)
        }
    except Exception as e:
        logger.error(f"Error analyzing existing units: {e}")
        return {
            'error': str(e)
        }

def create_migration_plan(property_analysis, unit_analysis):
    """Create a plan for migrating existing properties."""
    try:
        plan = {
            'property_mappings': [],
            'unit_mappings': [],
            'text_to_numeric_conversions': []
        }
        
        # Add mappings for matched properties
        for prop in property_analysis.get('matched', []):
            plan['property_mappings'].append({
                'property_name': prop['property_name'],
                'property_type_id': prop['property_type_id'],
                'standard_name': prop['standard_name'],
                'value_type': prop['value_type']
            })
        
        # Create suggested mappings for unmatched properties
        for prop in property_analysis.get('unmatched', []):
            # Standardize the name
            standardized_name = re.sub(r'[^a-z0-9_]', '', prop['property_name'].lower().replace(' ', '_'))
            
            # Guess the value type based on the name
            if any(term in standardized_name for term in ['temperature', 'point', 'weight', 'density', 
                                                         'pressure', 'concentration', 'rate', 'ratio',
                                                         'level', 'score', 'index', 'factor', 'constant']):
                suggested_type = 'numeric'
            elif any(term in standardized_name for term in ['formula', 'name', 'description', 'appearance',
                                                           'color', 'odor', 'smiles', 'inchi']):
                suggested_type = 'text'
            elif any(term in standardized_name for term in ['is_', 'has_', 'can_', 'should_', 'flag', 'active']):
                suggested_type = 'boolean'
            else:
                # Check for numeric terms
                if re.search(r'(temp|weight|point|density|concentration|level|score)', standardized_name):
                    suggested_type = 'numeric'
                else:
                    suggested_type = 'text'
            
            # Get a sample value to confirm the type
            sample_query = f"""
            SELECT
                COUNT(*) FILTER (WHERE numeric_value IS NOT NULL) AS num_count,
                COUNT(*) FILTER (WHERE text_value IS NOT NULL) AS text_count,
                COUNT(*) FILTER (WHERE boolean_value IS NOT NULL) AS bool_count
            FROM
                molecular_properties
            WHERE
                property_name = %s;
            """
            
            sample_counts = db_utils.execute_query(sample_query, (prop['property_name'],), cursor_factory=RealDictCursor)[0]
            
            # Determine actual type based on existing data
            if sample_counts['num_count'] > 0:
                actual_type = 'numeric'
            elif sample_counts['text_count'] > 0:
                actual_type = 'text'
            elif sample_counts['bool_count'] > 0:
                actual_type = 'boolean'
            else:
                actual_type = 'unknown'
            
            # If suggested and actual type differ, prefer actual
            if actual_type != 'unknown' and actual_type != suggested_type:
                suggested_type = actual_type
            
            plan['property_mappings'].append({
                'property_name': prop['property_name'],
                'standardized_name': standardized_name,
                'suggested_type': suggested_type,
                'actual_type': actual_type,
                'auto_create': True
            })
        
        # Add unit mappings based on unit analysis
        unit_patterns = {
            'temperature': ['°c', 'celsius', '°f', 'fahrenheit', 'kelvin', 'k'],
            'mass': ['g', 'kg', 'mg', 'gram', 'kilogram', 'milligram'],
            'molecular_weight': ['da', 'kda', 'dalton', 'molecular weight'],
            'volume': ['l', 'ml', 'µl', 'liter', 'milliliter', 'microliter'],
            'concentration': ['m', 'mm', 'µm', 'nm', 'mol', 'molar', '%', 'percent', 'ppm'],
            'pressure': ['pa', 'kpa', 'mpa', 'bar', 'atm', 'psi', 'mmhg'],
            'length': ['m', 'cm', 'mm', 'µm', 'nm', 'angstrom'],
            'time': ['s', 'min', 'h', 'hr', 'second', 'minute', 'hour', 'day'],
            'energy': ['j', 'kj', 'cal', 'kcal', 'ev'],
            'dimensionless': ['ph', 'unitless', 'ratio', 'score']
        }
        
        for prop_data in unit_analysis.get('properties_with_units', []):
            prop_name = prop_data['property_name']
            units = prop_data['units']
            primary_unit = max(units.items(), key=lambda x: x[1])[0] if units else None
            
            # Skip if no unit found
            if not primary_unit:
                continue
            
            # Determine unit category
            unit_category = None
            for category, patterns in unit_patterns.items():
                if any(pattern in primary_unit.lower() for pattern in patterns):
                    unit_category = category
                    break
            
            # Look up standard unit
            unit_query = """
            SELECT 
                id, 
                name,
                symbol,
                category
            FROM 
                units
            WHERE 
                LOWER(name) = LOWER(%s)
                OR LOWER(symbol) = LOWER(%s)
                OR category = %s;
            """
            
            unit_matches = db_utils.execute_query(
                unit_query, 
                (primary_unit, primary_unit, unit_category), 
                cursor_factory=RealDictCursor
            )
            
            best_unit_match = None
            if unit_matches:
                # First try exact match
                exact_matches = [u for u in unit_matches if u['name'].lower() == primary_unit.lower() or u['symbol'].lower() == primary_unit.lower()]
                if exact_matches:
                    best_unit_match = exact_matches[0]
                # Then try category match
                elif unit_category:
                    category_matches = [u for u in unit_matches if u['category'] == unit_category]
                    if category_matches:
                        best_unit_match = category_matches[0]
                # Fallback to first match
                else:
                    best_unit_match = unit_matches[0]
            
            # Add to text-to-numeric conversions if applicable
            plan['text_to_numeric_conversions'].append({
                'property_name': prop_name,
                'primary_unit': primary_unit,
                'unit_category': unit_category,
                'examples': prop_data['examples'],
                'unit_match': best_unit_match,
                'all_units': list(units.keys())
            })
        
        return plan
    except Exception as e:
        logger.error(f"Error creating migration plan: {e}")
        return {
            'error': str(e)
        }

def apply_migration_plan(plan, dry_run=True):
    """Apply the migration plan to the database."""
    results = {
        'property_types_created': 0,
        'property_types_linked': 0,
        'units_linked': 0,
        'text_to_numeric_conversions': 0,
        'errors': []
    }
    
    try:
        # Step 1: Create missing property types
        for mapping in plan['property_mappings']:
            if mapping.get('auto_create', False):
                if dry_run:
                    logger.info(f"Would create property type: {mapping['standardized_name']}")
                    results['property_types_created'] += 1
                    continue
                
                # Create the property type
                query = """
                INSERT INTO property_types (
                    name,
                    display_name,
                    description,
                    value_type
                )
                VALUES (
                    %s,
                    %s,
                    %s,
                    %s
                )
                RETURNING id;
                """
                
                display_name = ' '.join(word.capitalize() for word in mapping['standardized_name'].split('_'))
                
                try:
                    result = db_utils.execute_query(
                        query,
                        (
                            mapping['standardized_name'],
                            display_name,
                            f"Auto-generated from {mapping['property_name']}",
                            mapping['suggested_type']
                        )
                    )
                    
                    if result:
                        property_type_id = result[0][0]
                        mapping['property_type_id'] = property_type_id
                        results['property_types_created'] += 1
                        logger.info(f"Created property type: {mapping['standardized_name']} (ID: {property_type_id})")
                except Exception as e:
                    error = f"Error creating property type {mapping['standardized_name']}: {e}"
                    logger.error(error)
                    results['errors'].append(error)
        
        # Step 2: Link properties to property types
        for mapping in plan['property_mappings']:
            if 'property_type_id' in mapping:
                if dry_run:
                    logger.info(f"Would link property {mapping['property_name']} to property type ID {mapping['property_type_id']}")
                    results['property_types_linked'] += 1
                    continue
                
                # Link the property to the property type
                query = """
                UPDATE molecular_properties
                SET property_type_id = %s
                WHERE property_name = %s;
                """
                
                try:
                    db_utils.execute_query(
                        query,
                        (mapping['property_type_id'], mapping['property_name']),
                        fetch=False
                    )
                    
                    results['property_types_linked'] += 1
                    logger.info(f"Linked property {mapping['property_name']} to property type ID {mapping['property_type_id']}")
                except Exception as e:
                    error = f"Error linking property {mapping['property_name']}: {e}"
                    logger.error(error)
                    results['errors'].append(error)
        
        # Step 3: Process text-to-numeric conversions with units
        for conversion in plan['text_to_numeric_conversions']:
            if not conversion.get('unit_match'):
                error = f"No unit match found for {conversion['property_name']} with unit {conversion['primary_unit']}"
                logger.warning(error)
                results['errors'].append(error)
                continue
            
            if dry_run:
                logger.info(f"Would convert text values for {conversion['property_name']} to numeric with unit {conversion['unit_match']['name']}")
                results['text_to_numeric_conversions'] += 1
                continue
            
            # Find the corresponding property mapping
            property_mapping = next(
                (m for m in plan['property_mappings'] if m['property_name'] == conversion['property_name']),
                None
            )
            
            if not property_mapping or 'property_type_id' not in property_mapping:
                error = f"No property type found for {conversion['property_name']}"
                logger.warning(error)
                results['errors'].append(error)
                continue
            
            # Update the property type to use the matched unit
            update_property_type_query = """
            UPDATE property_types
            SET 
                value_type = 'numeric',
                default_unit = %s,
                unit_category = %s
            WHERE id = %s;
            """
            
            try:
                db_utils.execute_query(
                    update_property_type_query,
                    (
                        conversion['unit_match']['name'],
                        conversion['unit_match']['category'],
                        property_mapping['property_type_id']
                    ),
                    fetch=False
                )
                
                logger.info(f"Updated property type {property_mapping['property_type_id']} to use unit {conversion['unit_match']['name']}")
            except Exception as e:
                error = f"Error updating property type {property_mapping['property_type_id']}: {e}"
                logger.error(error)
                results['errors'].append(error)
                continue
            
            # Use the extract_unit_from_text function to convert text values to numeric
            conversion_query = """
            WITH extracted AS (
                SELECT
                    id,
                    (public.extract_unit_from_text(text_value)).*
                FROM
                    molecular_properties
                WHERE
                    property_name = %s
                    AND text_value IS NOT NULL
                    AND text_value ~ '[0-9]'
            )
            UPDATE molecular_properties mp
            SET
                numeric_value = e.extracted_value,
                unit_id = %s,
                text_value = NULL
            FROM
                extracted e
            WHERE
                mp.id = e.id
                AND e.extracted_value IS NOT NULL;
            """
            
            try:
                db_utils.execute_query(
                    conversion_query,
                    (conversion['property_name'], conversion['unit_match']['id']),
                    fetch=False
                )
                
                results['text_to_numeric_conversions'] += 1
                logger.info(f"Converted text values for {conversion['property_name']} to numeric with unit {conversion['unit_match']['name']}")
            except Exception as e:
                error = f"Error converting values for {conversion['property_name']}: {e}"
                logger.error(error)
                results['errors'].append(error)
        
        return results
    except Exception as e:
        logger.error(f"Error applying migration plan: {e}")
        results['errors'].append(str(e))
        return results

def save_report(
    reference_tables, 
    initial_data, 
    property_analysis, 
    unit_analysis, 
    plan, 
    migration_results=None, 
    filename_prefix="property_standardization_report"
):
    """Save the standardization report to a file."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{filename_prefix}_{timestamp}.json"
    
    # Calculate success metrics
    reference_success = all(v is True for k, v in reference_tables.items() if k != 'error') if reference_tables and 'error' not in reference_tables else False
    initial_data_success = all(v > 0 for k, v in initial_data.items() if k != 'error') if initial_data and 'error' not in initial_data else False
    property_analysis_success = 'total_properties' in property_analysis and property_analysis['total_properties'] > 0 if property_analysis else False
    unit_analysis_success = 'properties_with_units' in unit_analysis and len(unit_analysis['properties_with_units']) > 0 if unit_analysis else False
    plan_success = all(k in plan for k in ['property_mappings', 'unit_mappings', 'text_to_numeric_conversions']) if plan else False
    migration_success = migration_results is not None and len(migration_results.get('errors', [])) == 0 if migration_results else False
    
    report = {
        'timestamp': datetime.now().isoformat(),
        'reference_tables': reference_tables,
        'initial_data': initial_data,
        'property_analysis': property_analysis,
        'unit_analysis': unit_analysis,
        'migration_plan': plan,
        'migration_results': migration_results,
        'summary': {
            'reference_success': reference_success,
            'initial_data_success': initial_data_success,
            'property_analysis_success': property_analysis_success,
            'unit_analysis_success': unit_analysis_success,
            'plan_success': plan_success,
            'migration_success': migration_success,
            'overall_success': reference_success and initial_data_success and property_analysis_success and unit_analysis_success and plan_success and (migration_success if migration_results else True)
        }
    }
    
    try:
        with open(filename, 'w') as f:
            json.dump(report, f, indent=2)
        logger.info(f"Report saved to {filename}")
        return filename
    except Exception as e:
        logger.error(f"Error saving report: {e}")
        return None

def main():
    """Main function to run the script."""
    parser = argparse.ArgumentParser(description="Standardize property data formats and units in the molecular_properties table.")
    parser.add_argument("--analyze-only", action="store_true", help="Only analyze and create a migration plan without applying it")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be done without making changes")
    args = parser.parse_args()
    
    # Check database connection
    if not db_utils.test_connection():
        logger.error("Database connection failed. Exiting.")
        sys.exit(1)
    
    # Execute the migration
    if not args.analyze_only:
        logger.info("Executing migration to standardize property formats...")
        success = execute_migration()
        
        if not success:
            logger.error("Migration failed. Exiting.")
            sys.exit(1)
    
    # Check reference tables
    logger.info("Checking reference tables...")
    reference_tables = check_reference_tables()
    
    if 'error' in reference_tables:
        logger.error(f"Error checking reference tables: {reference_tables['error']}")
        sys.exit(1)
    
    # Check initial data
    logger.info("Checking initial data...")
    initial_data = check_initial_data()
    
    if 'error' in initial_data:
        logger.error(f"Error checking initial data: {initial_data['error']}")
        sys.exit(1)
    
    # Analyze existing properties
    logger.info("Analyzing existing properties...")
    property_analysis = analyze_existing_properties()
    
    if 'error' in property_analysis:
        logger.error(f"Error analyzing properties: {property_analysis['error']}")
        sys.exit(1)
    
    # Analyze existing units
    logger.info("Analyzing existing units in property values...")
    unit_analysis = analyze_existing_units()
    
    if 'error' in unit_analysis:
        logger.error(f"Error analyzing units: {unit_analysis['error']}")
        sys.exit(1)
    
    # Create migration plan
    logger.info("Creating migration plan...")
    plan = create_migration_plan(property_analysis, unit_analysis)
    
    if 'error' in plan:
        logger.error(f"Error creating migration plan: {plan['error']}")
        sys.exit(1)
    
    # Apply migration plan
    migration_results = None
    if not args.analyze_only:
        logger.info("Applying migration plan...")
        migration_results = apply_migration_plan(plan, dry_run=args.dry_run)
        
        if args.dry_run:
            logger.info("Dry run completed. No changes were made.")
        elif migration_results.get('errors', []):
            logger.warning("Migration completed with errors.")
        else:
            logger.info("Migration completed successfully.")
    
    # Save report
    report_file = save_report(
        reference_tables,
        initial_data,
        property_analysis,
        unit_analysis,
        plan,
        migration_results
    )
    
    if report_file:
        logger.info(f"Report saved to {report_file}")
        
        # Print summary
        print("\n\nProperty Standardization Summary:")
        print("=================================")
        print(f"Property Types Table: {'✓' if reference_tables.get('property_types_exists', False) else '✗'}")
        print(f"Units Table: {'✓' if reference_tables.get('units_exists', False) else '✗'}")
        print(f"Property Types Loaded: {initial_data.get('property_types_count', 0)}")
        print(f"Units Loaded: {initial_data.get('units_count', 0)}")
        print(f"Total Properties Analyzed: {property_analysis.get('total_properties', 0)}")
        print(f"Properties Matched to Types: {len(property_analysis.get('matched', []))}")
        print(f"Properties Needing New Types: {len(property_analysis.get('unmatched', []))}")
        print(f"Properties With Unit Patterns: {unit_analysis.get('properties_count', 0)}")
        
        if migration_results:
            print("\nMigration Results:")
            print(f"Property Types Created: {migration_results.get('property_types_created', 0)}")
            print(f"Properties Linked to Types: {migration_results.get('property_types_linked', 0)}")
            print(f"Text-to-Numeric Conversions: {migration_results.get('text_to_numeric_conversions', 0)}")
            print(f"Errors: {len(migration_results.get('errors', []))}")
            
            if args.dry_run:
                print("\nThis was a dry run. Execute the script without --dry-run to apply the changes.")
        
        # Success or not
        if migration_results and not args.dry_run:
            if migration_results.get('errors', []):
                print("\n⚠️ Migration completed with errors. See the report for details.")
                sys.exit(1)
            else:
                print("\n✅ Migration completed successfully!")
                sys.exit(0)
        else:
            print("\n✓ Analysis completed successfully!")
            sys.exit(0)
    else:
        logger.error("Failed to save report.")
        sys.exit(1)

if __name__ == "__main__":
    main()