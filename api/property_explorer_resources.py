#!/usr/bin/env python3
"""
Property Explorer API Resources

This module provides API endpoints for the Property Explorer interface,
allowing users to fetch and analyze property data from the database.
"""

def register_property_explorer_resources(api):
    """Register property explorer resources with the API"""
    # Register Property Explorer resources
    api.add_resource(DashboardStatsResource, '/stats/dashboard', endpoint='property_dashboard_stats')
    api.add_resource(PropertyTypeResource, '/properties/types', endpoint='property_types')
    api.add_resource(PropertyDataResource, '/properties/data', '/properties/data/<int:property_id>', endpoint='property_data')
    api.add_resource(PropertyStatisticsResource, '/properties/statistics', '/properties/statistics/<int:property_id>', endpoint='property_statistics')
    api.add_resource(PropertyCorrelationResource, '/properties/correlation', endpoint='property_correlation')
    api.add_resource(PropertyExportResource, '/properties/export', endpoint='property_export')
    
def register_property_explorer_docs(docs):
    """Register property explorer resources with the API documentation"""
    from api.api_docs import register_resource
    
    register_resource(docs, DashboardStatsResource, 'property_dashboard_stats')
    register_resource(docs, PropertyTypeResource, 'property_types') 
    register_resource(docs, PropertyDataResource, 'property_data')
    register_resource(docs, PropertyStatisticsResource, 'property_statistics')
    register_resource(docs, PropertyCorrelationResource, 'property_correlation')
    register_resource(docs, PropertyExportResource, 'property_export')

from flask import g, jsonify, request
from flask_restful import Resource
from marshmallow import Schema, fields
import numpy as np
from statistics import mean, median, stdev
from collections import defaultdict

from api.utils import get_supabase_client, token_required
from api.jwt_auth import jwt_required, get_current_user
from logging_enhanced import get_logger, log_with_context

logger = get_logger(__name__)

class PropertySchema(Schema):
    """Schema for property data validation"""
    id = fields.Int(required=True)
    name = fields.Str(required=True)
    type = fields.Str(required=True)
    unit = fields.Str(allow_none=True)
    description = fields.Str(allow_none=True)


class PropertyStatsSchema(Schema):
    """Schema for property statistics"""
    property_id = fields.Int(required=True)
    property_name = fields.Str(required=True)
    min = fields.Float(allow_none=True)
    max = fields.Float(allow_none=True)
    mean = fields.Float(allow_none=True)
    median = fields.Float(allow_none=True)
    std_dev = fields.Float(allow_none=True)
    count = fields.Int(required=True)
    distribution = fields.List(fields.Dict(), allow_none=True)


class DashboardStatsResource(Resource):
    """Resource for dashboard statistics"""
    @jwt_required
    def get(self):
        """Get dashboard statistics"""
        try:
            supabase = get_supabase_client()
            
            # Get total molecule count
            molecules_response = supabase.from_("molecules").select("id", count="exact").execute()
            total_molecules = molecules_response.count if hasattr(molecules_response, 'count') else 0
            
            # Get property types count
            property_types_response = supabase.from_("property_types").select("id", count="exact").execute()
            total_properties = property_types_response.count if hasattr(property_types_response, 'count') else 0
            
            # Get cryoprotectant count
            cryoprotectants_response = supabase.from_("molecules").select("id", count="exact").eq("is_cryoprotectant", True).execute()
            total_cryoprotectants = cryoprotectants_response.count if hasattr(cryoprotectants_response, 'count') else 0
            
            # Calculate data coverage - this is a simplified calculation
            properties_response = supabase.from_("molecular_properties").select("id", count="exact").execute()
            total_property_values = properties_response.count if hasattr(properties_response, 'count') else 0
            
            # Theoretical maximum would be total_molecules * total_properties
            theoretical_max = total_molecules * total_properties if total_molecules and total_properties else 1
            data_coverage = f"{round(total_property_values / theoretical_max * 100, 1)}%" if theoretical_max > 0 else "0%"
            
            return jsonify({
                "totalMolecules": total_molecules,
                "totalProperties": total_properties,
                "totalCryoprotectants": total_cryoprotectants,
                "dataCoverage": data_coverage
            })
            
        except Exception as e:
            log_with_context(
                logger, 'error',
                f"Error fetching dashboard stats: {str(e)}",
                context={
                    "event_type": "api_error",
                    "endpoint": "dashboard_stats",
                    "error_details": str(e)
                },
                exc_info=True
            )
            return jsonify({"error": "Failed to fetch dashboard statistics", "details": str(e)}), 500


class PropertyTypeResource(Resource):
    """Resource for property types"""
    @jwt_required
    def get(self):
        """Get all property types"""
        try:
            supabase = get_supabase_client()
            response = supabase.from_("property_types").select("*").execute()
            
            if not hasattr(response, 'data'):
                return jsonify({"error": "Failed to fetch property types"}), 500
                
            # Group properties by category
            property_types = defaultdict(list)
            for prop in response.data:
                category = prop.get("category", "Other")
                property_types[category].append({
                    "id": prop.get("id"),
                    "name": prop.get("name"),
                    "description": prop.get("description"),
                    "unit": prop.get("unit"),
                    "data_type": prop.get("data_type")
                })
            
            # Convert to list format
            result = [{"category": category, "properties": props} for category, props in property_types.items()]
            
            return jsonify(result)
            
        except Exception as e:
            log_with_context(
                logger, 'error',
                f"Error fetching property types: {str(e)}",
                context={
                    "event_type": "api_error",
                    "endpoint": "property_types",
                    "error_details": str(e)
                },
                exc_info=True
            )
            return jsonify({"error": "Failed to fetch property types", "details": str(e)}), 500


class PropertyDataResource(Resource):
    """Resource for property data"""
    @jwt_required
    def get(self, property_id=None):
        """Get property data for molecules"""
        try:
            supabase = get_supabase_client()
            
            # Handle single property request
            if property_id:
                # Validate property exists
                prop_response = supabase.from_("property_types").select("*").eq("id", property_id).execute()
                if not hasattr(prop_response, 'data') or len(prop_response.data) == 0:
                    return jsonify({"error": f"Property with ID {property_id} not found"}), 404
                
                # Get property data
                response = supabase.from_("molecular_properties").select(
                    "*, molecules(id, name, smiles)"
                ).eq("property_type_id", property_id).execute()
                
                if not hasattr(response, 'data'):
                    return jsonify({"error": f"Failed to fetch data for property {property_id}"}), 500
                
                return jsonify({
                    "property": prop_response.data[0],
                    "values": response.data
                })
            
            # Handle multi-property request
            else:
                # Get filter parameters
                property_ids = request.args.getlist('property_ids')
                limit = request.args.get('limit', 100, type=int)
                offset = request.args.get('offset', 0, type=int)
                
                if not property_ids:
                    return jsonify({"error": "No property_ids specified"}), 400
                
                # Convert property_ids to integers
                try:
                    property_ids = [int(pid) for pid in property_ids]
                except ValueError:
                    return jsonify({"error": "Invalid property_ids format"}), 400
                
                # Fetch molecules with the requested properties
                query = """
                SELECT 
                    m.id, m.name, m.smiles, m.formula, m.is_cryoprotectant,
                    jsonb_object_agg(pt.id::text, mp.value) as properties
                FROM 
                    molecules m
                INNER JOIN 
                    molecular_properties mp ON m.id = mp.molecule_id
                INNER JOIN 
                    property_types pt ON mp.property_type_id = pt.id
                WHERE 
                    pt.id = ANY($1)
                GROUP BY 
                    m.id, m.name, m.smiles, m.formula, m.is_cryoprotectant
                LIMIT $2 OFFSET $3
                """
                
                response = supabase.rpc(
                    "exec_sql", 
                    {"sql_query": query, "params": [property_ids, limit, offset]}
                ).execute()
                
                if not hasattr(response, 'data'):
                    return jsonify({"error": "Failed to fetch property data"}), 500
                
                # Get property information
                prop_response = supabase.from_("property_types").select("*").in_("id", property_ids).execute()
                properties = prop_response.data if hasattr(prop_response, 'data') else []
                
                return jsonify({
                    "properties": properties,
                    "molecules": response.data,
                    "total": len(response.data),  # In a real implementation, this would be a count query
                    "limit": limit,
                    "offset": offset
                })
                
        except Exception as e:
            log_with_context(
                logger, 'error',
                f"Error fetching property data: {str(e)}",
                context={
                    "event_type": "api_error",
                    "endpoint": "property_data",
                    "property_id": property_id,
                    "error_details": str(e)
                },
                exc_info=True
            )
            return jsonify({"error": "Failed to fetch property data", "details": str(e)}), 500


class PropertyStatisticsResource(Resource):
    """Resource for property statistics"""
    @jwt_required
    def get(self, property_id=None):
        """Get statistics for properties"""
        try:
            supabase = get_supabase_client()
            
            # Get filter parameters
            property_ids = request.args.getlist('property_ids')
            
            # If property_id is provided in URL, use that
            if property_id:
                property_ids = [property_id]
            
            # If no property_ids specified, return error
            if not property_ids:
                return jsonify({"error": "No property_ids specified"}), 400
            
            # Convert property_ids to integers
            try:
                property_ids = [int(pid) for pid in property_ids]
            except ValueError:
                return jsonify({"error": "Invalid property_ids format"}), 400
            
            # Get property information
            prop_response = supabase.from_("property_types").select("*").in_("id", property_ids).execute()
            properties = prop_response.data if hasattr(prop_response, 'data') else []
            
            # Create property ID to name mapping
            property_names = {prop["id"]: prop["name"] for prop in properties}
            
            # Initialize results
            statistics = []
            
            # Calculate statistics for each property
            for prop_id in property_ids:
                # Skip if property doesn't exist
                if prop_id not in property_names:
                    continue
                
                # Get property values
                response = supabase.from_("molecular_properties").select("value").eq("property_type_id", prop_id).execute()
                values = [item["value"] for item in response.data] if hasattr(response, 'data') else []
                
                # Skip if no values
                if not values:
                    statistics.append({
                        "property_id": prop_id,
                        "property_name": property_names[prop_id],
                        "min": None,
                        "max": None,
                        "mean": None,
                        "median": None,
                        "std_dev": None,
                        "count": 0,
                        "distribution": None
                    })
                    continue
                
                # Filter out non-numeric values for statistics
                numeric_values = [float(v) for v in values if isinstance(v, (int, float)) or (isinstance(v, str) and v.replace('.', '', 1).isdigit())]
                
                # Calculate statistics
                if numeric_values:
                    min_val = min(numeric_values)
                    max_val = max(numeric_values)
                    mean_val = mean(numeric_values)
                    median_val = median(numeric_values)
                    std_dev_val = stdev(numeric_values) if len(numeric_values) > 1 else 0
                    
                    # Calculate distribution for histogram (10 bins)
                    bins = 10
                    hist, bin_edges = np.histogram(numeric_values, bins=bins)
                    
                    distribution = [
                        {
                            "bin_start": float(bin_edges[i]),
                            "bin_end": float(bin_edges[i+1]),
                            "count": int(hist[i])
                        }
                        for i in range(bins)
                    ]
                else:
                    min_val = None
                    max_val = None
                    mean_val = None
                    median_val = None
                    std_dev_val = None
                    distribution = None
                
                statistics.append({
                    "property_id": prop_id,
                    "property_name": property_names[prop_id],
                    "min": min_val,
                    "max": max_val,
                    "mean": mean_val,
                    "median": median_val,
                    "std_dev": std_dev_val,
                    "count": len(values),
                    "distribution": distribution
                })
            
            return jsonify(statistics)
            
        except Exception as e:
            log_with_context(
                logger, 'error',
                f"Error calculating property statistics: {str(e)}",
                context={
                    "event_type": "api_error",
                    "endpoint": "property_statistics",
                    "property_ids": property_ids,
                    "error_details": str(e)
                },
                exc_info=True
            )
            return jsonify({"error": "Failed to calculate property statistics", "details": str(e)}), 500


class PropertyCorrelationResource(Resource):
    """Resource for property correlations"""
    @jwt_required
    def post(self):
        """Calculate correlation between properties"""
        try:
            data = request.get_json()
            if not data:
                return jsonify({"error": "No data provided"}), 400
            
            property_ids = data.get('property_ids', [])
            method = data.get('method', 'pearson')
            
            if not property_ids or len(property_ids) < 2:
                return jsonify({"error": "At least two property_ids must be specified"}), 400
            
            # Validate correlation method
            valid_methods = ['pearson', 'spearman', 'kendall']
            if method not in valid_methods:
                return jsonify({"error": f"Invalid correlation method. Must be one of: {', '.join(valid_methods)}"}), 400
            
            supabase = get_supabase_client()
            
            # Get property information
            prop_response = supabase.from_("property_types").select("*").in_("id", property_ids).execute()
            properties = prop_response.data if hasattr(prop_response, 'data') else []
            
            # Create property ID to name mapping
            property_names = {prop["id"]: prop["name"] for prop in properties}
            
            # For each pair of properties, get molecules that have both properties
            correlations = []
            
            for i, prop1_id in enumerate(property_ids):
                row = []
                for j, prop2_id in enumerate(property_ids):
                    # Diagonal is always 1.0
                    if i == j:
                        row.append({
                            "coefficient": 1.0,
                            "p_value": 0.0,
                            "sample_size": 0  # This will be updated later
                        })
                        continue
                    
                    # Get molecules with both properties
                    query = """
                    WITH prop1 AS (
                        SELECT molecule_id, value
                        FROM molecular_properties
                        WHERE property_type_id = $1
                    ),
                    prop2 AS (
                        SELECT molecule_id, value
                        FROM molecular_properties
                        WHERE property_type_id = $2
                    )
                    SELECT 
                        prop1.molecule_id,
                        prop1.value as value1,
                        prop2.value as value2
                    FROM 
                        prop1
                    INNER JOIN 
                        prop2 ON prop1.molecule_id = prop2.molecule_id
                    """
                    
                    response = supabase.rpc(
                        "exec_sql", 
                        {"sql_query": query, "params": [prop1_id, prop2_id]}
                    ).execute()
                    
                    pairs = response.data if hasattr(response, 'data') else []
                    
                    # Calculate correlation
                    if pairs:
                        # Filter out non-numeric values
                        numeric_pairs = [
                            (float(p["value1"]), float(p["value2"]))
                            for p in pairs
                            if (isinstance(p["value1"], (int, float)) or (isinstance(p["value1"], str) and p["value1"].replace('.', '', 1).isdigit())) and
                               (isinstance(p["value2"], (int, float)) or (isinstance(p["value2"], str) and p["value2"].replace('.', '', 1).isdigit()))
                        ]
                        
                        if numeric_pairs:
                            x = [p[0] for p in numeric_pairs]
                            y = [p[1] for p in numeric_pairs]
                            
                            # Calculate correlation coefficient based on method
                            coef = 0.0
                            p_value = 1.0
                            
                            if method == 'pearson':
                                from scipy.stats import pearsonr
                                coef, p_value = pearsonr(x, y)
                            elif method == 'spearman':
                                from scipy.stats import spearmanr
                                coef, p_value = spearmanr(x, y)
                            elif method == 'kendall':
                                from scipy.stats import kendalltau
                                coef, p_value = kendalltau(x, y)
                            
                            row.append({
                                "coefficient": float(coef),
                                "p_value": float(p_value),
                                "sample_size": len(numeric_pairs)
                            })
                            
                            # Update diagonal element's sample size
                            if i < j:  # Only update once
                                correlations[i][i]["sample_size"] = len(numeric_pairs)
                        else:
                            # No numeric pairs
                            row.append({
                                "coefficient": None,
                                "p_value": None,
                                "sample_size": 0
                            })
                    else:
                        # No pairs found
                        row.append({
                            "coefficient": None,
                            "p_value": None,
                            "sample_size": 0
                        })
                
                correlations.append(row)
            
            return jsonify({
                "properties": properties,
                "method": method,
                "correlations": correlations
            })
            
        except Exception as e:
            log_with_context(
                logger, 'error',
                f"Error calculating property correlations: {str(e)}",
                context={
                    "event_type": "api_error",
                    "endpoint": "property_correlations",
                    "error_details": str(e)
                },
                exc_info=True
            )
            return jsonify({"error": "Failed to calculate property correlations", "details": str(e)}), 500


class PropertyExportResource(Resource):
    """Resource for exporting property data"""
    @jwt_required
    def post(self):
        """Export property data"""
        try:
            data = request.get_json()
            if not data:
                return jsonify({"error": "No data provided"}), 400
            
            property_ids = data.get('property_ids', [])
            format = data.get('format', 'csv')
            dataset = data.get('dataset', 'all')
            include_structure = data.get('include_structure', True)
            
            if not property_ids:
                return jsonify({"error": "No property_ids specified"}), 400
            
            # Validate export format
            valid_formats = ['csv', 'excel', 'json', 'sdf']
            if format not in valid_formats:
                return jsonify({"error": f"Invalid export format. Must be one of: {', '.join(valid_formats)}"}), 400
            
            supabase = get_supabase_client()
            
            # Build filters based on dataset
            filters = []
            if dataset == 'cryoprotectants':
                filters.append("m.is_cryoprotectant = true")
            
            # We'd normally implement more sophisticated filtering here
            # For simplicity, we'll just create a basic query
            
            filter_clause = " AND ".join(filters) if filters else "TRUE"
            
            # Build query to get molecule data with properties
            query = f"""
            SELECT 
                m.id, m.name, m.smiles, m.formula, m.is_cryoprotectant,
                jsonb_object_agg(pt.id::text, mp.value) as properties
            FROM 
                molecules m
            INNER JOIN 
                molecular_properties mp ON m.id = mp.molecule_id
            INNER JOIN 
                property_types pt ON mp.property_type_id = pt.id
            WHERE 
                pt.id = ANY($1) AND {filter_clause}
            GROUP BY 
                m.id, m.name, m.smiles, m.formula, m.is_cryoprotectant
            LIMIT 1000
            """
            
            response = supabase.rpc(
                "exec_sql", 
                {"sql_query": query, "params": [property_ids]}
            ).execute()
            
            molecules = response.data if hasattr(response, 'data') else []
            
            # Get property information
            prop_response = supabase.from_("property_types").select("*").in_("id", property_ids).execute()
            properties = prop_response.data if hasattr(prop_response, 'data') else []
            
            # Generate a unique ID for this export
            from datetime import datetime
            import hashlib
            
            timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
            hash_input = f"{timestamp}_{','.join(map(str, property_ids))}_{dataset}_{format}"
            export_id = hashlib.md5(hash_input.encode()).hexdigest()[:10]
            
            # In a real implementation, we would write the exported data to a file
            # and provide a download link or initiate a download
            
            # For this demo, we'll just return a success message with export details
            return jsonify({
                "export_id": export_id,
                "filename": f"cryoprotect_export_{timestamp}.{format}",
                "format": format,
                "dataset": dataset,
                "molecule_count": len(molecules),
                "properties": properties,
                "include_structure": include_structure,
                "download_url": f"/api/v1/exports/{export_id}/download"  # This would be a real URL in production
            })
            
        except Exception as e:
            log_with_context(
                logger, 'error',
                f"Error exporting property data: {str(e)}",
                context={
                    "event_type": "api_error",
                    "endpoint": "property_export",
                    "error_details": str(e)
                },
                exc_info=True
            )
            return jsonify({"error": "Failed to export property data", "details": str(e)}), 500