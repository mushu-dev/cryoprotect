"""
CryoProtect Analyzer API - Export and Sharing Resources

This module contains API resources for exporting data in various formats and sharing results.
"""

import os
import io
import json
import uuid
import csv
import tempfile
import logging
from datetime import datetime
from typing import Dict, List, Any, Optional, Union, Tuple

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from reportlab.lib.pagesizes import letter, landscape
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, PageBreak
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.units import inch
import xlsxwriter

from flask import request, g, current_app, send_file, url_for, jsonify, render_template_string
from flask_restful import Resource, reqparse, abort
from flask_mail import Mail, Message
from marshmallow import Schema, fields as ma_fields, validate, ValidationError

from api.utils import get_supabase_client, token_required, handle_supabase_error, get_user_id
from api.models import BaseModel
from api.protocol_designer import ProtocolDesigner  # For protocol export

# Configure logging
logger = logging.getLogger(__name__)


# Schemas for request validation
class ExportRequestSchema(Schema):
    """
    Schema for data export requests.

    Supported data_type values:
        - molecules
        - mixtures
        - predictions
        - experiments
        - comparisons
        - protocols  # (NEW) Export stepwise protocols by protocol_id

    """
    format = ma_fields.String(required=True, validate=validate.OneOf(['csv', 'json', 'excel', 'pdf']))
    data_type = ma_fields.String(required=True, validate=validate.OneOf(['molecules', 'mixtures', 'predictions', 'experiments', 'comparisons', 'protocols']))
    id = ma_fields.String()  # Optional ID for specific item export
    include_related = ma_fields.Boolean(default=False)  # Include related data
    include_metadata = ma_fields.Boolean(default=True)  # Include metadata (creation date, user, etc.)
    filter_criteria = ma_fields.Dict(default={})  # Optional filter criteria


class VisualizationExportSchema(Schema):
    """Schema for visualization export requests."""
    format = ma_fields.String(required=True, validate=validate.OneOf(['png', 'svg', 'pdf']))
    chart_type = ma_fields.String(required=True, validate=validate.OneOf([
        'property_comparison', 'mixture_composition', 'prediction_accuracy',
        'custom', 'time_series', 'heatmap', 'radar', 'bubble'
    ]))
    data = ma_fields.Dict()  # Custom data for visualization
    width = ma_fields.Integer(default=800)
    height = ma_fields.Integer(default=600)
    title = ma_fields.String()
    style = ma_fields.String(default='default')
    dpi = ma_fields.Integer(default=100)  # Resolution for raster formats
    transparent = ma_fields.Boolean(default=False)  # Transparent background
    include_legend = ma_fields.Boolean(default=True)  # Include legend
    color_scheme = ma_fields.String(default='default')  # Color scheme


class ReportGenerationSchema(Schema):
    """Schema for report generation requests."""
    title = ma_fields.String(required=True)
    sections = ma_fields.List(ma_fields.Dict(), required=True)
    include_visualizations = ma_fields.Boolean(default=True)
    include_data_tables = ma_fields.Boolean(default=True)
    template_id = ma_fields.String()  # Optional template ID


class ShareRequestSchema(Schema):
    """Schema for sharing requests."""
    share_type = ma_fields.String(required=True, validate=validate.OneOf(['link', 'email', 'embed']))
    data_type = ma_fields.String(required=True, validate=validate.OneOf(['molecules', 'mixtures', 'predictions', 'experiments', 'comparisons', 'report', 'visualization', 'protocols']))
    id = ma_fields.String(required=True)  # ID of the item to share
    recipients = ma_fields.List(ma_fields.String())  # Email recipients if share_type is 'email'
    message = ma_fields.String()  # Optional message for email
    password_protected = ma_fields.Boolean(default=False)  # For link sharing
    password = ma_fields.String()  # Optional password for protected links
    expiration = ma_fields.Integer()  # Optional expiration time in seconds
    permissions = ma_fields.String(default='read', validate=validate.OneOf(['read', 'comment', 'edit']))  # Access permissions
    notify_on_access = ma_fields.Boolean(default=False)  # Notify owner when shared item is accessed
    track_analytics = ma_fields.Boolean(default=True)  # Track access analytics


# Helper functions for data export
def get_data_for_export(
    data_type: str,
    item_id: Optional[str] = None,
    include_related: bool = False,
    include_metadata: bool = True,
    filter_criteria: Dict[str, Any] = None
) -> List[Dict[str, Any]]:
    """
    Get data for export based on data type and optional ID.
    
    Args:
        data_type: Type of data to export (molecules, mixtures, etc.)
        item_id: Optional ID for specific item export
        include_related: Whether to include related data
        include_metadata: Whether to include metadata (creation date, user, etc.)
        filter_criteria: Optional filter criteria for the query
        
    Returns:
        List of data dictionaries for export
    """
    supabase = get_supabase_client()
    filter_criteria = filter_criteria or {}
    
    try:
        if data_type == 'molecules':
            query = supabase.table("molecule_with_properties").select("*")
            
            # Apply filters if provided
            for key, value in filter_criteria.items():
                if key == 'name_contains':
                    query = query.ilike('name', f'%{value}%')
                elif key == 'formula_contains':
                    query = query.ilike('molecular_formula', f'%{value}%')
                elif key == 'smiles_contains':
                    query = query.ilike('smiles', f'%{value}%')
                elif key == 'min_molecular_weight':
                    query = query.gte('molecular_weight', value)
                elif key == 'max_molecular_weight':
                    query = query.lte('molecular_weight', value)
                elif key == 'property_exists':
                    # This is a simplified approach - in a real implementation,
                    # you might need a more complex query for JSON properties
                    query = query.not_.is_('properties', None)
            
            if item_id:
                query = query.eq("id", item_id)
                
            response = query.execute()
            
        elif data_type == 'mixtures':
            query = supabase.table("mixture_with_components").select("*")
            
            # Apply filters if provided
            for key, value in filter_criteria.items():
                if key == 'name_contains':
                    query = query.ilike('name', f'%{value}%')
                elif key == 'description_contains':
                    query = query.ilike('description', f'%{value}%')
                elif key == 'component_contains':
                    # This is a simplified approach - in a real implementation,
                    # you might need a more complex query for components
                    query = query.ilike('components', f'%{value}%')
            
            if item_id:
                query = query.eq("id", item_id)
                
            response = query.execute()
            
            # Include related predictions and experiments if requested
            if include_related and item_id and response.data:
                # Get predictions
                pred_query = supabase.from_("predictions").select("""
                    id, mixture_id, property_type_id, numeric_value, text_value, boolean_value,
                    confidence, created_at,
                    property_types(name), calculation_methods(name)
                """).eq("mixture_id", item_id)
                
                # Apply filters to predictions if provided
                if 'prediction_min_confidence' in filter_criteria:
                    pred_query = pred_query.gte('confidence', filter_criteria['prediction_min_confidence'])
                
                pred_response = pred_query.execute()
                
                # Get experiments
                exp_query = supabase.from_("experiments").select("""
                    id, mixture_id, property_type_id, numeric_value, text_value, boolean_value,
                    experimental_conditions, date_performed, created_at,
                    property_types(name)
                """).eq("mixture_id", item_id)
                
                # Apply filters to experiments if provided
                if 'experiment_date_after' in filter_criteria:
                    exp_query = exp_query.gte('date_performed', filter_criteria['experiment_date_after'])
                if 'experiment_date_before' in filter_criteria:
                    exp_query = exp_query.lte('date_performed', filter_criteria['experiment_date_before'])
                
                exp_response = exp_query.execute()
                
                # Add related data to the response
                response.data[0]['predictions'] = pred_response.data if not pred_response.error else []
                response.data[0]['experiments'] = exp_response.data if not exp_response.error else []
                
        elif data_type == 'predictions':
            query = supabase.from_("predictions").select("""
                id, mixture_id, property_type_id, numeric_value, text_value, boolean_value,
                confidence, created_at,
                property_types(name), calculation_methods(name)
            """)
            
            # Apply filters if provided
            for key, value in filter_criteria.items():
                if key == 'min_confidence':
                    query = query.gte('confidence', value)
                elif key == 'property_type':
                    query = query.eq('property_type_id', value)
                elif key == 'calculation_method':
                    query = query.eq('calculation_method_id', value)
            
            if item_id:
                # For predictions, item_id is the mixture_id
                query = query.eq("mixture_id", item_id)
                
            response = query.execute()
            
        elif data_type == 'experiments':
            query = supabase.from_("experiments").select("""
                id, mixture_id, property_type_id, numeric_value, text_value, boolean_value,
                experimental_conditions, date_performed, created_at,
                property_types(name)
            """)
            
            # Apply filters if provided
            for key, value in filter_criteria.items():
                if key == 'date_after':
                    query = query.gte('date_performed', value)
                elif key == 'date_before':
                    query = query.lte('date_performed', value)
                elif key == 'property_type':
                    query = query.eq('property_type_id', value)
                elif key == 'conditions_contain':
                    query = query.ilike('experimental_conditions', f'%{value}%')
            
            if item_id:
                # For experiments, item_id is the mixture_id
                query = query.eq("mixture_id", item_id)
                
            response = query.execute()
            
        elif data_type == 'comparisons':
            if not item_id:
                abort(400, message="Mixture ID is required for comparisons export")
                
            # Get property types for comparison
            prop_response = supabase.table("property_types").select("*").execute()
            if prop_response.error:
                abort(500, message=f"Error fetching property types: {prop_response.error}")
                
            property_types = prop_response.data
            comparisons = []
            
            # Apply property type filter if provided
            if 'property_types' in filter_criteria and filter_criteria['property_types']:
                property_types = [pt for pt in property_types if pt['name'] in filter_criteria['property_types']]
            
            for prop in property_types:
                # Call the comparison function for each property
                comp_response = supabase.rpc(
                    "compare_prediction_with_experiment",
                    {"p_mixture_id": item_id, "p_property_name": prop["name"]}
                ).execute()
                
                if not comp_response.error and comp_response.data:
                    # Apply error threshold filter if provided
                    if 'max_percent_error' in filter_criteria:
                        if abs(comp_response.data.get('percent_error', 0)) <= filter_criteria['max_percent_error']:
                            comparisons.append(comp_response.data)
                    else:
                        comparisons.append(comp_response.data)
                    
            return comparisons
            
        elif data_type == 'protocols':
            # Export a protocol by protocol_id using ProtocolDesigner
            if not item_id:
                abort(400, message="Protocol ID is required for protocols export")
            
            # NOTE: ProtocolDesigner.get_saved_protocol is currently a stub.
            # When persistent storage is implemented, this should fetch from the database.
            protocol = ProtocolDesigner.get_saved_protocol(item_id)
            if not protocol or "error" in protocol:
                abort(404, message=f"Protocol not found for ID: {item_id}")
            
            # Apply filters if provided
            if filter_criteria and protocol:
                if 'step_contains' in filter_criteria:
                    # Filter protocol steps that contain the specified text
                    if 'steps' in protocol:
                        protocol['steps'] = [
                            step for step in protocol['steps']
                            if filter_criteria['step_contains'].lower() in json.dumps(step).lower()
                        ]
            
            # Return as a list for export compatibility
            return [protocol]
            
        else:
            abort(400, message=f"Invalid data type: {data_type}")
            
        # Process response for all data types except comparisons and protocols
        # which have already returned
        error_message, status_code = handle_supabase_error(response)
        if error_message:
            abort(status_code, message=error_message)
            
        # Add metadata if requested
        if include_metadata and response.data:
            for item in response.data:
                item['export_metadata'] = {
                    'exported_at': datetime.now().isoformat(),
                    'exported_by': get_user_id(),
                    'data_type': data_type,
                    'filters_applied': filter_criteria
                }
                
        return response.data or []
        
    except Exception as e:
        logger.error(f"Error exporting data: {str(e)}")
        abort(500, message=f"Error exporting data: {str(e)}")
    
    error_message, status_code = handle_supabase_error(response)
    if error_message:
        abort(status_code, message=error_message)
        
    return response.data


def generate_csv(data: List[Dict[str, Any]], filename: str) -> io.StringIO:
    """Generate CSV from data."""
    if not data:
        return io.StringIO("No data available")
        
    # Convert to DataFrame for easier handling
    df = pd.DataFrame(data)
    
    # Handle nested JSON columns by converting to string
    for col in df.columns:
        if df[col].apply(lambda x: isinstance(x, (dict, list))).any():
            df[col] = df[col].apply(lambda x: json.dumps(x) if isinstance(x, (dict, list)) else x)
    
    output = io.StringIO()
    df.to_csv(output, index=False)
    output.seek(0)
    return output


def generate_excel(data: List[Dict[str, Any]], filename: str) -> io.BytesIO:
    """Generate Excel from data."""
    if not data:
        # Create empty DataFrame with message
        df = pd.DataFrame({"Message": ["No data available"]})
    else:
        # Convert to DataFrame
        df = pd.DataFrame(data)
        
        # Handle nested JSON columns by creating separate sheets
        nested_columns = []
        for col in df.columns:
            if df[col].apply(lambda x: isinstance(x, (dict, list))).any():
                nested_columns.append(col)
        
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        # Write main data
        main_df = df.drop(columns=nested_columns, errors='ignore')
        main_df.to_excel(writer, sheet_name='Main Data', index=False)
        
        # Write nested data to separate sheets
        for col in nested_columns:
            if col in df.columns:
                # Create a new sheet for each nested column
                sheet_name = col[:31]  # Excel sheet names limited to 31 chars
                nested_data = []
                
                for idx, row in df.iterrows():
                    if isinstance(row[col], list):
                        for item in row[col]:
                            item_with_parent = item.copy() if isinstance(item, dict) else {"value": item}
                            item_with_parent["parent_id"] = row.get("id", idx)
                            nested_data.append(item_with_parent)
                    elif isinstance(row[col], dict):
                        item_with_parent = row[col].copy()
                        item_with_parent["parent_id"] = row.get("id", idx)
                        nested_data.append(item_with_parent)
                
                if nested_data:
                    nested_df = pd.DataFrame(nested_data)
                    nested_df.to_excel(writer, sheet_name=sheet_name, index=False)
    
    output.seek(0)
    return output


def generate_json(data: List[Dict[str, Any]], filename: str) -> io.StringIO:
    """
    Generate JSON from data with enhanced formatting and metadata.
    
    Args:
        data: List of data dictionaries to export
        filename: Base filename (without extension)
        
    Returns:
        StringIO object containing JSON data
    """
    try:
        # Create a wrapper object with metadata
        export_data = {
            "metadata": {
                "export_date": datetime.now().isoformat(),
                "record_count": len(data),
                "data_type": filename.split('_')[0] if '_' in filename else "data",
                "generator": "CryoProtect Analyzer"
            },
            "data": data
        }
        
        # Generate JSON with proper formatting
        output = io.StringIO()
        json.dump(export_data, output, indent=2, default=str, ensure_ascii=False)
        output.seek(0)
        
        # Log export
        logger.info(f"Generated JSON export with {len(data)} records")
        
        return output
    except Exception as e:
        logger.error(f"Error generating JSON: {str(e)}")
        raise


def generate_pdf(data: List[Dict[str, Any]], filename: str, data_type: str) -> io.BytesIO:
    """
    Generate PDF report from data with enhanced formatting and features.
    
    Args:
        data: List of data dictionaries to export
        filename: Base filename (without extension)
        data_type: Type of data being exported
        
    Returns:
        BytesIO object containing PDF data
    """
    try:
        buffer = io.BytesIO()
        
        # Use landscape for wider tables
        doc = SimpleDocTemplate(
            buffer,
            pagesize=landscape(letter),
            leftMargin=0.5*inch,
            rightMargin=0.5*inch,
            topMargin=0.5*inch,
            bottomMargin=0.5*inch
        )
        
        # Get styles and create custom styles
        styles = getSampleStyleSheet()
        title_style = ParagraphStyle(
            'CustomTitle',
            parent=styles['Heading1'],
            fontSize=16,
            alignment=1,  # Center alignment
            spaceAfter=12
        )
        subtitle_style = ParagraphStyle(
            'CustomSubtitle',
            parent=styles['Heading2'],
            fontSize=14,
            textColor=colors.navy,
            spaceAfter=6
        )
        section_style = ParagraphStyle(
            'SectionTitle',
            parent=styles['Heading3'],
            fontSize=12,
            textColor=colors.darkblue,
            spaceAfter=6
        )
        normal_style = ParagraphStyle(
            'CustomNormal',
            parent=styles['Normal'],
            fontSize=10,
            spaceAfter=6
        )
        
        elements = []
        
        # Add title with logo if available
        title = f"CryoProtect Analyzer - {data_type.capitalize()} Report"
        elements.append(Paragraph(title, title_style))
        
        # Add timestamp and metadata
        date_text = f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
        elements.append(Paragraph(date_text, normal_style))
        elements.append(Spacer(1, 12))
        
        # Add summary information
        summary_text = f"This report contains data for {len(data)} {data_type}."
        elements.append(Paragraph(summary_text, normal_style))
        elements.append(Spacer(1, 12))
        
        if not data:
            elements.append(Paragraph("No data available", normal_style))
        else:
            # For each item in data, create a section
            for i, item in enumerate(data):
                if i > 0:
                    elements.append(PageBreak())  # Start each item on a new page
                    
                # Item header
                item_name = item.get('name', f"Item {i+1}")
                elements.append(Paragraph(item_name, subtitle_style))
                elements.append(Spacer(1, 6))
                
                # Convert item to table data
                table_data = []
                
                # Handle different data types
                if data_type == 'molecules':
                    # Add basic molecule info
                    table_data.append(["Property", "Value"])
                    table_data.append(["ID", item.get('id', '')])
                    table_data.append(["Name", item.get('name', '')])
                    table_data.append(["Formula", item.get('molecular_formula', '')])
                    table_data.append(["SMILES", item.get('smiles', '')])
                    table_data.append(["PubChem CID", str(item.get('cid', ''))])
                    
                    # Add properties if available
                    if 'properties' in item and item['properties']:
                        elements.append(Spacer(1, 10))
                        elements.append(Paragraph("Molecular Properties", section_style))
                        elements.append(Spacer(1, 6))
                        
                        prop_table = []
                        prop_table.append(["Property", "Value"])
                        
                        for prop_name, prop_value in item['properties'].items():
                            if isinstance(prop_value, (dict, list)):
                                prop_value = json.dumps(prop_value)
                            prop_table.append([prop_name, str(prop_value)])
                        
                        if len(prop_table) > 1:
                            t = Table(prop_table, colWidths=[200, 300])
                            t.setStyle(TableStyle([
                                ('BACKGROUND', (0, 0), (1, 0), colors.darkblue),
                                ('TEXTCOLOR', (0, 0), (1, 0), colors.white),
                                ('ALIGN', (0, 0), (1, 0), 'CENTER'),
                                ('FONTNAME', (0, 0), (1, 0), 'Helvetica-Bold'),
                                ('BOTTOMPADDING', (0, 0), (1, 0), 12),
                                ('BACKGROUND', (0, 1), (-1, -1), colors.lightgrey),
                                ('GRID', (0, 0), (-1, -1), 1, colors.black),
                                ('VALIGN', (0, 0), (-1, -1), 'MIDDLE')
                            ]))
                            elements.append(t)
                
                elif data_type == 'mixtures':
                    # Add basic mixture info
                    table_data.append(["Property", "Value"])
                    table_data.append(["ID", item.get('id', '')])
                    table_data.append(["Name", item.get('name', '')])
                    table_data.append(["Description", item.get('description', '')])
                    
                    # Add components if available
                    if 'components' in item and item['components']:
                        elements.append(Spacer(1, 10))
                        elements.append(Paragraph("Components", section_style))
                        elements.append(Spacer(1, 6))
                        
                        comp_table = []
                        comp_table.append(["Molecule", "Concentration", "Unit"])
                        
                        for comp in item['components']:
                            comp_table.append([
                                comp.get('molecule_name', ''),
                                str(comp.get('concentration', '')),
                                comp.get('concentration_unit', '')
                            ])
                        
                        if len(comp_table) > 1:
                            t = Table(comp_table, colWidths=[200, 150, 150])
                            t.setStyle(TableStyle([
                                ('BACKGROUND', (0, 0), (2, 0), colors.darkblue),
                                ('TEXTCOLOR', (0, 0), (2, 0), colors.white),
                                ('ALIGN', (0, 0), (2, 0), 'CENTER'),
                                ('FONTNAME', (0, 0), (2, 0), 'Helvetica-Bold'),
                                ('BOTTOMPADDING', (0, 0), (2, 0), 12),
                                ('BACKGROUND', (0, 1), (-1, -1), colors.lightgrey),
                                ('GRID', (0, 0), (-1, -1), 1, colors.black),
                                ('VALIGN', (0, 0), (-1, -1), 'MIDDLE')
                            ]))
                            elements.append(t)
                    
                    # Add predictions if available
                    if 'predictions' in item and item['predictions']:
                        elements.append(Spacer(1, 10))
                        elements.append(Paragraph("Predictions", section_style))
                        elements.append(Spacer(1, 6))
                        
                        pred_table = []
                        pred_table.append(["Property", "Value", "Confidence", "Method"])
                        
                        for pred in item['predictions']:
                            # Determine the value based on data type
                            if pred.get('numeric_value') is not None:
                                value = str(pred.get('numeric_value'))
                            elif pred.get('text_value') is not None:
                                value = pred.get('text_value')
                            elif pred.get('boolean_value') is not None:
                                value = str(pred.get('boolean_value'))
                            else:
                                value = "N/A"
                                
                            pred_table.append([
                                pred.get('property_types', {}).get('name', ''),
                                value,
                                str(pred.get('confidence', '')),
                                pred.get('calculation_methods', {}).get('name', '')
                            ])
                        
                        if len(pred_table) > 1:
                            t = Table(pred_table, colWidths=[150, 150, 100, 100])
                            t.setStyle(TableStyle([
                                ('BACKGROUND', (0, 0), (3, 0), colors.darkblue),
                                ('TEXTCOLOR', (0, 0), (3, 0), colors.white),
                                ('ALIGN', (0, 0), (3, 0), 'CENTER'),
                                ('FONTNAME', (0, 0), (3, 0), 'Helvetica-Bold'),
                                ('BOTTOMPADDING', (0, 0), (3, 0), 12),
                                ('BACKGROUND', (0, 1), (-1, -1), colors.lightgrey),
                                ('GRID', (0, 0), (-1, -1), 1, colors.black),
                                ('VALIGN', (0, 0), (-1, -1), 'MIDDLE')
                            ]))
                            elements.append(t)
                    
                    # Add experiments if available
                    if 'experiments' in item and item['experiments']:
                        elements.append(Spacer(1, 10))
                        elements.append(Paragraph("Experiments", section_style))
                        elements.append(Spacer(1, 6))
                        
                        exp_table = []
                        exp_table.append(["Property", "Value", "Conditions", "Date"])
                        
                        for exp in item['experiments']:
                            # Determine the value based on data type
                            if exp.get('numeric_value') is not None:
                                value = str(exp.get('numeric_value'))
                            elif exp.get('text_value') is not None:
                                value = exp.get('text_value')
                            elif exp.get('boolean_value') is not None:
                                value = str(exp.get('boolean_value'))
                            else:
                                value = "N/A"
                                
                            exp_table.append([
                                exp.get('property_types', {}).get('name', ''),
                                value,
                                exp.get('experimental_conditions', ''),
                                exp.get('date_performed', '')
                            ])
                        
                        if len(exp_table) > 1:
                            t = Table(exp_table, colWidths=[150, 100, 150, 100])
                            t.setStyle(TableStyle([
                                ('BACKGROUND', (0, 0), (3, 0), colors.darkblue),
                                ('TEXTCOLOR', (0, 0), (3, 0), colors.white),
                                ('ALIGN', (0, 0), (3, 0), 'CENTER'),
                                ('FONTNAME', (0, 0), (3, 0), 'Helvetica-Bold'),
                                ('BOTTOMPADDING', (0, 0), (3, 0), 12),
                                ('BACKGROUND', (0, 1), (-1, -1), colors.lightgrey),
                                ('GRID', (0, 0), (-1, -1), 1, colors.black),
                                ('VALIGN', (0, 0), (-1, -1), 'MIDDLE')
                            ]))
                            elements.append(t)
                
                elif data_type in ['predictions', 'experiments']:
                    # Create a table for predictions or experiments
                    if data_type == 'predictions':
                        table_data.append(["Property", "Value", "Confidence", "Method"])
                    else:
                        table_data.append(["Property", "Value", "Conditions", "Date"])
                    
                    # Determine the value based on data type
                    if item.get('numeric_value') is not None:
                        value = str(item.get('numeric_value'))
                    elif item.get('text_value') is not None:
                        value = item.get('text_value')
                    elif item.get('boolean_value') is not None:
                        value = str(item.get('boolean_value'))
                    else:
                        value = "N/A"
                    
                    if data_type == 'predictions':
                        table_data.append([
                            item.get('property_types', {}).get('name', ''),
                            value,
                            str(item.get('confidence', '')),
                            item.get('calculation_methods', {}).get('name', '')
                        ])
                    else:
                        table_data.append([
                            item.get('property_types', {}).get('name', ''),
                            value,
                            item.get('experimental_conditions', ''),
                            item.get('date_performed', '')
                        ])
                
                elif data_type == 'comparisons':
                    # Create a table for comparisons
                    table_data.append(["Property", "Predicted", "Experimental", "Difference", "% Error"])
                    
                    table_data.append([
                        item.get('property_name', ''),
                        str(item.get('prediction_value', '')),
                        str(item.get('experiment_value', '')),
                        str(item.get('difference', '')),
                        str(item.get('percent_error', ''))
                    ])
                
                elif data_type == 'protocols':
                    # Create a table for protocol steps
                    table_data.append(["Step", "Description", "Parameters"])
                    
                    if 'steps' in item:
                        for step_idx, step in enumerate(item['steps']):
                            step_desc = step.get('description', f"Step {step_idx+1}")
                            step_params = json.dumps(step.get('parameters', {}))
                            table_data.append([str(step_idx+1), step_desc, step_params])
                
                # Add the table to the document if it has data
                if len(table_data) > 1:
                    col_count = len(table_data[0])
                    col_width = (doc.width - 36) / col_count
                    t = Table(table_data, colWidths=[col_width] * col_count)
                    t.setStyle(TableStyle([
                        ('BACKGROUND', (0, 0), (-1, 0), colors.darkblue),
                        ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
                        ('ALIGN', (0, 0), (-1, 0), 'CENTER'),
                        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
                        ('BACKGROUND', (0, 1), (-1, -1), colors.lightgrey),
                        ('GRID', (0, 0), (-1, -1), 1, colors.black),
                        ('VALIGN', (0, 0), (-1, -1), 'MIDDLE')
                    ]))
                    elements.append(t)
        
        # Add footer with page numbers
        def add_page_number(canvas, doc):
            canvas.saveState()
            canvas.setFont('Helvetica', 9)
            page_num = canvas.getPageNumber()
            text = f"Page {page_num}"
            canvas.drawRightString(doc.pagesize[0] - 30, 30, text)
            canvas.drawString(30, 30, "CryoProtect Analyzer")
            canvas.restoreState()
        
        # Build the PDF
        doc.build(elements, onFirstPage=add_page_number, onLaterPages=add_page_number)
        buffer.seek(0)
        
        # Log export
        logger.info(f"Generated PDF export with {len(data)} records")
        
        return buffer
    except Exception as e:
        logger.error(f"Error generating PDF: {str(e)}")
        raise
    
    if not data:
        elements.append(Paragraph("No data available", styles["Normal"]))
    else:
        # For each item in data, create a section
        for i, item in enumerate(data):
            if i > 0:
                elements.append(Spacer(1, 20))
                
            # Item header
            item_name = item.get('name', f"Item {i+1}")
            elements.append(Paragraph(item_name, styles["Heading2"]))
            elements.append(Spacer(1, 6))
            
            # Convert item to table data
            table_data = []
            
            # Handle different data types
            if data_type == 'molecules':
                # Add basic molecule info
                table_data.append(["Property", "Value"])
                table_data.append(["ID", item.get('id', '')])
                table_data.append(["Name", item.get('name', '')])
                table_data.append(["Formula", item.get('molecular_formula', '')])
                table_data.append(["SMILES", item.get('smiles', '')])
                table_data.append(["PubChem CID", str(item.get('cid', ''))])
                
                # Add properties if available
                if 'properties' in item and item['properties']:
                    elements.append(Spacer(1, 10))
                    elements.append(Paragraph("Molecular Properties", styles["Heading3"]))
                    elements.append(Spacer(1, 6))
                    
                    prop_table = []
                    prop_table.append(["Property", "Value"])
                    
                    for prop_name, prop_value in item['properties'].items():
                        if isinstance(prop_value, (dict, list)):
                            prop_value = json.dumps(prop_value)
                        prop_table.append([prop_name, str(prop_value)])
                    
                    if len(prop_table) > 1:
                        t = Table(prop_table, colWidths=[200, 300])
                        t.setStyle(TableStyle([
                            ('BACKGROUND', (0, 0), (1, 0), colors.grey),
                            ('TEXTCOLOR', (0, 0), (1, 0), colors.whitesmoke),
                            ('ALIGN', (0, 0), (1, 0), 'CENTER'),
                            ('FONTNAME', (0, 0), (1, 0), 'Helvetica-Bold'),
                            ('BOTTOMPADDING', (0, 0), (1, 0), 12),
                            ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
                            ('GRID', (0, 0), (-1, -1), 1, colors.black)
                        ]))
                        elements.append(t)
            
            elif data_type == 'mixtures':
                # Add basic mixture info
                table_data.append(["Property", "Value"])
                table_data.append(["ID", item.get('id', '')])
                table_data.append(["Name", item.get('name', '')])
                table_data.append(["Description", item.get('description', '')])
                
                # Add components if available
                if 'components' in item and item['components']:
                    elements.append(Spacer(1, 10))
                    elements.append(Paragraph("Components", styles["Heading3"]))
                    elements.append(Spacer(1, 6))
                    
                    comp_table = []
                    comp_table.append(["Molecule", "Concentration", "Unit"])
                    
                    for comp in item['components']:
                        comp_table.append([
                            comp.get('molecule_name', ''),
                            str(comp.get('concentration', '')),
                            comp.get('concentration_unit', '')
                        ])
                    
                    if len(comp_table) > 1:
                        t = Table(comp_table, colWidths=[200, 150, 150])
                        t.setStyle(TableStyle([
                            ('BACKGROUND', (0, 0), (2, 0), colors.grey),
                            ('TEXTCOLOR', (0, 0), (2, 0), colors.whitesmoke),
                            ('ALIGN', (0, 0), (2, 0), 'CENTER'),
                            ('FONTNAME', (0, 0), (2, 0), 'Helvetica-Bold'),
                            ('BOTTOMPADDING', (0, 0), (2, 0), 12),
                            ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
                            ('GRID', (0, 0), (-1, -1), 1, colors.black)
                        ]))
                        elements.append(t)
                
                # Add predictions if available
                if 'predictions' in item and item['predictions']:
                    elements.append(Spacer(1, 10))
                    elements.append(Paragraph("Predictions", styles["Heading3"]))
                    elements.append(Spacer(1, 6))
                    
                    pred_table = []
                    pred_table.append(["Property", "Value", "Confidence", "Method"])
                    
                    for pred in item['predictions']:
                        # Determine the value based on data type
                        if pred.get('numeric_value') is not None:
                            value = str(pred.get('numeric_value'))
                        elif pred.get('text_value') is not None:
                            value = pred.get('text_value')
                        elif pred.get('boolean_value') is not None:
                            value = str(pred.get('boolean_value'))
                        else:
                            value = "N/A"
                            
                        pred_table.append([
                            pred.get('property_types', {}).get('name', ''),
                            value,
                            str(pred.get('confidence', '')),
                            pred.get('calculation_methods', {}).get('name', '')
                        ])
                    
                    if len(pred_table) > 1:
                        t = Table(pred_table, colWidths=[150, 150, 100, 100])
                        t.setStyle(TableStyle([
                            ('BACKGROUND', (0, 0), (3, 0), colors.grey),
                            ('TEXTCOLOR', (0, 0), (3, 0), colors.whitesmoke),
                            ('ALIGN', (0, 0), (3, 0), 'CENTER'),
                            ('FONTNAME', (0, 0), (3, 0), 'Helvetica-Bold'),
                            ('BOTTOMPADDING', (0, 0), (3, 0), 12),
                            ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
                            ('GRID', (0, 0), (-1, -1), 1, colors.black)
                        ]))
                        elements.append(t)
                
                # Add experiments if available
                if 'experiments' in item and item['experiments']:
                    elements.append(Spacer(1, 10))
                    elements.append(Paragraph("Experiments", styles["Heading3"]))
                    elements.append(Spacer(1, 6))
                    
                    exp_table = []
                    exp_table.append(["Property", "Value", "Conditions", "Date"])
                    
                    for exp in item['experiments']:
                        # Determine the value based on data type
                        if exp.get('numeric_value') is not None:
                            value = str(exp.get('numeric_value'))
                        elif exp.get('text_value') is not None:
                            value = exp.get('text_value')
                        elif exp.get('boolean_value') is not None:
                            value = str(exp.get('boolean_value'))
                        else:
                            value = "N/A"
                            
                        exp_table.append([
                            exp.get('property_types', {}).get('name', ''),
                            value,
                            exp.get('experimental_conditions', ''),
                            exp.get('date_performed', '')
                        ])
                    
                    if len(exp_table) > 1:
                        t = Table(exp_table, colWidths=[150, 100, 150, 100])
                        t.setStyle(TableStyle([
                            ('BACKGROUND', (0, 0), (3, 0), colors.grey),
                            ('TEXTCOLOR', (0, 0), (3, 0), colors.whitesmoke),
                            ('ALIGN', (0, 0), (3, 0), 'CENTER'),
                            ('FONTNAME', (0, 0), (3, 0), 'Helvetica-Bold'),
                            ('BOTTOMPADDING', (0, 0), (3, 0), 12),
                            ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
                            ('GRID', (0, 0), (-1, -1), 1, colors.black)
                        ]))
                        elements.append(t)
            
            elif data_type in ['predictions', 'experiments']:
                # Create a table for predictions or experiments
                if data_type == 'predictions':
                    table_data.append(["Property", "Value", "Confidence", "Method"])
                else:
                    table_data.append(["Property", "Value", "Conditions", "Date"])
                
                # Determine the value based on data type
                if item.get('numeric_value') is not None:
                    value = str(item.get('numeric_value'))
                elif item.get('text_value') is not None:
                    value = item.get('text_value')
                elif item.get('boolean_value') is not None:
                    value = str(item.get('boolean_value'))
                else:
                    value = "N/A"
                
                if data_type == 'predictions':
                    table_data.append([
                        item.get('property_types', {}).get('name', ''),
                        value,
                        str(item.get('confidence', '')),
                        item.get('calculation_methods', {}).get('name', '')
                    ])
                else:
                    table_data.append([
                        item.get('property_types', {}).get('name', ''),
                        value,
                        item.get('experimental_conditions', ''),
                        item.get('date_performed', '')
                    ])
            
            elif data_type == 'comparisons':
                # Create a table for comparisons
                table_data.append(["Property", "Predicted", "Experimental", "Difference", "% Error"])
                
                table_data.append([
                    item.get('property_name', ''),
                    str(item.get('prediction_value', '')),
                    str(item.get('experiment_value', '')),
                    str(item.get('difference', '')),
                    str(item.get('percent_error', ''))
                ])
            
            # Add the table to the document if it has data
            if len(table_data) > 1:
                t = Table(table_data, colWidths=[120] * (len(table_data[0])))
                t.setStyle(TableStyle([
                    ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
                    ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                    ('ALIGN', (0, 0), (-1, 0), 'CENTER'),
                    ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                    ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
                    ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
                    ('GRID', (0, 0), (-1, -1), 1, colors.black)
                ]))
                elements.append(t)
    
    # Build the PDF
    doc.build(elements)
    buffer.seek(0)
    return buffer


# Helper functions for visualization export
def generate_visualization(
    chart_type: str,
    data: Dict[str, Any],
    format: str,
    width: int,
    height: int,
    title: str,
    style: str,
    dpi: int = 100,
    transparent: bool = False,
    include_legend: bool = True,
    color_scheme: str = 'default'
) -> io.BytesIO:
    """
    Generate visualization based on chart type and data with enhanced formatting and options.
    
    Args:
        chart_type: Type of chart to generate
        data: Data for the visualization
        format: Output format (png, svg, pdf)
        width: Width of the visualization in pixels
        height: Height of the visualization in pixels
        title: Title of the visualization
        style: Style of the visualization
        dpi: Resolution for raster formats
        transparent: Whether to use transparent background
        include_legend: Whether to include a legend
        color_scheme: Color scheme to use
        
    Returns:
        BytesIO object containing the visualization
    """
    try:
        # Set the style
        if style == 'dark':
            plt.style.use('dark_background')
        elif style == 'seaborn':
            plt.style.use('seaborn')
        elif style == 'ggplot':
            plt.style.use('ggplot')
        elif style == 'bmh':
            plt.style.use('bmh')
        elif style == 'fivethirtyeight':
            plt.style.use('fivethirtyeight')
        else:
            plt.style.use('default')
        
        # Set color scheme
        if color_scheme == 'viridis':
            colors = plt.cm.viridis.colors
        elif color_scheme == 'plasma':
            colors = plt.cm.plasma.colors
        elif color_scheme == 'inferno':
            colors = plt.cm.inferno.colors
        elif color_scheme == 'magma':
            colors = plt.cm.magma.colors
        elif color_scheme == 'cividis':
            colors = plt.cm.cividis.colors
        elif color_scheme == 'tab10':
            colors = plt.cm.tab10.colors
        elif color_scheme == 'tab20':
            colors = plt.cm.tab20.colors
        else:
            colors = None
        
        # Create figure with specified dimensions
        fig, ax = plt.subplots(figsize=(width/100, height/100), dpi=dpi)
        
        if chart_type == 'property_comparison':
            # Property comparison chart (bar chart)
            properties = data.get('properties', [])
            values = data.get('values', [])
            errors = data.get('errors', None)
            
            if not properties or not values or len(properties) != len(values):
                abort(400, message="Invalid data for property comparison chart")
            
            # Create bar chart with optional error bars
            if errors and len(errors) == len(values):
                bars = ax.bar(properties, values, yerr=errors, capsize=5,
                             color=colors[0:len(properties)] if colors else None)
            else:
                bars = ax.bar(properties, values,
                             color=colors[0:len(properties)] if colors else None)
            
            # Add value labels on top of bars
            if data.get('show_values', True):
                for bar in bars:
                    height = bar.get_height()
                    ax.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                           f'{height:.2f}', ha='center', va='bottom', fontsize=9)
            
            ax.set_xlabel('Properties', fontweight='bold')
            ax.set_ylabel('Values', fontweight='bold')
            
            # Rotate x-axis labels if there are many properties
            if len(properties) > 5:
                plt.xticks(rotation=45, ha='right')
            
            # Add grid lines
            ax.grid(axis='y', linestyle='--', alpha=0.7)
            
        elif chart_type == 'mixture_composition':
            # Mixture composition chart (pie chart)
            components = data.get('components', [])
            
            if not components:
                abort(400, message="Invalid data for mixture composition chart")
            
            labels = [comp.get('name', f"Component {i+1}") for i, comp in enumerate(components)]
            values = [comp.get('concentration', 0) for comp in components]
            
            # Optional explode to emphasize specific components
            explode = data.get('explode', None)
            if not explode:
                explode = [0.1 if i == data.get('highlight_index', -1) else 0 for i in range(len(components))]
            
            # Create pie chart
            wedges, texts, autotexts = ax.pie(
                values,
                labels=None if data.get('legend_only', False) else labels,
                autopct='%1.1f%%' if data.get('show_percentages', True) else None,
                startangle=data.get('start_angle', 90),
                explode=explode,
                shadow=data.get('shadow', False),
                colors=colors[0:len(components)] if colors else None
            )
            
            # Customize text
            for autotext in autotexts:
                autotext.set_fontsize(9)
                autotext.set_fontweight('bold')
            
            ax.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle
            
            # Add legend if requested
            if include_legend or data.get('legend_only', False):
                ax.legend(wedges, labels, title="Components",
                         loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))
            
        elif chart_type == 'prediction_accuracy':
            # Prediction accuracy chart (scatter plot with line)
            predictions = data.get('predictions', [])
            experiments = data.get('experiments', [])
            labels = data.get('labels', [])
            
            if not predictions or not experiments:
                abort(400, message="Invalid data for prediction accuracy chart")
            
            # Create a perfect prediction line
            min_val = min(min(predictions), min(experiments))
            max_val = max(max(predictions), max(experiments))
            buffer_val = (max_val - min_val) * 0.1  # Add 10% buffer
            line_vals = [min_val - buffer_val, max_val + buffer_val]
            
            # Plot perfect prediction line
            ax.plot(line_vals, line_vals, 'k--', alpha=0.7, linewidth=2, label='Perfect Prediction')
            
            # Plot actual predictions vs experiments
            if labels and len(labels) == len(predictions):
                # Use labels for scatter points
                for i, (pred, exp, label) in enumerate(zip(predictions, experiments, labels)):
                    ax.scatter(pred, exp, alpha=0.8, label=label,
                              color=colors[i % len(colors)] if colors else None)
                
                if include_legend:
                    ax.legend(title="Data Points", loc="best")
            else:
                # Simple scatter plot
                sc = ax.scatter(predictions, experiments, alpha=0.8, s=80,
                               color=colors[0] if colors else None)
                
                # Add error bands if provided
                if data.get('error_bands', False):
                    error_margin = data.get('error_margin', 0.1)  # 10% error by default
                    upper_line = [(x * (1 + error_margin)) for x in line_vals]
                    lower_line = [(x * (1 - error_margin)) for x in line_vals]
                    ax.plot(line_vals, upper_line, 'r--', alpha=0.3)
                    ax.plot(line_vals, lower_line, 'r--', alpha=0.3)
                    ax.fill_between(line_vals, lower_line, upper_line, alpha=0.1, color='red')
            
            ax.set_xlabel('Predicted Values', fontweight='bold')
            ax.set_ylabel('Experimental Values', fontweight='bold')
            
            # Add grid
            ax.grid(True, linestyle='--', alpha=0.6)
            
            # Add R² value if requested
            if data.get('show_r2', True):
                from scipy import stats
                slope, intercept, r_value, p_value, std_err = stats.linregress(predictions, experiments)
                r_squared = r_value ** 2
                ax.text(0.05, 0.95, f'R² = {r_squared:.4f}', transform=ax.transAxes,
                       fontsize=10, fontweight='bold', bbox=dict(facecolor='white', alpha=0.7))
            
        elif chart_type == 'time_series':
            # Time series chart
            times = data.get('times', [])
            series = data.get('series', [])
            series_names = data.get('series_names', [])
            
            if not times or not series or not all(len(s) == len(times) for s in series):
                abort(400, message="Invalid data for time series chart")
            
            # Plot each series
            for i, s in enumerate(series):
                name = series_names[i] if i < len(series_names) else f"Series {i+1}"
                ax.plot(times, s, label=name, linewidth=2, marker=data.get('marker', 'o'),
                       markersize=data.get('marker_size', 5),
                       color=colors[i % len(colors)] if colors else None)
            
            ax.set_xlabel('Time', fontweight='bold')
            ax.set_ylabel('Value', fontweight='bold')
            
            # Format x-axis as dates if specified
            if data.get('date_format'):
                from matplotlib import dates as mdates
                ax.xaxis.set_major_formatter(mdates.DateFormatter(data.get('date_format')))
                fig.autofmt_xdate()
            
            # Add grid
            ax.grid(True, linestyle='--', alpha=0.6)
            
            if include_legend and len(series) > 1:
                ax.legend(title=data.get('legend_title', 'Series'), loc='best')
            
        elif chart_type == 'heatmap':
            # Heatmap chart
            matrix = data.get('matrix', [])
            x_labels = data.get('x_labels', [])
            y_labels = data.get('y_labels', [])
            
            if not matrix:
                abort(400, message="Invalid data for heatmap chart")
            
            # Create heatmap
            cmap = data.get('colormap', 'viridis')
            im = ax.imshow(matrix, cmap=cmap)
            
            # Add colorbar
            cbar = fig.colorbar(im, ax=ax)
            cbar.set_label(data.get('colorbar_label', 'Value'))
            
            # Set ticks and labels
            if x_labels:
                ax.set_xticks(range(len(x_labels)))
                ax.set_xticklabels(x_labels)
                plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
            
            if y_labels:
                ax.set_yticks(range(len(y_labels)))
                ax.set_yticklabels(y_labels)
            
            # Add text annotations if requested
            if data.get('show_values', True):
                for i in range(len(matrix)):
                    for j in range(len(matrix[i])):
                        text_color = 'white' if im.norm(matrix[i][j]) > 0.5 else 'black'
                        ax.text(j, i, f"{matrix[i][j]:.2f}", ha="center", va="center", color=text_color)
            
            ax.set_xlabel(data.get('x_label', 'X Axis'), fontweight='bold')
            ax.set_ylabel(data.get('y_label', 'Y Axis'), fontweight='bold')
            
        elif chart_type == 'radar':
            # Radar chart (polar plot)
            categories = data.get('categories', [])
            values = data.get('values', [])
            
            if not categories or not values or len(categories) != len(values):
                abort(400, message="Invalid data for radar chart")
            
            # Close the loop for the radar chart
            categories = categories + [categories[0]]
            values = values + [values[0]]
            
            # Calculate angles for each category
            angles = np.linspace(0, 2*np.pi, len(categories), endpoint=True)
            
            # Plot radar chart
            ax = plt.subplot(111, polar=True)
            ax.plot(angles, values, 'o-', linewidth=2,
                   color=colors[0] if colors else None)
            ax.fill(angles, values, alpha=0.25,
                   color=colors[0] if colors else None)
            
            # Set category labels
            ax.set_thetagrids(angles[:-1] * 180/np.pi, categories[:-1])
            
            # Add value labels if requested
            if data.get('show_values', True):
                for angle, value, category in zip(angles, values, categories):
                    ax.text(angle, value*1.1, f"{value:.2f}",
                           ha='center', va='center', fontweight='bold')
            
            # Customize grid
            ax.set_rlabel_position(0)
            if data.get('grid_step'):
                plt.yticks(
                    np.arange(0, max(values)*1.1, data.get('grid_step')),
                    color="grey", size=7
                )
            
            # Add second series if provided
            if 'values2' in data:
                values2 = data['values2'] + [data['values2'][0]]  # Close the loop
                ax.plot(angles, values2, 'o-', linewidth=2,
                       color=colors[1] if colors and len(colors) > 1 else None)
                ax.fill(angles, values2, alpha=0.25,
                       color=colors[1] if colors and len(colors) > 1 else None)
                
                if include_legend:
                    ax.legend([data.get('label1', 'Series 1'), data.get('label2', 'Series 2')])
            
        elif chart_type == 'bubble':
            # Bubble chart (scatter plot with varying bubble sizes)
            x_values = data.get('x_values', [])
            y_values = data.get('y_values', [])
            sizes = data.get('sizes', [])
            labels = data.get('labels', [])
            
            if not x_values or not y_values or not sizes or len(x_values) != len(y_values) or len(x_values) != len(sizes):
                abort(400, message="Invalid data for bubble chart")
            
            # Scale sizes appropriately
            size_scale = data.get('size_scale', 1000)
            scaled_sizes = [s * size_scale for s in sizes]
            
            # Create bubble chart
            scatter = ax.scatter(
                x_values, y_values, s=scaled_sizes, alpha=0.6,
                c=range(len(x_values)) if colors is None else None,
                cmap=data.get('colormap', 'viridis') if colors is None else None
            )
            
            # Add labels if provided
            if labels and len(labels) == len(x_values) and data.get('show_labels', False):
                for i, (x, y, label) in enumerate(zip(x_values, y_values, labels)):
                    ax.annotate(label, (x, y), xytext=(5, 5), textcoords='offset points')
            
            ax.set_xlabel(data.get('x_label', 'X Axis'), fontweight='bold')
            ax.set_ylabel(data.get('y_label', 'Y Axis'), fontweight='bold')
            
            # Add grid
            ax.grid(True, linestyle='--', alpha=0.6)
            
            # Add size legend if requested
            if data.get('size_legend', False):
                # Create size legend
                size_legend_values = data.get('size_legend_values', [min(sizes), (min(sizes)+max(sizes))/2, max(sizes)])
                size_legend_labels = [f"{v:.2f}" for v in size_legend_values]
                
                # Create proxy artists for the legend
                legend_sizes = [s * size_scale for s in size_legend_values]
                legend_handles = [plt.scatter([], [], s=s, color='gray', alpha=0.6) for s in legend_sizes]
                
                # Add the legend
                ax.legend(legend_handles, size_legend_labels,
                         title=data.get('size_legend_title', 'Size'),
                         loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0)
            
        elif chart_type == 'custom':
            # Custom chart based on provided data
            chart_subtype = data.get('chart_subtype', 'bar')
            x_data = data.get('x_data', [])
            y_data = data.get('y_data', [])
            
            if not x_data or not y_data or len(x_data) != len(y_data):
                abort(400, message="Invalid data for custom chart")
            
            if chart_subtype == 'bar':
                bars = ax.bar(x_data, y_data,
                             color=colors[0:len(x_data)] if colors else None)
                
                # Add value labels on top of bars if requested
                if data.get('show_values', True):
                    for bar in bars:
                        height = bar.get_height()
                        ax.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                               f'{height:.2f}', ha='center', va='bottom', fontsize=9)
                
                # Rotate x-axis labels if there are many categories
                if len(x_data) > 5:
                    plt.xticks(rotation=45, ha='right')
                
                # Add grid lines
                ax.grid(axis='y', linestyle='--', alpha=0.7)
                
            elif chart_subtype == 'line':
                ax.plot(x_data, y_data, marker=data.get('marker', 'o'),
                       linewidth=data.get('line_width', 2), markersize=data.get('marker_size', 5),
                       color=colors[0] if colors else None)
                
                # Add grid
                ax.grid(True, linestyle='--', alpha=0.6)
                
            elif chart_subtype == 'scatter':
                ax.scatter(x_data, y_data, s=data.get('point_size', 50), alpha=0.7,
                          c=data.get('colors', range(len(x_data))),
                          cmap=data.get('colormap', 'viridis'))
                
                # Add trend line if requested
                if data.get('trend_line', False):
                    z = np.polyfit(x_data, y_data, 1)
                    p = np.poly1d(z)
                    ax.plot(x_data, p(x_data), "r--", alpha=0.8)
                
                # Add grid
                ax.grid(True, linestyle='--', alpha=0.6)
                
            elif chart_subtype == 'pie':
                wedges, texts, autotexts = ax.pie(
                    y_data,
                    labels=None if data.get('legend_only', False) else x_data,
                    autopct='%1.1f%%' if data.get('show_percentages', True) else None,
                    startangle=data.get('start_angle', 90),
                    shadow=data.get('shadow', False),
                    colors=colors[0:len(y_data)] if colors else None
                )
                
                # Customize text
                for autotext in autotexts:
                    autotext.set_fontsize(9)
                    autotext.set_fontweight('bold')
                
                ax.axis('equal')
                
                # Add legend if requested
                if include_legend or data.get('legend_only', False):
                    ax.legend(wedges, x_data, title=data.get('legend_title', 'Categories'),
                             loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))
                
            elif chart_subtype == 'area':
                ax.fill_between(x_data, y_data, alpha=0.5,
                               color=colors[0] if colors else None)
                ax.plot(x_data, y_data, 'o-', linewidth=2,
                       color=colors[0] if colors else None)
                
                # Add grid
                ax.grid(True, linestyle='--', alpha=0.6)
                
            elif chart_subtype == 'horizontal_bar':
                bars = ax.barh(x_data, y_data,
                              color=colors[0:len(x_data)] if colors else None)
                
                # Add value labels
                if data.get('show_values', True):
                    for bar in bars:
                        width = bar.get_width()
                        ax.text(width + 0.01, bar.get_y() + bar.get_height()/2.,
                               f'{width:.2f}', va='center', fontsize=9)
                
                # Add grid lines
                ax.grid(axis='x', linestyle='--', alpha=0.7)
                
            else:
                abort(400, message=f"Unsupported chart subtype: {chart_subtype}")
            
            ax.set_xlabel(data.get('x_label', 'X Axis'), fontweight='bold')
            ax.set_ylabel(data.get('y_label', 'Y Axis'), fontweight='bold')
        
        else:
            abort(400, message=f"Unsupported chart type: {chart_type}")
        
        # Set title with custom formatting
        if title:
            ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
        
        # Add a text watermark if specified
        if data.get('watermark'):
            fig.text(0.5, 0.5, data.get('watermark'), fontsize=40, color='gray',
                    ha='center', va='center', alpha=0.2, rotation=30)
        
        # Adjust layout
        plt.tight_layout()
        
        # Save to buffer in requested format
        buffer = io.BytesIO()
        
        if format == 'svg':
            plt.savefig(buffer, format='svg', bbox_inches='tight', transparent=transparent)
        elif format == 'pdf':
            plt.savefig(buffer, format='pdf', bbox_inches='tight', transparent=transparent)
        else:  # Default to PNG
            plt.savefig(buffer, format='png', bbox_inches='tight', transparent=transparent, dpi=dpi)
        
        plt.close(fig)
        buffer.seek(0)
        
        # Log export
        logger.info(f"Generated {format} visualization of type {chart_type}")
        
        return buffer
    except Exception as e:
        logger.error(f"Error generating visualization: {str(e)}")
        raise


# Helper functions for report generation
def generate_report(title: str, sections: List[Dict[str, Any]], include_visualizations: bool, include_data_tables: bool, template_id: Optional[str] = None) -> io.BytesIO:
    """
    Generate a comprehensive report with sections, visualizations, and data tables.
    
    Args:
        title: Report title
        sections: List of section dictionaries with content
        include_visualizations: Whether to include visualizations
        include_data_tables: Whether to include data tables
        template_id: Optional template ID for styling
        
    Returns:
        BytesIO object containing the report PDF
    """
    buffer = io.BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=letter)
    styles = getSampleStyleSheet()
    elements = []
    
    # Add title
    title_style = styles["Heading1"]
    elements.append(Paragraph(title, title_style))
    elements.append(Spacer(1, 12))
    
    # Add timestamp
    date_style = styles["Normal"]
    date_text = f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
    elements.append(Paragraph(date_text, date_style))
    elements.append(Spacer(1, 20))
    
    # Process each section
    for section in sections:
        section_title = section.get('title', 'Untitled Section')
        section_content = section.get('content', '')
        section_data = section.get('data', [])
        section_visualization = section.get('visualization', {})
        
        # Add section title
        section_title_style = styles["Heading2"]
        elements.append(Paragraph(section_title, section_title_style))
        elements.append(Spacer(1, 6))
        
        # Add section content
        if section_content:
            content_style = styles["Normal"]
            elements.append(Paragraph(section_content, content_style))
            elements.append(Spacer(1, 10))
        
        # Add visualization if requested and available
        if include_visualizations and section_visualization:
            try:
                chart_type = section_visualization.get('chart_type', 'bar')
                chart_data = section_visualization.get('data', {})
                chart_title = section_visualization.get('title', '')
                chart_format = section_visualization.get('format', 'png')
                chart_width = section_visualization.get('width', 800)
                chart_height = section_visualization.get('height', 600)
                chart_style = section_visualization.get('style', 'default')
                
                # Generate the visualization using the existing function
                viz_buffer = generate_visualization(
                    chart_type=chart_type,
                    data=chart_data,
                    format=chart_format,
                    width=chart_width,
                    height=chart_height,
                    title=chart_title,
                    style=chart_style
                )
                
                # Create a temporary file to store the image
                with tempfile.NamedTemporaryFile(delete=False, suffix=f'.{chart_format}') as tmp:
                    tmp.write(viz_buffer.getvalue())
                    tmp_path = tmp.name
                
                # Add the image to the report
                img = Image(tmp_path, width=chart_width/2, height=chart_height/2)
                elements.append(img)
                elements.append(Spacer(1, 10))
                
                # Store the generated visualization path
                section_visualization['generated'] = tmp_path
            except Exception as e:
                section_visualization['error'] = str(e)
                # Add error message to the report
                error_style = ParagraphStyle('ErrorStyle', parent=styles['Normal'], textColor=colors.red)
                elements.append(Paragraph(f"Error generating visualization: {str(e)}", error_style))
                elements.append(Spacer(1, 10))
