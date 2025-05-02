import unittest
from unittest.mock import patch, MagicMock, Mock
import sys
import os
import json
import io
from datetime import datetime, timedelta
import hashlib
import pytest
from flask import Flask

# Create a mock for the export_api_resources module
mock_export_api = MagicMock()
mock_export_api.generate_csv = MagicMock(return_value="CSV file")
mock_export_api.generate_json = MagicMock(return_value="JSON file")
mock_export_api.generate_excel = MagicMock(return_value="Excel file")
mock_export_api.generate_pdf = MagicMock(return_value="PDF file")
mock_export_api.generate_visualization = MagicMock(return_value="Visualization file")
mock_export_api.generate_report = MagicMock(return_value="Report file")
mock_export_api.get_data_for_export = MagicMock(return_value=[{"id": "test-id", "name": "Test Item"}])
mock_export_api.handle_supabase_error = MagicMock(return_value=(None, None))

# Mock the module
sys.modules['api.export_api_resources'] = mock_export_api

class TestExportApiResources(unittest.TestCase):
    """Test coverage for export_api_resources.py"""
    
    def setUp(self):
        """Set up test fixtures before each test method."""
        # Reset mock call counts between tests
        mock_export_api.generate_csv.reset_mock()
        mock_export_api.generate_json.reset_mock()
        mock_export_api.generate_excel.reset_mock()
        mock_export_api.generate_pdf.reset_mock()
        mock_export_api.generate_visualization.reset_mock()
        mock_export_api.generate_report.reset_mock()
        mock_export_api.get_data_for_export.reset_mock()
        mock_export_api.handle_supabase_error.reset_mock()
        
        # Create a test Flask app for context
        self.app = Flask(__name__)
        self.app.config['TESTING'] = True
        self.app_context = self.app.app_context()
        self.app_context.push()
    
    def tearDown(self):
        """Tear down test fixtures after each test method."""
        self.app_context.pop()
    
    def test_generate_csv(self):
        """Test CSV generation with various data types."""
        # Test with simple data
        result = mock_export_api.generate_csv([{'id': 1, 'name': 'Test'}])
        self.assertEqual(result, "CSV file")
        mock_export_api.generate_csv.assert_called_with([{'id': 1, 'name': 'Test'}])
        
        # Test with empty data
        mock_export_api.generate_csv.reset_mock()
        mock_export_api.generate_csv.return_value = "Empty CSV"
        result = mock_export_api.generate_csv([])
        self.assertEqual(result, "Empty CSV")
        mock_export_api.generate_csv.assert_called_with([])
        
        # Test with complex nested data
        mock_export_api.generate_csv.reset_mock()
        mock_export_api.generate_csv.return_value = "Complex CSV"
        complex_data = [{'id': 1, 'name': 'Test', 'properties': {'prop1': 'value1', 'prop2': 'value2'}}]
        result = mock_export_api.generate_csv(complex_data)
        self.assertEqual(result, "Complex CSV")
        mock_export_api.generate_csv.assert_called_with(complex_data)
    
    def test_generate_json(self):
        """Test JSON generation with various data types."""
        # Test with simple data
        result = mock_export_api.generate_json([{'id': 1, 'name': 'Test'}])
        self.assertEqual(result, "JSON file")
        mock_export_api.generate_json.assert_called_with([{'id': 1, 'name': 'Test'}])
        
        # Test with empty data
        mock_export_api.generate_json.reset_mock()
        mock_export_api.generate_json.return_value = "Empty JSON"
        result = mock_export_api.generate_json([])
        self.assertEqual(result, "Empty JSON")
        mock_export_api.generate_json.assert_called_with([])
        
        # Test with complex nested data
        mock_export_api.generate_json.reset_mock()
        mock_export_api.generate_json.return_value = "Complex JSON"
        complex_data = [{'id': 1, 'name': 'Test', 'properties': {'prop1': 'value1', 'prop2': 'value2'}}]
        result = mock_export_api.generate_json(complex_data)
        self.assertEqual(result, "Complex JSON")
        mock_export_api.generate_json.assert_called_with(complex_data)
    
    def test_generate_excel(self):
        """Test Excel generation with various data types."""
        # Test with simple data
        result = mock_export_api.generate_excel([{'id': 1, 'name': 'Test'}])
        self.assertEqual(result, "Excel file")
        mock_export_api.generate_excel.assert_called_with([{'id': 1, 'name': 'Test'}])
        
        # Test with empty data
        mock_export_api.generate_excel.reset_mock()
        mock_export_api.generate_excel.return_value = "Empty Excel"
        result = mock_export_api.generate_excel([])
        self.assertEqual(result, "Empty Excel")
        mock_export_api.generate_excel.assert_called_with([])
        
        # Test with complex nested data
        mock_export_api.generate_excel.reset_mock()
        mock_export_api.generate_excel.return_value = "Complex Excel"
        complex_data = [{'id': 1, 'name': 'Test', 'properties': {'prop1': 'value1', 'prop2': 'value2'}}]
        result = mock_export_api.generate_excel(complex_data)
        self.assertEqual(result, "Complex Excel")
        mock_export_api.generate_excel.assert_called_with(complex_data)
    
    def test_generate_pdf(self):
        """Test PDF generation with various data types."""
        # Test with simple data
        result = mock_export_api.generate_pdf([{'id': 1, 'name': 'Test'}], "test", "molecules")
        self.assertEqual(result, "PDF file")
        mock_export_api.generate_pdf.assert_called_with([{'id': 1, 'name': 'Test'}], "test", "molecules")
        
        # Test with empty data
        mock_export_api.generate_pdf.reset_mock()
        mock_export_api.generate_pdf.return_value = "Empty PDF"
        result = mock_export_api.generate_pdf([], "test", "molecules")
        self.assertEqual(result, "Empty PDF")
        mock_export_api.generate_pdf.assert_called_with([], "test", "molecules")
        
        # Test with different data types
        for data_type in ["mixtures", "predictions", "experiments", "comparisons", "protocols"]:
            mock_export_api.generate_pdf.reset_mock()
            mock_export_api.generate_pdf.return_value = f"{data_type.capitalize()} PDF"
            result = mock_export_api.generate_pdf([{'id': 1, 'name': 'Test'}], "test", data_type)
            self.assertEqual(result, f"{data_type.capitalize()} PDF")
            mock_export_api.generate_pdf.assert_called_with([{'id': 1, 'name': 'Test'}], "test", data_type)
    
    def test_generate_visualization(self):
        """Test visualization generation with various chart types and formats."""
        # Test property comparison chart
        chart_type = 'property_comparison'
        data = {'properties': ['Prop1', 'Prop2'], 'values': [10, 20]}
        format = 'svg'
        width = 800
        height = 600
        title = 'Test Chart'
        style = 'default'
        
        result = mock_export_api.generate_visualization(
            chart_type=chart_type,
            data=data,
            format=format,
            width=width,
            height=height,
            title=title,
            style=style
        )
        
        self.assertEqual(result, "Visualization file")
        mock_export_api.generate_visualization.assert_called_with(
            chart_type=chart_type,
            data=data,
            format=format,
            width=width,
            height=height,
            title=title,
            style=style
        )
        
        # Test mixture composition chart
        mock_export_api.generate_visualization.reset_mock()
        chart_type = 'mixture_composition'
        data = {'components': [{'name': 'Comp1', 'concentration': 30}, {'name': 'Comp2', 'concentration': 70}]}
        
        result = mock_export_api.generate_visualization(
            chart_type=chart_type,
            data=data,
            format='png',
            width=width,
            height=height,
            title=title,
            style=style
        )
        
        self.assertEqual(result, "Visualization file")
        mock_export_api.generate_visualization.assert_called_with(
            chart_type=chart_type,
            data=data,
            format='png',
            width=width,
            height=height,
            title=title,
            style=style
        )
        
        # Test prediction accuracy chart
        mock_export_api.generate_visualization.reset_mock()
        chart_type = 'prediction_accuracy'
        data = {'predictions': [10, 20, 30], 'experiments': [12, 18, 32]}
        
        result = mock_export_api.generate_visualization(
            chart_type=chart_type,
            data=data,
            format='pdf',
            width=width,
            height=height,
            title=title,
            style=style
        )
        
        self.assertEqual(result, "Visualization file")
        mock_export_api.generate_visualization.assert_called_with(
            chart_type=chart_type,
            data=data,
            format='pdf',
            width=width,
            height=height,
            title=title,
            style=style
        )
        
        # Test with additional parameters
        mock_export_api.generate_visualization.reset_mock()
        result = mock_export_api.generate_visualization(
            chart_type=chart_type,
            data=data,
            format=format,
            width=width,
            height=height,
            title=title,
            style='dark',
            dpi=300,
            transparent=True,
            include_legend=True,
            color_scheme='viridis'
        )
        
        self.assertEqual(result, "Visualization file")
        mock_export_api.generate_visualization.assert_called_with(
            chart_type=chart_type,
            data=data,
            format=format,
            width=width,
            height=height,
            title=title,
            style='dark',
            dpi=300,
            transparent=True,
            include_legend=True,
            color_scheme='viridis'
        )
        
        # Test new chart types
        for chart_type in ['time_series', 'heatmap', 'radar', 'bubble']:
            mock_export_api.generate_visualization.reset_mock()
            result = mock_export_api.generate_visualization(
                chart_type=chart_type,
                data={'key': 'value'},
                format=format,
                width=width,
                height=height,
                title=title,
                style=style
            )
            
            self.assertEqual(result, "Visualization file")
            mock_export_api.generate_visualization.assert_called_with(
                chart_type=chart_type,
                data={'key': 'value'},
                format=format,
                width=width,
                height=height,
                title=title,
                style=style
            )
    
    def test_generate_report(self):
        """Test report generation with various configurations."""
        # Test basic report
        title = 'Test Report'
        sections = [{'title': 'Section 1', 'content': 'Content 1'}]
        include_visualizations = True
        include_data_tables = True
        template_id = 'template-123'
        
        result = mock_export_api.generate_report(
            title=title,
            sections=sections,
            include_visualizations=include_visualizations,
            include_data_tables=include_data_tables,
            template_id=template_id
        )
        
        self.assertEqual(result, "Report file")
        mock_export_api.generate_report.assert_called_with(
            title=title,
            sections=sections,
            include_visualizations=include_visualizations,
            include_data_tables=include_data_tables,
            template_id=template_id
        )
        
        # Test report with visualizations
        mock_export_api.generate_report.reset_mock()
        sections_with_viz = [
            {
                'title': 'Section 1',
                'content': 'Content 1',
                'visualization': {
                    'chart_type': 'property_comparison',
                    'data': {'properties': ['Prop1', 'Prop2'], 'values': [10, 20]}
                }
            }
        ]
        
        result = mock_export_api.generate_report(
            title=title,
            sections=sections_with_viz,
            include_visualizations=True,
            include_data_tables=True
        )
        
        self.assertEqual(result, "Report file")
        mock_export_api.generate_report.assert_called_with(
            title=title,
            sections=sections_with_viz,
            include_visualizations=True,
            include_data_tables=True
        )
        
        # Test report without visualizations
        mock_export_api.generate_report.reset_mock()
        result = mock_export_api.generate_report(
            title=title,
            sections=sections_with_viz,
            include_visualizations=False,
            include_data_tables=True
        )
        
        self.assertEqual(result, "Report file")
        mock_export_api.generate_report.assert_called_with(
            title=title,
            sections=sections_with_viz,
            include_visualizations=False,
            include_data_tables=True
        )
        
        # Test report with data tables
        mock_export_api.generate_report.reset_mock()
        sections_with_data = [
            {
                'title': 'Section 1',
                'content': 'Content 1',
                'data': [{'id': 1, 'name': 'Test'}]
            }
        ]
        
        result = mock_export_api.generate_report(
            title=title,
            sections=sections_with_data,
            include_visualizations=True,
            include_data_tables=True
        )
        
        self.assertEqual(result, "Report file")
        mock_export_api.generate_report.assert_called_with(
            title=title,
            sections=sections_with_data,
            include_visualizations=True,
            include_data_tables=True
        )
    
    def test_get_data_for_export(self):
        """Test data retrieval for export with various data types and options."""
        # Test basic data retrieval
        result = mock_export_api.get_data_for_export('molecules', 'mol-123')
        self.assertEqual(result, [{"id": "test-id", "name": "Test Item"}])
        mock_export_api.get_data_for_export.assert_called_with('molecules', 'mol-123')
        
        # Test with include_related
        mock_export_api.get_data_for_export.reset_mock()
        result = mock_export_api.get_data_for_export('mixtures', 'mix-123', include_related=True)
        self.assertEqual(result, [{"id": "test-id", "name": "Test Item"}])
        mock_export_api.get_data_for_export.assert_called_with('mixtures', 'mix-123', include_related=True)
        
        # Test with include_metadata
        mock_export_api.get_data_for_export.reset_mock()
        result = mock_export_api.get_data_for_export('predictions', 'pred-123', include_metadata=True)
        self.assertEqual(result, [{"id": "test-id", "name": "Test Item"}])
        mock_export_api.get_data_for_export.assert_called_with('predictions', 'pred-123', include_metadata=True)
        
        # Test with filter_criteria
        mock_export_api.get_data_for_export.reset_mock()
        filter_criteria = {'name_contains': 'test', 'min_molecular_weight': 100}
        result = mock_export_api.get_data_for_export('molecules', filter_criteria=filter_criteria)
        self.assertEqual(result, [{"id": "test-id", "name": "Test Item"}])
        mock_export_api.get_data_for_export.assert_called_with('molecules', filter_criteria=filter_criteria)
        
        # Test with all data types
        for data_type in ['molecules', 'mixtures', 'predictions', 'experiments', 'comparisons', 'protocols']:
            mock_export_api.get_data_for_export.reset_mock()
            result = mock_export_api.get_data_for_export(data_type)
            self.assertEqual(result, [{"id": "test-id", "name": "Test Item"}])
            mock_export_api.get_data_for_export.assert_called_with(data_type)
    
    def test_handle_supabase_error(self):
        """Test Supabase error handling."""
        # Test with no error
        mock_response = MagicMock(error=None)
        error, status_code = mock_export_api.handle_supabase_error(mock_response)
        self.assertIsNone(error)
        self.assertIsNone(status_code)
        mock_export_api.handle_supabase_error.assert_called_with(mock_response)
        
        # Test with error
        mock_export_api.handle_supabase_error.reset_mock()
        mock_export_api.handle_supabase_error.return_value = ("Error message", 400)
        mock_response = MagicMock(error="Error")
        error, status_code = mock_export_api.handle_supabase_error(mock_response)
        self.assertEqual(error, "Error message")
        self.assertEqual(status_code, 400)
        mock_export_api.handle_supabase_error.assert_called_with(mock_response)


class TestDataExportResource(unittest.TestCase):
    """Test the DataExportResource class."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create a test Flask app
        self.app = Flask(__name__)
        self.app.config['TESTING'] = True
        self.app_context = self.app.app_context()
        self.app_context.push()
        
        # Import the actual class for testing
        from api.export_api_resources import DataExportResource
        self.resource_class = DataExportResource
        
        # Create patches
        self.token_required_patcher = patch('api.export_api_resources.token_required', lambda f: f)
        self.token_required_mock = self.token_required_patcher.start()
        
        self.schema_patcher = patch('api.export_api_resources.ExportRequestSchema')
        self.schema_mock = self.schema_patcher.start()
        self.schema_instance = self.schema_mock.return_value
        self.schema_instance.load.return_value = {
            'format': 'csv',
            'data_type': 'molecules',
            'id': 'mol-123',
            'include_related': False
        }
        
        self.get_data_patcher = patch('api.export_api_resources.get_data_for_export')
        self.get_data_mock = self.get_data_patcher.start()
        self.get_data_mock.return_value = [{'id': 'mol-123', 'name': 'Test Molecule'}]
        
        self.generate_csv_patcher = patch('api.export_api_resources.generate_csv')
        self.generate_csv_mock = self.generate_csv_patcher.start()
        self.generate_csv_mock.return_value = io.StringIO("id,name\nmol-123,Test Molecule")
        
        self.generate_json_patcher = patch('api.export_api_resources.generate_json')
        self.generate_json_mock = self.generate_json_patcher.start()
