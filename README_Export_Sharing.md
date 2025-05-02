# CryoProtect Analyzer - Export and Sharing Features

This document provides information on how to use the export and sharing features of the CryoProtect Analyzer application.

## Table of Contents

1. [Overview](#overview)
2. [Data Export](#data-export)
3. [Visualization Export](#visualization-export)
4. [Report Generation](#report-generation)
5. [Sharing Results](#sharing-results)
6. [API Reference](#api-reference)
7. [Database Schema](#database-schema)
8. [Troubleshooting](#troubleshooting)

## Overview

The CryoProtect Analyzer provides comprehensive export and sharing capabilities that allow users to:

- Export data in various formats (CSV, JSON, Excel, PDF)
- Export visualizations (PNG, SVG, PDF)
- Generate comprehensive reports
- Share results via links, email, or embedding

These features are designed to enhance collaboration and make it easier to incorporate CryoProtect Analyzer results into external documents and workflows.

## Data Export

### Supported Data Types

The following data types can be exported:

- Molecules
- Mixtures
- Predictions
- Experiments
- Comparisons

### Supported Formats

The following formats are supported for data export:

- CSV (Comma-Separated Values)
- JSON (JavaScript Object Notation)
- Excel (XLSX)
- PDF (Portable Document Format)

### How to Export Data

1. Navigate to the page containing the data you want to export (e.g., Molecules, Mixtures, etc.)
2. Click the "Export" button
3. Select the export format (CSV, JSON, Excel, PDF)
4. Configure any format-specific options
5. Click "Export" to download the file

### Including Related Data

When exporting mixtures, you can choose to include related data such as predictions and experiments. This option is available in the export dialog.

## Visualization Export

### Supported Chart Types

The following chart types can be exported:

- Property Comparison (Bar Chart)
- Mixture Composition (Pie Chart)
- Prediction Accuracy (Scatter Plot)
- Custom Charts

### Supported Formats

The following formats are supported for visualization export:

- PNG (Portable Network Graphics)
- SVG (Scalable Vector Graphics)
- PDF (Portable Document Format)

### How to Export Visualizations

1. Navigate to the page containing the visualization you want to export
2. Click the "Export Visualization" button
3. Select the export format (PNG, SVG, PDF)
4. Configure any format-specific options (width, height, etc.)
5. Click "Export" to download the file

## Report Generation

### Report Sections

Reports can include the following sections:

- Title and Introduction
- Data Tables
- Visualizations
- Analysis and Conclusions
- References

### How to Generate Reports

1. Navigate to the "Reports" page
2. Click "New Report"
3. Enter a title for the report
4. Add sections to the report
5. Configure each section with content, data, and visualizations
6. Click "Generate Report" to create and download the PDF

### Report Templates

You can save report templates for future use. To create a template:

1. Configure a report as desired
2. Click "Save as Template"
3. Enter a name and description for the template
4. Click "Save"

To use a template:

1. Click "New Report"
2. Select "Use Template"
3. Choose a template from the list
4. Modify as needed
5. Click "Generate Report"

## Sharing Results

### Sharing Methods

The following sharing methods are available:

- Link Sharing: Generate a URL that can be shared with others
- Email Sharing: Send results directly via email
- Embedding: Generate HTML code to embed results in external websites

### Security Options

The following security options are available for shared content:

- Password Protection: Require a password to access shared content
- Expiration: Set an expiration date for shared content
- Access Tracking: Track who has accessed shared content

### How to Share Results

1. Navigate to the page containing the data you want to share
2. Click the "Share" button
3. Select the sharing method (Link, Email, Embed)
4. Configure security options if desired
5. Click "Share" to generate the link, send the email, or get the embed code

## API Reference

### Export Endpoints

- `POST /api/v1/export`: Export data in various formats
- `POST /api/v1/export/visualization`: Export visualizations
- `POST /api/v1/export/report`: Generate reports

### Sharing Endpoints

- `POST /api/v1/share`: Share results
- `GET /api/v1/share/:id`: Access shared content
- `POST /api/v1/share/:id`: Access password-protected shared content

### Request and Response Examples

#### Export Data

Request:
```json
{
  "format": "csv",
  "data_type": "mixtures",
  "id": "123e4567-e89b-12d3-a456-426614174000",
  "include_related": true
}
```

Response: CSV file download

#### Export Visualization

Request:
```json
{
  "format": "png",
  "chart_type": "mixture_composition",
  "data": {
    "components": [
      {"name": "Component 1", "concentration": 30},
      {"name": "Component 2", "concentration": 70}
    ]
  },
  "width": 800,
  "height": 600,
  "title": "Mixture Composition",
  "style": "default"
}
```

Response: PNG file download

#### Generate Report

Request:
```json
{
  "title": "CryoProtect Analysis Report",
  "sections": [
    {
      "title": "Introduction",
      "content": "This report presents an analysis of cryoprotectant mixtures..."
    },
    {
      "title": "Results",
      "content": "The following results were obtained...",
      "data": [...],
      "visualization": {
        "chart_type": "property_comparison",
        "data": {...}
      }
    }
  ],
  "include_visualizations": true,
  "include_data_tables": true
}
```

Response: PDF file download

#### Share Results

Request:
```json
{
  "share_type": "link",
  "data_type": "mixtures",
  "id": "123e4567-e89b-12d3-a456-426614174000",
  "password_protected": true,
  "password": "securepassword",
  "expiration": 86400
}
```

Response:
```json
{
  "share_id": "abcdef12-3456-7890-abcd-ef1234567890",
  "share_url": "https://example.com/share/abcdef12-3456-7890-abcd-ef1234567890",
  "password_protected": true,
  "expiration": "2023-04-16T18:17:40Z"
}
```

## Database Schema

The export and sharing functionality uses the following database tables:

- `shares`: Stores information about shared items
- `report_templates`: Stores report templates
- `visualization_templates`: Stores visualization templates
- `export_logs`: Logs export operations
- `share_access_logs`: Logs access to shared items

For detailed schema information, see the [database migration file](migrations/004_export_sharing_schema.sql).

## Troubleshooting

### Common Issues

#### Export fails with "No data available"

- Ensure that the data you're trying to export exists
- Check that you have permission to access the data
- Verify that the data type and ID are correct

#### Visualization export produces blank image

- Ensure that the chart data is properly formatted
- Check that the chart type is supported
- Try adjusting the width and height parameters

#### Report generation fails

- Ensure that all required sections have content
- Check that any referenced data exists
- Verify that visualizations are properly configured

#### Shared link is not accessible

- Ensure that the share has not expired
- Check that the password is correct (if password-protected)
- Verify that the share ID is correct

### Getting Help

If you encounter issues with the export and sharing features, please:

1. Check this documentation for solutions
2. Look for error messages in the browser console
3. Check the server logs for more detailed error information
4. Contact support with a detailed description of the issue