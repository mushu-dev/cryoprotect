# Property Explorer User Guide

The CryoProtect Property Explorer is a powerful tool that allows you to analyze, visualize, and compare molecular properties across the entire database. This guide will walk you through using the various features of the Property Explorer.

## 1. Getting Started

Access the Property Explorer by clicking on the "Property Explorer" link in the main navigation menu. The tool provides four main sections:

1. **Property Explorer** - The main dashboard for exploring and visualizing properties
2. **Comparison Tool** - For side-by-side comparison of molecules
3. **Correlation Analysis** - For analyzing relationships between different properties
4. **Export Data** - For exporting property data in various formats

## 2. Exploring Properties

The Property Explorer dashboard displays key statistics about the molecular database, including:

- Total molecules in the database
- Number of available properties
- Identified cryoprotectants
- Overall data coverage

### 2.1 Filtering Properties

Use the filter panel on the left side to narrow down the molecules displayed:

- **Property Types** - Select categories of properties (Physical, Topological, etc.)
- **Value Ranges** - Filter molecules based on specific property value ranges
- **Classification** - Focus on specific molecules like verified cryoprotectants
- **Structure** - Search by molecular substructure using SMARTS patterns

Click "Apply Filters" to update the visualization and results.

### 2.2 Visualizing Properties

The Property Visualization panel allows you to:

1. Select which properties to display using the checkboxes
2. Choose a visualization type:
   - **Histogram** - Shows the distribution of property values
   - **Scatter Plot** - Compares two properties against each other
   - **Box Plot** - Displays statistical summaries of property distributions
   - **Radar Chart** - Multi-dimensional visualization of several properties

The visualizations are interactive - hover over data points for more information.

### 2.3 Molecule Results

The Molecule Results grid shows molecules matching your current filters:

- Click on a molecule name to view its detailed information
- Use the "View" button to open the molecular structure viewer
- Use the "Compare" button to add a molecule to the comparison tool
- Search within results using the search box
- Export your filtered results using the "Export Results" button

## 3. Comparing Molecules

The Comparison Tool allows for detailed side-by-side comparison of up to four molecules:

1. Select molecules using the dropdown or drag them from the results grid
2. Choose which properties to include in the comparison
3. Select a comparison visualization type:
   - **Side-by-Side** - Direct numerical comparison
   - **Radar Chart** - Visual comparison across multiple properties
   - **Difference Analysis** - Highlights differences between molecules

Click "Generate Comparison" to create the visualization.

## 4. Correlation Analysis

The Correlation Analysis tool helps identify relationships between different molecular properties:

1. Select properties to include in the correlation analysis
2. Choose a correlation method:
   - **Pearson** - Linear correlation (default)
   - **Spearman** - Rank-based correlation
   - **Kendall** - Another rank correlation measure

The tool generates:
- A correlation matrix showing the strength of relationships between properties
- Detailed analysis with statistical significance information
- Insights about key correlations

Use the dataset selector to focus on specific molecule groups.

## 5. Exporting Data

The Export Data tool allows you to export property data in various formats:

1. Select which dataset to export:
   - All molecules
   - Cryoprotectants only
   - Current filtered set
   - Custom selection

2. Choose which properties to include in the export

3. Select an export format:
   - CSV (for spreadsheet applications)
   - Excel
   - JSON (for programmatic use)
   - SDF (for chemical software)

The preview grid shows a sample of the data to be exported.

## 6. Tips for Effective Use

- **Start broad, then refine**: Begin with minimal filters, then add constraints as you explore
- **Use correlation analysis** to identify properties that might predict cryoprotection efficacy
- **Combine property ranges** to find molecules with specific profiles
- **Export results** for further analysis in specialized software
- **Save important comparisons** for reference in your research notes

## 7. API Access

All property explorer data is also available via the API, allowing programmatic access:

- `/api/v1/stats/dashboard` - Get database statistics
- `/api/v1/properties/types` - Get property types
- `/api/v1/properties/data` - Access property data
- `/api/v1/properties/statistics` - Get statistical summaries of properties
- `/api/v1/properties/correlation` - Calculate correlations between properties
- `/api/v1/properties/export` - Generate data exports

For API documentation, refer to the OpenAPI documentation at `/api/v1/docs`.

---

For further assistance or to report issues, please contact the CryoProtect support team.