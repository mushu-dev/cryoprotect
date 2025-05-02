#!/usr/bin/env python3
"""
Data Quality Report Generator

This script generates a comprehensive report on the quality and completeness
of data in the CryoProtect database, based on verification results from
verify_imported_data.py.
"""

import os
import sys
import json
import logging
import argparse
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
from typing import Dict, List, Any, Optional, Tuple

# Set up logging
os.makedirs('logs', exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('logs/data_quality_report.log')
    ]
)

logger = logging.getLogger(__name__)

# Status constants
STATUS_SUCCESS = "SUCCESS"
STATUS_WARNING = "WARNING"
STATUS_FAILURE = "FAILURE"

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Generate data quality report")
    parser.add_argument("input_file", help="Input JSON file with verification results")
    parser.add_argument("--output", default="reports", help="Output directory for report")
    parser.add_argument("--template", help="Path to report template file")
    parser.add_argument("--format", choices=["md", "html", "both"], default="md",
                       help="Output format for report")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    return parser.parse_args()

class DataQualityReport:
    """Data quality report generator"""
    
    def __init__(self, input_file: str, output_dir: str = "reports", template_file: Optional[str] = None):
        """Initialize report generator"""
        self.input_file = input_file
        self.output_dir = output_dir
        self.template_file = template_file
        self.verification_data = {}
        self.plots_dir = os.path.join(output_dir, "plots")
        
        # Create output directories if they don't exist
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(self.plots_dir, exist_ok=True)
        
    def load_verification_data(self) -> bool:
        """Load verification data from input file"""
        try:
            logger.info(f"Loading verification data from {self.input_file}")
            with open(self.input_file, 'r') as f:
                self.verification_data = json.load(f)
            
            # Basic validation of verification data
            required_keys = ["timestamp", "data_source", "assessments"]
            for key in required_keys:
                if key not in self.verification_data:
                    logger.error(f"Missing required key '{key}' in verification data")
                    return False
                    
            logger.info(f"Successfully loaded verification data for {self.verification_data.get('data_source', 'unknown')} source")
            return True
            
        except FileNotFoundError:
            logger.error(f"Input file not found: {self.input_file}")
            return False
        except json.JSONDecodeError:
            logger.error(f"Invalid JSON in input file: {self.input_file}")
            return False
        except Exception as e:
            logger.error(f"Error loading verification data: {e}")
            return False
    
    def generate_plots(self) -> Dict[str, str]:
        """Generate plots for the report"""
        plot_files = {}
        
        try:
            # 1. Property completeness plot
            if "property_completeness" in self.verification_data and self.verification_data["property_completeness"]:
                property_data = self.verification_data["property_completeness"]
                # Filter out non-property entries like 'average_coverage'
                properties = {k: v for k, v in property_data.items() if isinstance(v, dict) and "coverage_percent" in v}
                
                if properties:
                    # Create bar chart of property coverage
                    plt.figure(figsize=(10, 6))
                    property_names = list(properties.keys())
                    coverage_values = [prop["coverage_percent"] for prop in properties.values()]
                    
                    # Sort by coverage percentage
                    sorted_indices = np.argsort(coverage_values)
                    sorted_names = [property_names[i] for i in sorted_indices]
                    sorted_values = [coverage_values[i] for i in sorted_indices]
                    
                    # Create color gradient based on coverage values
                    colors = []
                    for value in sorted_values:
                        if value >= 90:
                            colors.append('green')
                        elif value >= 70:
                            colors.append('orange')
                        else:
                            colors.append('red')
                    
                    plt.barh(sorted_names, sorted_values, color=colors)
                    plt.xlabel('Coverage (%)')
                    plt.ylabel('Property')
                    plt.title('Property Completeness')
                    plt.xlim(0, 100)
                    plt.grid(axis='x', linestyle='--', alpha=0.7)
                    
                    # Add value labels to the bars
                    for i, v in enumerate(sorted_values):
                        plt.text(v + 1, i, f"{v:.1f}%", va='center')
                    
                    # Save the plot
                    plot_file = os.path.join(self.plots_dir, "property_completeness.png")
                    plt.tight_layout()
                    plt.savefig(plot_file)
                    plt.close()
                    
                    plot_files["property_completeness"] = plot_file
                    logger.info(f"Generated property completeness plot: {plot_file}")
            
            # 2. Null values plot
            if "null_checks" in self.verification_data and self.verification_data["null_checks"]:
                null_data = self.verification_data["null_checks"]
                
                if "molecules" in null_data:
                    molecules_data = null_data["molecules"]
                    # Extract percentage values
                    labels = []
                    values = []
                    
                    for key, value in molecules_data.items():
                        if key.endswith("_percent"):
                            field_name = key.replace("_percent", "").replace("null_", "")
                            labels.append(field_name)
                            values.append(value)
                    
                    if labels and values:
                        plt.figure(figsize=(8, 6))
                        
                        # Create color gradient based on null percentages (higher is worse)
                        colors = []
                        for value in values:
                            if value <= 5:
                                colors.append('green')
                            elif value <= 10:
                                colors.append('orange')
                            else:
                                colors.append('red')
                        
                        plt.bar(labels, values, color=colors)
                        plt.ylabel('Null Values (%)')
                        plt.title('Null Values in Molecule Fields')
                        plt.grid(axis='y', linestyle='--', alpha=0.7)
                        
                        # Add value labels to the bars
                        for i, v in enumerate(values):
                            plt.text(i, v + 0.5, f"{v:.1f}%", ha='center')
                        
                        # Save the plot
                        plot_file = os.path.join(self.plots_dir, "null_values.png")
                        plt.tight_layout()
                        plt.savefig(plot_file)
                        plt.close()
                        
                        plot_files["null_values"] = plot_file
                        logger.info(f"Generated null values plot: {plot_file}")
            
            # 3. Assessment summary plot
            if "assessments" in self.verification_data:
                assessments = self.verification_data["assessments"]
                
                if assessments:
                    # Count status types
                    status_counts = {
                        STATUS_SUCCESS: 0,
                        STATUS_WARNING: 0,
                        STATUS_FAILURE: 0,
                        "ERROR": 0
                    }
                    
                    for check, status in assessments.items():
                        if status in status_counts:
                            status_counts[status] += 1
                        else:
                            status_counts["ERROR"] += 1
                    
                    # Create pie chart
                    plt.figure(figsize=(8, 8))
                    labels = []
                    sizes = []
                    colors = []
                    
                    for status, count in status_counts.items():
                        if count > 0:
                            labels.append(f"{status} ({count})")
                            sizes.append(count)
                            if status == STATUS_SUCCESS:
                                colors.append('green')
                            elif status == STATUS_WARNING:
                                colors.append('orange')
                            elif status == STATUS_FAILURE:
                                colors.append('red')
                            else:
                                colors.append('gray')
                    
                    if sizes:
                        plt.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', 
                                startangle=90, shadow=True)
                        plt.axis('equal')
                        plt.title('Assessment Results')
                        
                        # Save the plot
                        plot_file = os.path.join(self.plots_dir, "assessment_summary.png")
                        plt.tight_layout()
                        plt.savefig(plot_file)
                        plt.close()
                        
                        plot_files["assessment_summary"] = plot_file
                        logger.info(f"Generated assessment summary plot: {plot_file}")
            
            return plot_files
            
        except Exception as e:
            logger.error(f"Error generating plots: {e}")
            return plot_files
    
    def generate_markdown_report(self, plot_files: Dict[str, str]) -> str:
        """Generate a Markdown report from verification data"""
        try:
            # Extract data for the report
            data_source = self.verification_data.get("data_source", "Unknown")
            timestamp = self.verification_data.get("timestamp", datetime.now().isoformat())
            try:
                report_date = datetime.fromisoformat(timestamp).strftime("%Y-%m-%d %H:%M:%S")
            except:
                report_date = timestamp
                
            overall_assessment = self.verification_data.get("overall_assessment", "INCOMPLETE")
            assessments = self.verification_data.get("assessments", {})
            
            # Start building the report
            report = []
            
            # Title and header
            report.append(f"# Data Quality Report: {data_source}")
            report.append(f"\nGenerated on: {report_date}")
            report.append("\n## Executive Summary")
            
            # Overall assessment
            status_emoji = "✅" if overall_assessment == STATUS_SUCCESS else "⚠️" if overall_assessment == STATUS_WARNING else "❌"
            report.append(f"\nOverall data quality assessment: **{overall_assessment}** {status_emoji}")
            
            # Basic statistics
            counts = self.verification_data.get("counts", {})
            if counts:
                report.append("\n### Database Statistics")
                report.append("\n| Metric | Count |")
                report.append("| ------ | ----- |")
                for metric, count in counts.items():
                    report.append(f"| {metric.replace('_', ' ').title()} | {count:,} |")
            
            # Assessment summary plot
            if "assessment_summary" in plot_files:
                report.append("\n### Assessment Summary")
                report.append("\nThe following chart shows the distribution of assessment results:")
                report.append(f"\n![Assessment Summary](plots/{os.path.basename(plot_files['assessment_summary'])})")
            
            # Detailed assessments
            report.append("\n## Detailed Assessments")
            
            for check, status in assessments.items():
                check_name = check.replace("_", " ").title()
                status_emoji = "✅" if status == STATUS_SUCCESS else "⚠️" if status == STATUS_WARNING else "❌"
                report.append(f"\n### {check_name}: {status} {status_emoji}")
                
                # Add specific details for each check type
                if check == "counts":
                    report.append("\nDatabase count verification results:")
                    counts_data = self.verification_data.get("counts", {})
                    if counts_data:
                        report.append("\n| Metric | Count |")
                        report.append("| ------ | ----- |")
                        for metric, count in counts_data.items():
                            report.append(f"| {metric.replace('_', ' ').title()} | {count:,} |")
                
                elif check == "null_checks":
                    report.append("\nNull value verification results:")
                    null_data = self.verification_data.get("null_checks", {})
                    
                    if "molecules" in null_data:
                        report.append("\n#### Molecules Table")
                        molecules_data = null_data["molecules"]
                        report.append("\n| Field | Null Count | Null Percentage |")
                        report.append("| ----- | ---------- | --------------- |")
                        
                        for field in ["name", "smiles", "chembl_id"]:
                            null_count = molecules_data.get(f"null_{field}", "N/A")
                            null_percent = molecules_data.get(f"null_{field}_percent", "N/A")
                            if null_percent != "N/A":
                                null_percent = f"{null_percent:.2f}%"
                            report.append(f"| {field.title()} | {null_count} | {null_percent} |")
                    
                    if "properties" in null_data:
                        report.append("\n#### Properties Table")
                        properties_data = null_data["properties"]
                        report.append("\n| Field | Null Count | Null Percentage |")
                        report.append("| ----- | ---------- | --------------- |")
                        
                        for field in ["value", "property_type"]:
                            null_count = properties_data.get(f"null_{field}", "N/A")
                            null_percent = properties_data.get(f"null_{field}_percent", "N/A")
                            if null_percent != "N/A":
                                null_percent = f"{null_percent:.2f}%"
                            report.append(f"| {field.title()} | {null_count} | {null_percent} |")
                    
                    # Add null values plot if available
                    if "null_values" in plot_files:
                        report.append("\n#### Null Values Visualization")
                        report.append(f"\n![Null Values](plots/{os.path.basename(plot_files['null_values'])})")
                
                elif check == "consistency_checks":
                    report.append("\nData consistency verification results:")
                    consistency_data = self.verification_data.get("consistency_checks", {})
                    
                    if consistency_data:
                        report.append("\n| Issue | Count |")
                        report.append("| ----- | ----- |")
                        
                        issues = {
                            "molecules_without_properties": "Molecules without properties",
                            "orphaned_properties": "Orphaned properties",
                            "duplicate_molecules": "Duplicate molecules"
                        }
                        
                        for key, label in issues.items():
                            if key in consistency_data:
                                report.append(f"| {label} | {consistency_data[key]} |")
                
                elif check == "property_completeness":
                    report.append("\nProperty completeness verification results:")
                    property_data = self.verification_data.get("property_completeness", {})
                    
                    # Extract average coverage if available
                    avg_coverage = property_data.get("average_coverage", "N/A")
                    if avg_coverage != "N/A":
                        report.append(f"\nAverage property coverage: **{avg_coverage:.2f}%**")
                    
                    # Filter out non-property entries like 'average_coverage'
                    properties = {k: v for k, v in property_data.items() if isinstance(v, dict) and "coverage_percent" in v}
                    
                    if properties:
                        report.append("\n| Property | Coverage | Molecules with Property | Total Molecules |")
                        report.append("| -------- | -------- | ----------------------- | --------------- |")
                        
                        # Sort properties by coverage percentage (descending)
                        sorted_properties = sorted(
                            properties.items(), 
                            key=lambda x: x[1].get("coverage_percent", 0), 
                            reverse=True
                        )
                        
                        for prop_name, prop_data in sorted_properties:
                            coverage = prop_data.get("coverage_percent", "N/A")
                            if coverage != "N/A":
                                coverage = f"{coverage:.2f}%"
                            
                            molecules_with_prop = prop_data.get("molecules_with_property", "N/A")
                            total_molecules = prop_data.get("total_molecules", "N/A")
                            
                            report.append(f"| {prop_name} | {coverage} | {molecules_with_prop} | {total_molecules} |")
                    
                    # Add property completeness plot if available
                    if "property_completeness" in plot_files:
                        report.append("\n#### Property Completeness Visualization")
                        report.append(f"\n![Property Completeness](plots/{os.path.basename(plot_files['property_completeness'])})")
                
                elif check == "reference_compounds":
                    report.append("\nReference compound verification results:")
                    ref_data = self.verification_data.get("reference_compounds", {})
                    
                    # Extract summary if available
                    if "summary" in ref_data:
                        summary = ref_data["summary"]
                        total = summary.get("total", "N/A")
                        present = summary.get("present", "N/A")
                        missing = summary.get("missing", "N/A")
                        success_rate = summary.get("success_rate", "N/A")
                        
                        if success_rate != "N/A":
                            success_rate = f"{success_rate*100:.2f}%"
                        
                        report.append(f"\n- Total reference compounds: {total}")
                        report.append(f"- Present in database: {present}")
                        report.append(f"- Missing from database: {missing}")
                        report.append(f"- Success rate: {success_rate}")
                        
                        # List missing compounds if there are any
                        if missing and missing != "N/A" and missing > 0:
                            report.append("\n#### Missing Reference Compounds")
                            
                            missing_compounds = []
                            for cid, details in ref_data.items():
                                if cid != "summary" and not details.get("present", False):
                                    missing_compounds.append(cid)
                            
                            if missing_compounds:
                                report.append("\n```")
                                for i, cid in enumerate(missing_compounds):
                                    report.append(cid + (", " if (i+1) % 5 != 0 and i < len(missing_compounds)-1 else "\n" if (i+1) % 5 == 0 else ""))
                                report.append("```")
            
            # Recommendations section
            report.append("\n## Recommendations")
            
            # Generate recommendations based on assessment results
            recommendations = []
            
            if assessments.get("counts") == STATUS_FAILURE:
                recommendations.append("- **Critical**: Investigate the low molecule count. The database may be missing significant data.")
            
            if assessments.get("null_checks") in [STATUS_WARNING, STATUS_FAILURE]:
                recommendations.append("- **High Priority**: Address the high percentage of null values in key fields.")
                null_data = self.verification_data.get("null_checks", {})
                if "molecules" in null_data:
                    molecules_data = null_data["molecules"]
                    for field in ["name", "smiles", "chembl_id"]:
                        null_percent = molecules_data.get(f"null_{field}_percent", 0)
                        if null_percent > 10:
                            recommendations.append(f"  - Fix missing {field} values (currently {null_percent:.2f}% null)")
            
            if assessments.get("consistency_checks") in [STATUS_WARNING, STATUS_FAILURE]:
                recommendations.append("- **Medium Priority**: Resolve data consistency issues:")
                consistency_data = self.verification_data.get("consistency_checks", {})
                if consistency_data.get("molecules_without_properties", 0) > 0:
                    recommendations.append(f"  - Add properties for {consistency_data['molecules_without_properties']} molecules that currently have none")
                if consistency_data.get("orphaned_properties", 0) > 0:
                    recommendations.append(f"  - Remove or fix {consistency_data['orphaned_properties']} orphaned properties")
                if consistency_data.get("duplicate_molecules", 0) > 0:
                    recommendations.append(f"  - Deduplicate {consistency_data['duplicate_molecules']} duplicate molecule entries")
            
            if assessments.get("property_completeness") in [STATUS_WARNING, STATUS_FAILURE]:
                recommendations.append("- **Medium Priority**: Improve property completeness:")
                property_data = self.verification_data.get("property_completeness", {})
                # Filter out non-property entries like 'average_coverage'
                properties = {k: v for k, v in property_data.items() if isinstance(v, dict) and "coverage_percent" in v}
                low_coverage_props = []
                for prop_name, prop_data in properties.items():
                    coverage = prop_data.get("coverage_percent", 100)
                    if coverage < 70:
                        low_coverage_props.append((prop_name, coverage))
                
                # Sort by coverage (ascending)
                low_coverage_props.sort(key=lambda x: x[1])
                
                for prop_name, coverage in low_coverage_props[:3]:  # Show top 3 lowest coverage properties
                    recommendations.append(f"  - Increase coverage for {prop_name} (currently only {coverage:.2f}%)")
            
            if assessments.get("reference_compounds") in [STATUS_WARNING, STATUS_FAILURE]:
                ref_data = self.verification_data.get("reference_compounds", {})
                if "summary" in ref_data:
                    summary = ref_data["summary"]
                    missing = summary.get("missing", 0)
                    if missing > 0:
                        recommendations.append(f"- **High Priority**: Import the {missing} missing reference compounds")
            
            # Add recommendations to report
            if recommendations:
                for recommendation in recommendations:
                    report.append(recommendation)
            else:
                report.append("\nNo specific recommendations at this time. Data quality is good.")
            
            # Next steps section
            report.append("\n## Next Steps")
            
            if overall_assessment == STATUS_SUCCESS:
                report.append("\nThe data quality is good. Recommended next steps:")
                report.append("1. Continue with regular data quality monitoring")
                report.append("2. Consider enhancing the database with additional properties")
                report.append("3. Implement automated data quality checks in the data pipeline")
            elif overall_assessment == STATUS_WARNING:
                report.append("\nThe data quality has some issues that should be addressed. Recommended next steps:")
                report.append("1. Address the recommendations listed above")
                report.append("2. Re-run the verification after fixes are implemented")
                report.append("3. Implement data quality gates in the import process")
            else:  # FAILURE or other
                report.append("\nThe data quality has critical issues that must be addressed. Recommended next steps:")
                report.append("1. Immediately address the critical recommendations listed above")
                report.append("2. Investigate the root causes of data quality issues")
                report.append("3. Re-run the verification after each fix to track progress")
                report.append("4. Review and potentially revise the data import process")
            
            # Join all parts of the report
            return "\n".join(report)
            
        except Exception as e:
            logger.error(f"Error generating Markdown report: {e}")
            return f"# Error Generating Report\n\nAn error occurred while generating the report: {e}"
    
    def generate_html_report(self, markdown_content: str) -> str:
        """Convert Markdown report to HTML"""
        try:
            import markdown
            
            # Convert Markdown to HTML
            html_content = markdown.markdown(markdown_content, extensions=['tables'])
            
            # Wrap in basic HTML template
            html_report = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Data Quality Report</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            line-height: 1.6;
            color: #333;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
        }}
        h1, h2, h3 {{
            color: #2c3e50;
        }}
        table {{
            border-collapse: collapse;
            width: 100%;
            margin-bottom: 20px;
        }}
        th, td {{
            border: 1px solid #ddd;
            padding: 8px;
            text-align: left;
        }}
        th {{
            background-color: #f2f2f2;
        }}
        tr:nth-child(even) {{
            background-color: #f9f9f9;
        }}
        img {{
            max-width: 100%;
            height: auto;
        }}
        .success {{
            color: green;
        }}
        .warning {{
            color: orange;
        }}
        .failure {{
            color: red;
        }}
    </style>
</head>
<body>
    {html_content}
</body>
</html>
"""
            return html_report
            
        except ImportError:
            logger.warning("Python-Markdown package not installed. HTML report generation skipped.")
            return ""
        except Exception as e:
            logger.error(f"Error generating HTML report: {e}")
            return ""
    
    def save_report(self, content: str, format: str) -> str:
        """Save report to file"""
        try:
            # Generate filename based on data source and timestamp
            data_source = self.verification_data.get("data_source", "unknown").lower()
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            
            filename = f"{data_source}_quality_report_{timestamp}.{format}"
            filepath = os.path.join(self.output_dir, filename)
            
            with open(filepath, 'w', encoding='utf-8') as f:
                f.write(content)
                
            logger.info(f"Report saved to {filepath}")
            return filepath
            
        except Exception as e:
            logger.error(f"Error saving report: {e}")
            return ""
    
    def generate_report(self, format: str = "md") -> Tuple[bool, List[str]]:
        """Generate and save the report"""
        try:
            # Load verification data
            if not self.load_verification_data():
                return False, []
            
            # Generate plots
            plot_files = self.generate_plots()
            
            # Generate Markdown report
            markdown_content = self.generate_markdown_report(plot_files)
            
            # Save reports
            saved_files = []
            
            if format in ["md", "both"]:
                md_file = self.save_report(markdown_content, "md")
                if md_file:
                    saved_files.append(md_file)
            
            if format in ["html", "both"]:
                html_content = self.generate_html_report(markdown_content)
                if html_content:
                    html_file = self.save_report(html_content, "html")
                    if html_file:
                        saved_files.append(html_file)
            
            return len(saved_files) > 0, saved_files
            
        except Exception as e:
            logger.error(f"Error generating report: {e}")
            return False, []

def main():
    """Main function"""
    # Parse command line arguments
    args = parse_arguments()
    
    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Create report generator
    report_generator = DataQualityReport(
        input_file=args.input_file,
        output_dir=args.output,
        template_file=args.template
    )
    
    # Generate report
    success, saved_files = report_generator.generate_report(args.format)
    
    if success:
        print("\n" + "=" * 60)
        print("Data Quality Report Generation Complete")
        print("=" * 60)
        print(f"Generated {len(saved_files)} report file(s):")
        for file in saved_files:
            print(f"- {file}")
        print("=" * 60)
        return 0
    else:
        print("\n" + "=" * 60)
        print("Data Quality Report Generation Failed")
        print("=" * 60)
        print("Check the logs for more information.")
        print("=" * 60)
        return 1

if __name__ == "__main__":
    sys.exit(main())