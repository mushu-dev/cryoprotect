#!/usr/bin/env python3
"""
CryoProtect - Visualize ChEMBL Import Performance

This script generates visualization charts for ChEMBL import performance data.
It reads the JSON report file and creates performance charts.

Usage:
    python visualize_chembl_import.py [--report-file report.json]
"""

import os
import sys
import json
import argparse
import logging
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Try to import visualization libraries
try:
    import matplotlib.pyplot as plt
    import numpy as np
    HAS_VISUALIZATION = True
except ImportError:
    logger.warning("Matplotlib not available. Install with: pip install matplotlib")
    HAS_VISUALIZATION = False

def load_report(filename):
    """Load report data from JSON file."""
    try:
        with open(filename, 'r') as f:
            return json.load(f)
    except Exception as e:
        logger.error(f"Failed to load report: {str(e)}")
        return None

def visualize_performance(report_data, output_dir=None):
    """Generate performance visualizations from report data."""
    if not HAS_VISUALIZATION:
        logger.error("Cannot generate visualizations: matplotlib not available")
        return False
    
    if not output_dir:
        output_dir = os.path.dirname(os.path.abspath(__file__))
    
    try:
        # Extract performance metrics if available
        if "performance_metrics" not in report_data:
            logger.error("No performance metrics found in report")
            return False
        
        metrics = report_data["performance_metrics"]
        if not metrics:
            logger.error("Performance metrics list is empty")
            return False
        
        # Create output directory if needed
        os.makedirs(output_dir, exist_ok=True)
        
        # Extract data for plotting
        timestamps = [m["timestamp"] for m in metrics]
        overall_rates = [m["overall_rate"] for m in metrics]
        current_rates = [m["current_rate"] for m in metrics]
        success_rates = [m["success_rate"] * 100 for m in metrics]  # Convert to percentage
        processed = [m["processed"] for m in metrics]
        
        # Plot 1: Import Rate Over Time
        plt.figure(figsize=(10, 6))
        plt.plot(timestamps, overall_rates, 'b-', label='Overall Rate')
        plt.plot(timestamps, current_rates, 'r-', label='Current Rate')
        plt.xlabel('Time (seconds)')
        plt.ylabel('Molecules per Second')
        plt.title('ChEMBL Import Rate Over Time')
        plt.legend()
        plt.grid(True)
        rate_chart_path = os.path.join(output_dir, "chembl_import_rate.png")
        plt.savefig(rate_chart_path)
        logger.info(f"Import rate chart saved to {rate_chart_path}")
        
        # Plot 2: Success Rate Over Time
        plt.figure(figsize=(10, 6))
        plt.plot(processed, success_rates, 'g-')
        plt.xlabel('Molecules Processed')
        plt.ylabel('Success Rate (%)')
        plt.title('ChEMBL Import Success Rate')
        plt.grid(True)
        success_chart_path = os.path.join(output_dir, "chembl_import_success_rate.png")
        plt.savefig(success_chart_path)
        logger.info(f"Success rate chart saved to {success_chart_path}")
        
        # Plot 3: Cumulative Progress
        plt.figure(figsize=(10, 6))
        plt.plot(timestamps, processed, 'b-')
        
        # Add predicted line
        if len(timestamps) > 2 and len(processed) > 2:
            # Simple linear prediction
            x = np.array(timestamps)
            y = np.array(processed)
            
            # Get last point
            last_t = timestamps[-1]
            last_p = processed[-1]
            
            # Calculate average rate from the last 30% of the data
            start_idx = int(len(timestamps) * 0.7)
            if start_idx < len(timestamps):
                avg_rate = (processed[-1] - processed[start_idx]) / (timestamps[-1] - timestamps[start_idx])
                
                # Project to 5000 molecules
                if avg_rate > 0:
                    target_molecules = report_data["test_info"]["molecule_count"]
                    time_to_completion = (target_molecules - last_p) / avg_rate
                    
                    # Plot prediction line
                    pred_x = [last_t, last_t + time_to_completion]
                    pred_y = [last_p, target_molecules]
                    plt.plot(pred_x, pred_y, 'r--', label='Prediction')
                    
                    # Add annotation
                    plt.annotate(f'Est. completion: {(last_t + time_to_completion)/60:.1f} min',
                                xy=(last_t + time_to_completion/2, (last_p + target_molecules)/2),
                                xytext=(0, 30), textcoords='offset points',
                                arrowprops=dict(arrowstyle='->'))
        
        plt.xlabel('Time (seconds)')
        plt.ylabel('Molecules Processed')
        plt.title('ChEMBL Import Progress')
        plt.grid(True)
        progress_chart_path = os.path.join(output_dir, "chembl_import_progress.png")
        plt.savefig(progress_chart_path)
        logger.info(f"Progress chart saved to {progress_chart_path}")
        
        # Plot 4: Summary Chart
        plt.figure(figsize=(12, 8))
        
        # Test info
        molecule_count = report_data["test_info"]["molecule_count"]
        success_count = report_data["import_stats"]["success_count"]
        error_count = report_data["import_stats"]["error_count"]
        
        # Create pie chart
        plt.subplot(2, 2, 1)
        plt.pie([success_count, error_count], 
                labels=['Success', 'Error'],
                autopct='%1.1f%%',
                colors=['#4CAF50', '#F44336'])
        plt.title('Import Results')
        
        # Plot rate distribution
        plt.subplot(2, 2, 2)
        plt.hist(current_rates, bins=10, color='#2196F3')
        plt.xlabel('Molecules per Second')
        plt.ylabel('Frequency')
        plt.title('Rate Distribution')
        
        # Plot progression over time
        plt.subplot(2, 2, 3)
        plt.plot(timestamps, [p/molecule_count*100 for p in processed], 'g-')
        plt.xlabel('Time (seconds)')
        plt.ylabel('Completion (%)')
        plt.title('Import Progression')
        plt.grid(True)
        
        # Add performance summary as text
        plt.subplot(2, 2, 4)
        total_time = report_data["import_stats"]["total_time"]
        avg_rate = report_data["import_stats"]["molecules_per_second"]
        
        text = (
            f"Total Molecules: {molecule_count}\n"
            f"Success Rate: {success_count/molecule_count*100:.1f}%\n"
            f"Total Time: {total_time:.1f} seconds\n"
            f"Average Rate: {avg_rate:.2f} mol/sec\n\n"
        )
        
        if "performance_summary" in report_data:
            summary = report_data["performance_summary"]
            text += (
                f"Min Rate: {summary['min_rate']:.2f} mol/sec\n"
                f"Max Rate: {summary['max_rate']:.2f} mol/sec\n"
                f"Rate Stability: {summary['rate_stability']*100:.1f}%\n"
            )
            
        if "analysis" in report_data:
            analysis = report_data["analysis"]
            text += (
                f"\nTime for 5000 molecules: {analysis['estimated_full_import_time']/60:.1f} min\n"
                f"Time per 1000 molecules: {analysis['time_per_1000_molecules']:.1f} sec\n"
            )
        
        plt.text(0.1, 0.5, text, fontsize=10, verticalalignment='center')
        plt.title('Performance Summary')
        plt.axis('off')
        
        # Add super title
        plt.suptitle(f'ChEMBL Import Performance ({Path(output_dir).name})', fontsize=16)
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        
        # Save the figure
        summary_chart_path = os.path.join(output_dir, "chembl_import_summary.png")
        plt.savefig(summary_chart_path)
        logger.info(f"Summary chart saved to {summary_chart_path}")
        
        return True
    
    except Exception as e:
        logger.error(f"Failed to generate visualizations: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

def main():
    parser = argparse.ArgumentParser(description="Visualize ChEMBL import performance")
    parser.add_argument("--report-file", help="Path to the report JSON file")
    parser.add_argument("--output-dir", help="Directory to save visualizations")
    args = parser.parse_args()
    
    # Find the latest report if not specified
    report_file = args.report_file
    if not report_file:
        # Look for report files
        import glob
        report_files = glob.glob("chembl_import_report_*.json")
        if report_files:
            # Sort by modification time
            report_file = sorted(report_files, key=os.path.getmtime)[-1]
            logger.info(f"Using latest report file: {report_file}")
        else:
            logger.error("No report file specified and no report files found")
            return 1
    
    # Load report data
    report_data = load_report(report_file)
    if not report_data:
        return 1
    
    # Determine output directory
    output_dir = args.output_dir
    if not output_dir:
        output_dir = os.path.join(os.path.dirname(report_file), "visualizations")
    
    # Generate visualizations
    if not HAS_VISUALIZATION:
        logger.error("Cannot generate visualizations: matplotlib not available")
        logger.info("Install matplotlib with: pip install matplotlib")
        return 1
    
    if visualize_performance(report_data, output_dir):
        logger.info(f"All visualizations saved to {output_dir}")
        return 0
    else:
        return 1

if __name__ == "__main__":
    sys.exit(main())