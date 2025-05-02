"""
Script to measure the coverage of the mixture_analysis.py module.
"""

import os
import sys
import coverage
import subprocess

def measure_coverage():
    """Measure the coverage of the mixture_analysis.py module."""
    # Create a coverage object
    cov = coverage.Coverage(source=['api.mixture_analysis'])
    
    # Start measuring coverage
    cov.start()
    
    try:
        # Import the module
        from api.mixture_analysis import (
            MixtureProperty, MixtureCompatibility, MixtureSynergy, 
            MixtureOptimization, MixtureRecommendation
        )
        
        # Create mock data
        components = [
            {"molecule_id": "mol-1", "concentration": 60, "concentration_unit": "%"},
            {"molecule_id": "mol-2", "concentration": 40, "concentration_unit": "%"}
        ]
        
        # Try to exercise some of the code
        try:
            # MixtureProperty
            MixtureProperty.predict_mixture_properties(components)
        except:
            pass
        
        try:
            # MixtureCompatibility
            MixtureCompatibility.analyze_compatibility(components)
        except:
            pass
        
        try:
            # MixtureSynergy
            MixtureSynergy.analyze_synergy("mix-123")
        except:
            pass
        
        try:
            # MixtureOptimization
            MixtureOptimization.optimize_composition("mix-123")
        except:
            pass
        
        try:
            # MixtureRecommendation
            MixtureRecommendation.analyze_mixture("mix-123")
        except:
            pass
        
    except Exception as e:
        print(f"Error importing or using the module: {e}")
    
    # Stop measuring coverage
    cov.stop()
    
    # Generate a report
    cov.save()
    
    # Print the report
    print("\nCoverage Report:")
    cov.report()
    
    # Generate HTML report
    cov.html_report(directory='reports/mixture_analysis_coverage')
    
    print(f"\nHTML report generated in reports/mixture_analysis_coverage")
    
    return cov


if __name__ == "__main__":
    # Create reports directory if it doesn't exist
    os.makedirs('reports/mixture_analysis_coverage', exist_ok=True)
    
    # Measure coverage
    cov = measure_coverage()
    
    # Run the existing tests to see if they improve coverage
    print("\nRunning existing tests...")
    try:
        subprocess.run(
            ["python", "-m", "pytest", "tests/test_mixture_analysis.py", "-v"],
            check=False
        )
    except Exception as e:
        print(f"Error running tests: {e}")
    
    # Generate a final report
    print("\nFinal Coverage Report:")
    cov.report()
    
    # Generate HTML report
    cov.html_report(directory='reports/mixture_analysis_coverage')
    
    print(f"\nFinal HTML report generated in reports/mixture_analysis_coverage")