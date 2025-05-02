"""
Script to measure coverage for the comparisons module.
"""

import sys
import os
import coverage
import unittest
import importlib.util

# Start coverage measurement
cov = coverage.Coverage(source=['api.comparisons'])
cov.start()

# Import our test module
spec = importlib.util.spec_from_file_location("test_comparisons", "tests/test_comparisons.py")
test_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(test_module)

# Run the tests
suite = unittest.TestLoader().loadTestsFromModule(test_module)
result = unittest.TextTestRunner().run(suite)

# Stop coverage measurement and report
cov.stop()
cov.save()

print("\nCoverage Report:")
cov.report()

# Write coverage report to a file
with open('reports/comparisons_coverage_summary.md', 'w') as f:
    f.write("# Coverage Report for api/comparisons.py\n\n")
    f.write("## Summary\n\n")
    
    # Get the coverage data
    total_lines = 0
    covered_lines = 0
    
    for filename in cov.get_data().measured_files():
        if 'comparisons.py' in filename:
            analysis = cov._analyze(filename)
            total_lines = len(analysis.statements)
            missing_lines = len(analysis.missing)
            covered_lines = total_lines - missing_lines
            coverage_percent = (covered_lines / total_lines) * 100 if total_lines > 0 else 0
            
            f.write(f"- **File**: {os.path.basename(filename)}\n")
            f.write(f"- **Lines**: {total_lines}\n")
            f.write(f"- **Covered**: {covered_lines}\n")
            f.write(f"- **Missing**: {missing_lines}\n")
            f.write(f"- **Coverage**: {coverage_percent:.2f}%\n\n")
            
            f.write("## Missing Lines\n\n")
            if missing_lines > 0:
                f.write("```\n")
                for line in sorted(analysis.missing):
                    f.write(f"{line}\n")
                f.write("```\n")
            else:
                f.write("No missing lines! Full coverage achieved.\n")

print(f"\nDetailed coverage report written to reports/comparisons_coverage_summary.md")