# Task 1.3: Organize Report Files

## Objective
Create a structured reports directory and organize all report files.

## Context
Report files are scattered throughout the codebase with timestamp suffixes. This creates confusion and makes it difficult to find the most recent reports. A proper organization system will improve maintainability and clarity.

## Acceptance Criteria
- Structured reports directory created with appropriate subdirectories
- All report files moved to relevant directories
- .gitignore updated to exclude archival reports
- Latest reports preserved for reference
- Documentation updated to reflect new report locations

## Implementation Steps

1. Create reports directory structure:
   ```bash
   mkdir -p reports/api/verification
   mkdir -p reports/rls/reports
   mkdir -p reports/rls/results
   mkdir -p reports/database
   mkdir -p reports/performance
   mkdir -p reports/archives/api
   mkdir -p reports/archives/rls
   mkdir -p reports/archives/database
   ```

2. Move API verification reports to appropriate directories:
   ```bash
   # Find and move API verification reports
   find . -name "api_verification_standalone_*.json" -exec mv {} ./reports/archives/api/ \;
   
   # Copy the latest API verification report as reference
   cp $(ls -t ./reports/archives/api/api_verification_standalone_*.json | head -1) ./reports/api/verification/latest.json
   ```

3. Move RLS verification reports to appropriate directories:
   ```bash
   # Find and move RLS reports
   find . -name "rls_verification_report_*.json" -exec mv {} ./reports/archives/rls/reports/ \;
   find . -name "rls_verification_results_*.json" -exec mv {} ./reports/archives/rls/results/ \;
   
   # Copy the latest RLS reports as reference
   cp $(ls -t ./reports/archives/rls/reports/rls_verification_report_*.json | head -1) ./reports/rls/reports/latest.json
   cp $(ls -t ./reports/archives/rls/results/rls_verification_results_*.json | head -1) ./reports/rls/results/latest.json
   ```

4. Move schema standardization reports:
   ```bash
   mkdir -p reports/archives/schema
   find . -name "schema_standardization_report_*.json" -exec mv {} ./reports/archives/schema/ \;
   ```

5. Move performance reports:
   ```bash
   mkdir -p reports/archives/performance
   find . -name "performance_improvements_*.json" -exec mv {} ./reports/archives/performance/ \;
   find . -name "performance_improvements_report_*.txt" -exec mv {} ./reports/archives/performance/ \;
   ```

6. Create README files for each directory explaining the purpose:
   ```bash
   cat > reports/README.md << 'EOF'
   # CryoProtect Reports
   
   This directory contains organized reports from various verification and testing processes.
   
   ## Structure
   
   - `api/` - API-related reports
   - `rls/` - Row Level Security verification reports
   - `database/` - Database verification reports
   - `performance/` - Performance testing reports
   - `archives/` - Historical reports (not tracked in git)
   
   ## Latest Reports
   
   The latest version of each report type is stored as `latest.json` in its respective directory.
   EOF
   
   # Create README for archives
   cat > reports/archives/README.md << 'EOF'
   # Report Archives
   
   This directory contains historical reports that are not tracked in git.
   These are preserved for reference but should not be committed to the repository.
   EOF
   ```

7. Update .gitignore to exclude archive directories:
   ```bash
   echo "# Ignore archived reports" >> .gitignore
   echo "reports/archives/" >> .gitignore
   ```

8. Create a script to regenerate latest reports from archives if needed:
   ```bash
   cat > reports/update_latest_reports.sh << 'EOF'
   #!/bin/bash
   # Script to update latest.json files from archives
   
   # API verification
   cp $(ls -t ./archives/api/api_verification_standalone_*.json 2>/dev/null | head -1) ./api/verification/latest.json 2>/dev/null
   
   # RLS reports
   cp $(ls -t ./archives/rls/reports/rls_verification_report_*.json 2>/dev/null | head -1) ./rls/reports/latest.json 2>/dev/null
   cp $(ls -t ./archives/rls/results/rls_verification_results_*.json 2>/dev/null | head -1) ./rls/results/latest.json 2>/dev/null
   
   echo "Latest reports updated from archives"
   EOF
   
   chmod +x reports/update_latest_reports.sh
   ```

9. Commit the changes:
   ```bash
   git add reports/
   git add .gitignore
   git commit -m "Organize report files into structured directory"
   ```

## Files to Modify
- Create reports directory structure
- Move all verification report files to appropriate directories
- .gitignore (add exclusion patterns)
- Create README files for documentation

## Verification
1. Check that reports are organized in their correct directories:
   ```bash
   find reports/ -name "*.json" | sort
   ```

2. Verify .gitignore is updated:
   ```bash
   grep "reports/archives" .gitignore
   ```

3. Test that archived reports are not tracked:
   ```bash
   touch reports/archives/test.json
   git status # Should not show test.json
   rm reports/archives/test.json # Clean up test file
   ```

4. Verify latest reference reports are accessible and tracked:
   ```bash
   ls -la reports/api/verification/latest.json
   ls -la reports/rls/reports/latest.json
   ls -la reports/rls/results/latest.json
   ```

## Notes for Roo Code Agent
- Be careful not to lose any report files during the reorganization
- Make sure the latest report copies are valid json files before committing
- If no reports are found for a particular type, create empty placeholder directories
- Update any scripts that output reports to use the new directory structure