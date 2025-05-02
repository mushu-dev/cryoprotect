# Task 1.2: Consolidate Team Models

## Objective
Consolidate fragmented team model files into a single coherent module.

## Context
The team model functionality is split across 7 files creating confusion and maintenance difficulties:
- team_models.py (3,746 bytes)
- team_models_combined.py (32,152 bytes)
- team_models_part2.py (3,227 bytes)
- team_models_part3.py (8,158 bytes)
- team_models_part4.py (6,673 bytes)
- team_models_part5.py (5,163 bytes)
- team_models_part6.py (5,175 bytes)

## Acceptance Criteria
- All team model functionality exists in a single file
- No functionality is lost during consolidation
- Imports are updated in dependent files
- Redundant files are removed

## Implementation Steps

1. Create a backup of the original team_models.py:
   ```bash
   cp api/team_models.py api/team_models.py.backup
   ```

2. Verify that team_models_combined.py contains all functionality from the part files:
   - Analyze team_models_combined.py content
   - Compare with part files to ensure all classes, methods, and functions are included
   - Check for any unique functionality in the part files that might be missing

3. Replace team_models.py with team_models_combined.py:
   ```bash
   cp api/team_models_combined.py api/team_models.py
   ```

4. Find all files that import from the team models parts and update imports:
   ```bash
   grep -r "from api.team_models_part" . --include="*.py"
   ```
   For each file found, modify imports to point to the consolidated file.

5. Test that functionality works with the consolidated file by running relevant tests:
   ```bash
   python tests/test_team_models.py
   python tests/test_team_resources.py
   ```

6. Once tests pass, remove the redundant files:
   ```bash
   git rm api/team_models_combined.py
   git rm api/team_models_part2.py
   git rm api/team_models_part3.py
   git rm api/team_models_part4.py
   git rm api/team_models_part5.py
   git rm api/team_models_part6.py
   ```

7. Commit the changes:
   ```bash
   git add api/team_models.py
   git commit -m "Consolidate team models into a single file"
   ```

## Files to Modify
- api/team_models.py (replace content with team_models_combined.py)
- Any files that import from team_models_part*.py files
- Remove:
  - api/team_models_combined.py
  - api/team_models_part2.py
  - api/team_models_part3.py
  - api/team_models_part4.py
  - api/team_models_part5.py
  - api/team_models_part6.py

## Verification
1. Run team-related tests to ensure functionality is preserved:
   ```bash
   python tests/test_team_models.py
   python tests/test_team_resources.py
   ```
2. Check for any import errors or runtime errors that might indicate missing functionality
3. Verify each functionality previously provided by the separate files now works with the consolidated file

## Notes for Roo Code Agent
- Pay special attention to imports and dependencies between the models
- Some classes might reference other models or utilities
- If functionality tests don't exist, create basic test cases before consolidation to ensure functionality is preserved
- The team_models.py file should be comprehensive and maintain the same API as before