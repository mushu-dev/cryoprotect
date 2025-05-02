# Team Models Consolidation Plan

## Current Situation

The team models functionality is currently split across 7 files:

```
team_models.py             -  3,746 bytes
team_models_combined.py    - 32,152 bytes
team_models_part2.py       -  3,227 bytes
team_models_part3.py       -  8,158 bytes
team_models_part4.py       -  6,673 bytes
team_models_part5.py       -  5,163 bytes
team_models_part6.py       -  5,175 bytes
```

This causes confusion and makes maintenance difficult.

## Recommendation

Keep `team_models_combined.py` as it appears to be the most complete version, and remove the fragmented part files.

## Immediate Actions

1. Verify that `team_models_combined.py` contains all functionality from other files
2. Update any imports in other files to reference `team_models_combined.py`
3. Remove redundant files

## Commands

```bash
# Rename team_models_combined.py to team_models.py (after backing up the original)
mv ./api/team_models.py ./api/team_models.py.original
mv ./api/team_models_combined.py ./api/team_models.py

# Remove part files after verification
git rm ./api/team_models_part2.py
git rm ./api/team_models_part3.py
git rm ./api/team_models_part4.py
git rm ./api/team_models_part5.py
git rm ./api/team_models_part6.py
git rm ./api/team_models.py.original

# Commit changes
git commit -m "Consolidate team models into a single file"
```