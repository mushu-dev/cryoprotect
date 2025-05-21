#!/bin/bash
# Script to archive redundant ChEMBL import files

# Create backup directory
BACKUP_DIR="/home/mushu/Projects/CryoProtect/archive-chembl-files"
mkdir -p "$BACKUP_DIR"
echo "Created backup directory: $BACKUP_DIR"

# List of files to move to backup
FILES_TO_BACKUP=(
  "/home/mushu/Projects/CryoProtect/ChEMBL_Integrated_Import.py"
  "/home/mushu/Projects/CryoProtect/import_chembl_simplified.py"
  "/home/mushu/Projects/CryoProtect/import_full_chembl.py"
  "/home/mushu/Projects/CryoProtect/run_chembl_full_import.sh"
  "/home/mushu/Projects/CryoProtect/run_chembl_import.sh"
  "/home/mushu/Projects/CryoProtect/integrated_chembl_import_fix.py"
  "/home/mushu/Projects/CryoProtect/fix_chembl_import_data.py"
  "/home/mushu/Projects/CryoProtect/verify_chembl_import_data.py"
)

# Move files to backup directory
for file in "${FILES_TO_BACKUP[@]}"; do
  if [ -f "$file" ]; then
    echo "Moving $file to backup directory..."
    mv "$file" "$BACKUP_DIR/"
  else
    echo "Warning: $file does not exist, skipping"
  fi
done

echo ""
echo "Cleanup complete. The following files have been moved to $BACKUP_DIR:"
ls -la "$BACKUP_DIR"

echo ""
echo "The following new files have been created to replace them:"
echo "1. unified_chembl_import.py - Main unified import script"
echo "2. run_unified_chembl_import.sh - Shell script to run the unified import"
echo "3. UNIFIED_CHEMBL_IMPORT.md - Documentation for the unified import"
echo "4. chembl_cleanup_plan.md - Cleanup plan documentation"

echo ""
echo "If you need to restore any files, you can find them in the backup directory."