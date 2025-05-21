# ChEMBL Import Code Cleanup Plan

This document outlines the plan for consolidating the ChEMBL import code and removing unnecessary files.

## New Unified Files

These files constitute the new unified ChEMBL import solution:

1. `unified_chembl_import.py` - Main unified import script
2. `run_unified_chembl_import.sh` - Shell script to run the unified import
3. `CHEMBL_IMPORT_VERIFICATION.md` - Documentation for the import process

## Files to Keep (Still Useful)

These files provide valuable functionality that is not redundant with the unified import:

1. `test_chembl_with_db_public.py` - Tests for the ChEMBL import with db_public
2. `test_chembl_integrated_import.py` - Tests for the integrated import
3. `chembl_pubchem_resolver.py` - Used as a reference for the PubChem resolver implementation
4. `enhanced_property_calculator.py` - Used as a reference for the property calculator implementation

## Files to Remove (Redundant)

These files can be safely removed as their functionality is now incorporated into the unified import:

1. `ChEMBL_Integrated_Import.py` - Replaced by unified_chembl_import.py
2. `import_chembl_simplified.py` - Replaced by unified_chembl_import.py
3. `import_full_chembl.py` - Replaced by unified_chembl_import.py
4. `run_chembl_full_import.sh` - Replaced by run_unified_chembl_import.sh
5. `run_chembl_import.sh` - Replaced by run_unified_chembl_import.sh
6. `integrated_chembl_import_fix.py` - Integrated into unified_chembl_import.py
7. `fix_chembl_import_data.py` - Integrated into unified_chembl_import.py
8. `verify_chembl_import_data.py` - Integrated into unified_chembl_import.py

## Additional Cleanup Recommendations

1. Review test files to ensure they work with the unified import
2. Update documentation to refer to the new unified import script
3. Consider adding more test coverage for the unified script

## Cleanup Command

To execute the cleanup, run:

```bash
# Create backup directory
mkdir -p /home/mushu/Projects/CryoProtect/archive-chembl-files

# Move redundant files to backup
mv /home/mushu/Projects/CryoProtect/ChEMBL_Integrated_Import.py \
   /home/mushu/Projects/CryoProtect/import_chembl_simplified.py \
   /home/mushu/Projects/CryoProtect/import_full_chembl.py \
   /home/mushu/Projects/CryoProtect/run_chembl_full_import.sh \
   /home/mushu/Projects/CryoProtect/run_chembl_import.sh \
   /home/mushu/Projects/CryoProtect/integrated_chembl_import_fix.py \
   /home/mushu/Projects/CryoProtect/fix_chembl_import_data.py \
   /home/mushu/Projects/CryoProtect/verify_chembl_import_data.py \
   /home/mushu/Projects/CryoProtect/archive-chembl-files/
```

The files will be archived rather than deleted to ensure that no important code is lost. This allows for easy recovery if needed.