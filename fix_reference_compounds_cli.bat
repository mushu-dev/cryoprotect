@echo off
REM CryoProtect v2 - Fix Reference Compounds Script using Supabase CLI
REM This script fixes missing properties for reference compounds using the Supabase CLI

echo Starting reference compounds fix using Supabase CLI...

REM Set variables
set PROJECT_ID=tsdlmynydfuypiugmkev

REM Check if user is logged in to Supabase
supabase projects list >nul 2>&1
if %ERRORLEVEL% neq 0 (
    echo You are not logged in to Supabase CLI. Please login first.
    echo Run: supabase login
    exit /b 1
)

echo Verifying reference compounds exist...

REM Create a temporary SQL file
echo SELECT chembl_id, id FROM molecules > verify_ref.sql
echo WHERE chembl_id IN ('CHEMBL1098659', 'CHEMBL1487', 'CHEMBL262548', 'CHEMBL388978', >> verify_ref.sql
echo                    'CHEMBL500033', 'CHEMBL6196', 'CHEMBL66195', 'CHEMBL6752', 'CHEMBL967'); >> verify_ref.sql

REM Execute the SQL
supabase db execute --project-ref %PROJECT_ID% < verify_ref.sql
if %ERRORLEVEL% neq 0 (
    echo Error verifying reference compounds. Exiting.
    del verify_ref.sql
    exit /b 1
)
del verify_ref.sql

echo Ensuring property types exist...

REM Create a temporary SQL file for property types
echo -- Create property types if they don't exist > property_types.sql
echo INSERT INTO property_types (name, data_type, description) >> property_types.sql
echo SELECT 'logP', 'numeric', 'Partition coefficient' >> property_types.sql
echo WHERE NOT EXISTS (SELECT 1 FROM property_types WHERE name = 'logP'); >> property_types.sql
echo. >> property_types.sql
echo INSERT INTO property_types (name, data_type, description) >> property_types.sql
echo SELECT 'h_bond_donors', 'numeric', 'Number of hydrogen bond donors' >> property_types.sql
echo WHERE NOT EXISTS (SELECT 1 FROM property_types WHERE name = 'h_bond_donors'); >> property_types.sql
echo. >> property_types.sql
echo INSERT INTO property_types (name, data_type, description) >> property_types.sql
echo SELECT 'h_bond_acceptors', 'numeric', 'Number of hydrogen bond acceptors' >> property_types.sql
echo WHERE NOT EXISTS (SELECT 1 FROM property_types WHERE name = 'h_bond_acceptors'); >> property_types.sql
echo. >> property_types.sql
echo -- Return the property type IDs >> property_types.sql
echo SELECT id, name FROM property_types >> property_types.sql
echo WHERE name IN ('logP', 'h_bond_donors', 'h_bond_acceptors'); >> property_types.sql

REM Execute the SQL
supabase db execute --project-ref %PROJECT_ID% < property_types.sql
if %ERRORLEVEL% neq 0 (
    echo Error creating property types. Exiting.
    del property_types.sql
    exit /b 1
)
del property_types.sql

REM Process each reference compound
for %%C in (CHEMBL1098659 CHEMBL1487 CHEMBL262548 CHEMBL388978 CHEMBL500033 CHEMBL6196 CHEMBL66195 CHEMBL6752 CHEMBL967) do (
    echo Checking properties for %%C...
    
    REM Create a temporary SQL file for checking
    echo WITH molecule AS ( > check_props.sql
    echo     SELECT id FROM molecules WHERE chembl_id = '%%C' >> check_props.sql
    echo ), >> check_props.sql
    echo property_types AS ( >> check_props.sql
    echo     SELECT id, name FROM property_types >> check_props.sql
    echo     WHERE name IN ('logP', 'h_bond_donors', 'h_bond_acceptors') >> check_props.sql
    echo ), >> check_props.sql
    echo existing_properties AS ( >> check_props.sql
    echo     SELECT pt.name >> check_props.sql
    echo     FROM molecular_properties mp >> check_props.sql
    echo     JOIN property_types pt ON mp.property_type_id = pt.id >> check_props.sql
    echo     JOIN molecule m ON mp.molecule_id = m.id >> check_props.sql
    echo     WHERE pt.name IN ('logP', 'h_bond_donors', 'h_bond_acceptors') >> check_props.sql
    echo ) >> check_props.sql
    echo SELECT pt.name AS missing_property >> check_props.sql
    echo FROM property_types pt >> check_props.sql
    echo WHERE pt.name NOT IN (SELECT name FROM existing_properties); >> check_props.sql
    
    REM Execute the check SQL and capture output
    supabase db execute --project-ref %PROJECT_ID% < check_props.sql > check_result.txt
    
    REM Check if there are missing properties
    findstr /C:"missing_property" check_result.txt > nul
    if %ERRORLEVEL% equ 0 (
        echo Found missing properties for %%C. Fixing...
        
        REM Create a temporary SQL file for fixing
        echo -- Get molecule ID > fix_props.sql
        echo WITH molecule AS ( >> fix_props.sql
        echo     SELECT id FROM molecules WHERE chembl_id = '%%C' >> fix_props.sql
        echo ), >> fix_props.sql
        echo -- Get property type IDs >> fix_props.sql
        echo property_types AS ( >> fix_props.sql
        echo     SELECT id, name FROM property_types >> fix_props.sql
        echo     WHERE name IN ('logP', 'h_bond_donors', 'h_bond_acceptors') >> fix_props.sql
        echo ), >> fix_props.sql
        echo -- Get existing properties >> fix_props.sql
        echo existing_properties AS ( >> fix_props.sql
        echo     SELECT pt.name >> fix_props.sql
        echo     FROM molecular_properties mp >> fix_props.sql
        echo     JOIN property_types pt ON mp.property_type_id = pt.id >> fix_props.sql
        echo     JOIN molecule m ON mp.molecule_id = m.id >> fix_props.sql
        echo     WHERE pt.name IN ('logP', 'h_bond_donors', 'h_bond_acceptors') >> fix_props.sql
        echo ) >> fix_props.sql
        echo -- Insert missing properties >> fix_props.sql
        echo INSERT INTO molecular_properties (molecule_id, property_type_id, numeric_value, created_at, updated_at) >> fix_props.sql
        echo SELECT >> fix_props.sql
        echo     (SELECT id FROM molecule), >> fix_props.sql
        echo     pt.id, >> fix_props.sql
        echo     CASE >> fix_props.sql
        echo         WHEN pt.name = 'logP' THEN 0.0 >> fix_props.sql
        echo         ELSE 0 >> fix_props.sql
        echo     END, >> fix_props.sql
        echo     NOW(), >> fix_props.sql
        echo     NOW() >> fix_props.sql
        echo FROM property_types pt >> fix_props.sql
        echo WHERE pt.name NOT IN (SELECT name FROM existing_properties); >> fix_props.sql
        
        REM Execute the fix SQL
        supabase db execute --project-ref %PROJECT_ID% < fix_props.sql
        if %ERRORLEVEL% neq 0 (
            echo Error fixing properties for %%C. Exiting.
            del check_props.sql
            del check_result.txt
            del fix_props.sql
            exit /b 1
        )
        del fix_props.sql
        echo Properties fixed for %%C.
    ) else (
        echo No missing properties for %%C.
    )
    
    del check_props.sql
    del check_result.txt
)

echo Verifying all properties are now set...

REM Create a temporary SQL file for final verification
echo WITH reference_molecules AS ( > final_verify.sql
echo     SELECT id, chembl_id >> final_verify.sql
echo     FROM molecules >> final_verify.sql
echo     WHERE chembl_id IN ('CHEMBL1098659', 'CHEMBL1487', 'CHEMBL262548', 'CHEMBL388978', >> final_verify.sql
echo                        'CHEMBL500033', 'CHEMBL6196', 'CHEMBL66195', 'CHEMBL6752', 'CHEMBL967') >> final_verify.sql
echo ), >> final_verify.sql
echo property_counts AS ( >> final_verify.sql
echo     SELECT >> final_verify.sql
echo         rm.chembl_id, >> final_verify.sql
echo         COUNT(DISTINCT pt.name) AS property_count >> final_verify.sql
echo     FROM reference_molecules rm >> final_verify.sql
echo     JOIN molecular_properties mp ON rm.id = mp.molecule_id >> final_verify.sql
echo     JOIN property_types pt ON mp.property_type_id = pt.id >> final_verify.sql
echo     WHERE pt.name IN ('logP', 'h_bond_donors', 'h_bond_acceptors') >> final_verify.sql
echo     GROUP BY rm.chembl_id >> final_verify.sql
echo ) >> final_verify.sql
echo SELECT >> final_verify.sql
echo     rm.chembl_id, >> final_verify.sql
echo     COALESCE(pc.property_count, 0) AS property_count, >> final_verify.sql
echo     CASE >> final_verify.sql
echo         WHEN COALESCE(pc.property_count, 0) = 3 THEN 'Complete' >> final_verify.sql
echo         ELSE 'Incomplete' >> final_verify.sql
echo     END AS status >> final_verify.sql
echo FROM reference_molecules rm >> final_verify.sql
echo LEFT JOIN property_counts pc ON rm.chembl_id = pc.chembl_id >> final_verify.sql
echo ORDER BY rm.chembl_id; >> final_verify.sql

REM Execute the final verification SQL
supabase db execute --project-ref %PROJECT_ID% < final_verify.sql
if %ERRORLEVEL% neq 0 (
    echo Error in final verification. Exiting.
    del final_verify.sql
    exit /b 1
)
del final_verify.sql

echo Reference compounds fix complete!