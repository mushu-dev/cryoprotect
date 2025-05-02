@echo off
echo CryoProtect v2 - Supabase Connection Fix Utility
echo =============================================
echo.
echo This utility will:
echo 1. Run diagnostics to identify Supabase connection issues
echo 2. Attempt to fix DNS and connection problems
echo.
echo NOTE: For some fixes, you may need to run this batch file as Administrator.
echo.
pause

echo.
echo Running diagnostic...
python supabase_connection_diagnostic.py
echo.

set /p choice=Do you want to fix DNS and connection issues? (y/n): 
if /i "%choice%"=="y" (
    echo.
    echo Attempting to fix DNS and connection issues...
    python fix_supabase_dns.py
) else (
    echo.
    echo Fix skipped.
)

echo.
echo Connection diagnosis and fix complete.
echo If issues persist, please contact support or check README_Supabase.md for more information.
echo.
pause