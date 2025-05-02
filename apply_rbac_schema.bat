@echo off
echo Applying RBAC schema migration...
python apply_rbac_schema.py
if %ERRORLEVEL% NEQ 0 (
    echo Error applying RBAC schema migration
    exit /b %ERRORLEVEL%
)
echo RBAC schema migration applied successfully