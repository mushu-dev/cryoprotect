@echo off
echo CryoProtect v2 - Create exec_sql Function
echo =======================================
echo.

REM Check if Python is installed
where python >nul 2>nul
if %ERRORLEVEL% neq 0 (
    echo Error: Python is not installed or not in PATH.
    echo Please install Python 3.8 or higher and try again.
    exit /b 1
)

REM Check if .env file exists
if not exist .env (
    echo Warning: .env file not found.
    echo Creating a template .env file. Please edit it with your Supabase credentials.
    echo SUPABASE_URL=your_supabase_url> .env
    echo SUPABASE_KEY=your_supabase_key>> .env
    echo SUPABASE_USER=your_supabase_user_email>> .env
    echo SUPABASE_PASSWORD=your_supabase_user_password>> .env
    echo.
    echo Template .env file created. Please edit it and run this script again.
    exit /b 1
)

REM Check if required packages are installed
echo Checking required packages...
pip show python-dotenv supabase >nul 2>nul
if %ERRORLEVEL% neq 0 (
    echo Installing required packages...
    pip install python-dotenv supabase
)

REM Create a temporary Python script to execute the SQL
echo import os> create_function_temp.py
echo from dotenv import load_dotenv>> create_function_temp.py
echo from supabase import create_client>> create_function_temp.py
echo.>> create_function_temp.py
echo load_dotenv()>> create_function_temp.py
echo.>> create_function_temp.py
echo SUPABASE_URL = os.getenv("SUPABASE_URL")>> create_function_temp.py
echo SUPABASE_KEY = os.getenv("SUPABASE_KEY")>> create_function_temp.py
echo.>> create_function_temp.py
echo if not SUPABASE_URL or not SUPABASE_KEY:>> create_function_temp.py
echo     raise ValueError("SUPABASE_URL and SUPABASE_KEY must be set in .env file")>> create_function_temp.py
echo.>> create_function_temp.py
echo supabase = create_client(SUPABASE_URL, SUPABASE_KEY)>> create_function_temp.py
echo.>> create_function_temp.py
echo with open("create_exec_sql_function.sql", "r") as f:>> create_function_temp.py
echo     sql = f.read()>> create_function_temp.py
echo.>> create_function_temp.py
echo print("Executing SQL to create exec_sql function...")>> create_function_temp.py
echo response = supabase.rpc("exec_sql", {"query": sql}).execute()>> create_function_temp.py
echo.>> create_function_temp.py
echo print("Function created successfully!")>> create_function_temp.py
echo print("You can now use the exec_sql function in your verification scripts.")>> create_function_temp.py

REM Run the temporary Python script
echo.
echo Creating exec_sql function in Supabase...
python create_function_temp.py

REM Check if the script ran successfully
if %ERRORLEVEL% neq 0 (
    echo.
    echo Error: Failed to create exec_sql function.
    echo Please check your Supabase credentials and try again.
    echo You can also create the function manually using the SQL in create_exec_sql_function.sql.
    del create_function_temp.py
    exit /b 1
)

REM Clean up
del create_function_temp.py

echo.
echo exec_sql function created successfully!
echo You can now run the verification script.
echo.