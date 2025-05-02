@echo off
:: Setup script for automated backup scheduling on Windows using Task Scheduler

:: Get the absolute path to the project directory
set "PROJECT_DIR=%~dp0.."
set "BACKUP_SCRIPT=%PROJECT_DIR%\backup\backup_manager.py"
set "CONFIG_FILE=%PROJECT_DIR%\backup\backup_config.json"
set "LOG_FILE=%PROJECT_DIR%\logs\backup.log"

:: Ensure the logs directory exists
if not exist "%PROJECT_DIR%\logs" mkdir "%PROJECT_DIR%\logs"

:: Check if config file exists
if not exist "%CONFIG_FILE%" (
    echo Config file not found. Creating from template...
    if exist "%PROJECT_DIR%\backup\backup_config.json.template" (
        copy "%PROJECT_DIR%\backup\backup_config.json.template" "%CONFIG_FILE%"
        echo Created config file from template. Please edit %CONFIG_FILE% to customize your backup settings.
    ) else (
        echo Error: Template config file not found.
        exit /b 1
    )
)

:: Parse the config file to get the schedule
:: This is a simplified version - in a real implementation, you would use a JSON parser
for /f "tokens=2 delims=:," %%a in ('findstr /C:"\"daily\":" "%CONFIG_FILE%"') do (
    set "DAILY_TIME=%%a"
    set "DAILY_TIME=!DAILY_TIME:"=!"
    set "DAILY_TIME=!DAILY_TIME: =!"
)

for /f "tokens=2 delims=:," %%a in ('findstr /C:"\"weekly\":" "%CONFIG_FILE%"') do (
    set "WEEKLY_DAY=%%a"
    set "WEEKLY_DAY=!WEEKLY_DAY:"=!"
    set "WEEKLY_DAY=!WEEKLY_DAY: =!"
)

for /f "tokens=2 delims=:," %%a in ('findstr /C:"\"monthly\":" "%CONFIG_FILE%"') do (
    set "MONTHLY_DAY=%%a"
    set "MONTHLY_DAY=!MONTHLY_DAY: =!"
)

:: Parse the daily time
for /f "tokens=1,2 delims=:" %%a in ("%DAILY_TIME%") do (
    set "HOUR=%%a"
    set "MINUTE=%%b"
)

:: Create the scheduled tasks
echo Creating scheduled tasks...

:: Daily backup task
schtasks /create /tn "CryoProtect v2 Daily Backup" /tr "cmd /c cd /d %PROJECT_DIR% && python -m backup.backup_manager backup --config %CONFIG_FILE% --type daily >> %LOG_FILE% 2>&1" /sc daily /st %DAILY_TIME% /f

:: Weekly backup task
schtasks /create /tn "CryoProtect v2 Weekly Backup" /tr "cmd /c cd /d %PROJECT_DIR% && python -m backup.backup_manager backup --config %CONFIG_FILE% --type weekly >> %LOG_FILE% 2>&1" /sc weekly /d %WEEKLY_DAY% /st %DAILY_TIME% /f

:: Monthly backup task
schtasks /create /tn "CryoProtect v2 Monthly Backup" /tr "cmd /c cd /d %PROJECT_DIR% && python -m backup.backup_manager backup --config %CONFIG_FILE% --type monthly >> %LOG_FILE% 2>&1" /sc monthly /d %MONTHLY_DAY% /st %DAILY_TIME% /f

echo Backup scheduler setup complete.
echo Daily backup: Every day at %DAILY_TIME%
echo Weekly backup: Every %WEEKLY_DAY% at %DAILY_TIME%
echo Monthly backup: Every month on day %MONTHLY_DAY% at %DAILY_TIME%
echo Logs will be written to: %LOG_FILE%

pause