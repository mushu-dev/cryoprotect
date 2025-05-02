@echo off
echo Running CryoProtect Supabase Database Backup...
python create_database_backup.py %*
echo Backup completed.
pause