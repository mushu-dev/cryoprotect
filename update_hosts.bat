@echo off
echo CryoProtect v2 - Update Hosts File
echo This script will update your hosts file with Supabase DNS entries

echo Adding entries to hosts file...
echo 172.64.149.246    tsdlmynydfuypiugmkev.supabase.co >> %windir%\System32\drivers\etc\hosts
echo 2600:1f18:2e13:9d04:d2da:b6ad:4834:140d    db.tsdlmynydfuypiugmkev.supabase.co >> %windir%\System32\drivers\etc\hosts

echo Done! Your hosts file has been updated.
echo Testing connection...
ping -n 1 tsdlmynydfuypiugmkev.supabase.co
ping -n 1 db.tsdlmynydfuypiugmkev.supabase.co

echo.
echo Run 'python supabase_connection_diagnostic.py' to verify the connection
pause
