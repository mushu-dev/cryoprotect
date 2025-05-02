@echo off
set DATE=%date:~-4%-%date:~4,2%-%date:~7,2%
echo Creating backup of memory bank...
if not exist ".roo\memory\archive" mkdir .roo\memory\archive
copy .roo\memory\activeContext.md .roo\memory\archive\activeContext-%DATE%.md
echo # Active Development Context > .roo\memory\activeContext.md
echo. >> .roo\memory\activeContext.md
echo ## Current Task >> .roo\memory\activeContext.md
echo Task continues from previous context >> .roo\memory\activeContext.md
echo. >> .roo\memory\activeContext.md
echo ## Key Files >> .roo\memory\activeContext.md
echo - [List of relevant files with line numbers] >> .roo\memory\activeContext.md
echo. >> .roo\memory\activeContext.md
echo ## Implementation Notes >> .roo\memory\activeContext.md
echo [Important implementation details] >> .roo\memory\activeContext.md
echo Memory bank rotation complete!
