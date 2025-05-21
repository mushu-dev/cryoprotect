#!/bin/bash

# Try the pipx installation first
if [ -x "/home/mushu/.local/share/pipx/venvs/gh-cli/bin/gh" ]; then
    /home/mushu/.local/share/pipx/venvs/gh-cli/bin/gh "$@"
    exit $?
# Then try the system installation
elif [ -x "/usr/bin/gh" ]; then
    /usr/bin/gh "$@"
    exit $?
else
    echo "ERROR: GitHub CLI not found in expected locations"
    echo "Searched: /home/mushu/.local/share/pipx/venvs/gh-cli/bin/gh, /usr/bin/gh"
    exit 1
fi
