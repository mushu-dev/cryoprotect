#!/bin/bash

# Try all possible locations
if [ -x /home/mushu/.local/bin/gh ]; then
    /home/mushu/.local/bin/gh "$@"
elif [ -x /usr/bin/gh ]; then
    /usr/bin/gh "$@"
else
    echo "GitHub CLI not found in any location"
    exit 1
fi
