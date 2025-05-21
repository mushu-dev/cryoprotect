#!/bin/bash
# Run the property completion script in the background

echo "Starting property completion in the background..."
nohup python complete_missing_properties_fixed.py > property_completion.log 2>&1 &
echo "Background job started with PID $!"
echo "You can monitor progress with: tail -f property_completion.log"