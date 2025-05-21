#!/bin/bash
# Run the property completion script in the background

# Navigate to the project directory
cd /home/mushu/Projects/CryoProtect

# Run the script in the background
nohup python3 complete_missing_properties_independent.py --batch-size 100 > property_completion.log 2>&1 &

echo "Started property completion in the background. Check property_completion.log for progress."
echo "Process ID: $!"