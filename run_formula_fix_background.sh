#!/bin/bash
# Run the formula fixing script in the background

# Navigate to the project directory
cd /home/mushu/Projects/CryoProtect

# Run the script in the background
nohup python3 fix_missing_formulas.py --batch-size 100 > formula_fix.log 2>&1 &

echo "Started formula fixing in the background. Check formula_fix.log for progress."
echo "Process ID: $!"