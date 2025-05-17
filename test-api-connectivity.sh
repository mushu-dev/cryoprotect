#!/bin/bash
# Simple script to test API connectivity to the Heroku backend

echo "Testing API connectivity to Heroku backend..."
echo "Endpoint: https://cryoprotect-8030e4025428.herokuapp.com/api/v1/health/connectivity"
echo ""

curl -s https://cryoprotect-8030e4025428.herokuapp.com/api/v1/health/connectivity | json_pp

echo ""
echo "Testing API connectivity using simplified endpoint..."
echo "Endpoint: https://cryoprotect-8030e4025428.herokuapp.com/v1/health/connectivity"
echo ""

curl -s https://cryoprotect-8030e4025428.herokuapp.com/v1/health/connectivity | json_pp

exit 0