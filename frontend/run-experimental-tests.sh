#!/bin/bash

# Run experimental data enhancement feature tests

echo "Running Experimental Data Enhancement tests..."
echo "====================================="

# Run basic user flow and data integrity tests first
echo "Step 1: Running user flow tests..."
npm run test:user-flow

# Run experimental data feature tests
echo "Step 2: Running experimental data feature tests..."
npm run test:experimental-data

# Run data analysis tests
echo "Step 3: Running data analysis tests..."
npm run test:data-analysis

# Generate and show test report
echo "====================================="
echo "Generating test report..."
npm run test:e2e:report

echo "All experimental data enhancement tests completed!"