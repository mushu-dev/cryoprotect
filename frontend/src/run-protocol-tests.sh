#!/bin/bash

echo "Running Protocol Version Comparison Tests..."
npm test -- --testPathPattern=tests/protocol-version-comparison.test.js

echo ""
echo "Running Protocol Detail Page Tests..."
npm test -- --testPathPattern=tests/protocol-detail.test.js

echo ""
echo "All tests complete!"