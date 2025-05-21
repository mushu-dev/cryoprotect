#!/bin/bash

# Set the environment for testing
export NODE_ENV=test

# Run the basic API tests
npx jest --no-config tests/api/basic-api-test.js

exit $?