#!/bin/bash
# Simple script to test PostgreSQL connectivity

echo "Testing PostgreSQL connectivity..."
export PGPASSWORD=postgres
psql -h localhost -p 5433 -U postgres -d postgres -c "\l" || echo "Connection failed"