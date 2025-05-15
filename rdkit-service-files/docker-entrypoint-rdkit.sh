#!/bin/bash
set -e

# Activate conda environment
source /opt/conda/etc/profile.d/conda.sh
conda activate cryoprotect

# Execute the provided command
exec "$@"