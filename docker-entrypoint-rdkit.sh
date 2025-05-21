#!/bin/bash
set -e

# Activate conda environment
source /opt/conda/etc/profile.d/conda.sh
conda activate cryoprotect

# Run the specified command
exec "$@"