#!/bin/bash
source $(conda info --base)/etc/profile.d/conda.sh
conda activate cryoprotect
python tests/run_tests.py "$@"