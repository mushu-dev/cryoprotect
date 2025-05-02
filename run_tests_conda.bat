@echo off
call conda activate cryoprotect
python tests/run_tests.py %*