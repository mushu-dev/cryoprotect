#!/bin/bash
# Run commands in the CryoProtect container with conda environment activated

# Check if container exists
if ! podman container exists CryoProtect; then
    echo "Error: CryoProtect container does not exist."
    echo "Please run create_cryoprotect_conda_container.sh first."
    exit 1
fi

# Check if container is running
if ! podman container inspect CryoProtect --format '{{.State.Running}}' | grep -q "true"; then
    echo "Starting CryoProtect container..."
    podman start CryoProtect
fi

# Execute the command in the container with conda environment activated
if [ $# -eq 0 ]; then
    # If no arguments, start an interactive shell with conda environment activated
    echo "Starting interactive shell in container with conda environment..."
    podman exec -it CryoProtect bash -c "
        # Try to initialize conda if not already done
        conda init bash &>/dev/null || true
        source ~/.bashrc &>/dev/null || true
        # Try to activate with conda or use the direct path to the environment's Python
        echo 'Activating conda environment...'
        conda activate cryoprotect &>/dev/null
        if [ $? -ne 0 ]; then
            echo 'Using direct path to conda environment instead'
            export PATH=/opt/conda/envs/cryoprotect/bin:$PATH
        fi
        bash
    "
else
    # Otherwise, run the specified command with conda environment activated
    echo "Running command in container with conda environment: $@"
    podman exec -it CryoProtect bash -c "
        # Try using direct path to the environment's Python first
        export PATH=/opt/conda/envs/cryoprotect/bin:$PATH
        $*
    "
fi
