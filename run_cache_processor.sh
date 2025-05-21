#!/bin/bash
# Run the cache invalidation processor as a background service

# Set working directory to the script directory
cd "$(dirname "$0")"

# Create log directory if it doesn't exist
mkdir -p logs

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    echo "Error: Python 3 is required but could not be found"
    exit 1
fi

# Check if Redis is running
if ! command -v redis-cli &> /dev/null || ! redis-cli ping &> /dev/null; then
    echo "Warning: Redis server may not be running"
    echo "Cache invalidation will continue to work, but no caching will occur"
fi

# Check if the processor is already running
if pgrep -f "python3.*cache_processor.py" > /dev/null; then
    echo "Cache processor is already running"
    exit 0
fi

# Run processor in daemon mode
echo "Starting cache invalidation processor..."
python3 -m database.cache_processor --daemon --interval 5 > logs/cache_processor.log 2>&1 &

# Save PID to file
echo $! > cache_processor.pid
echo "Cache processor started with PID $!"

# Set up a stopping script if it doesn't exist
if [ ! -f stop_cache_processor.sh ]; then
    cat > stop_cache_processor.sh << 'EOF'
#!/bin/bash
# Stop the cache invalidation processor

# Set working directory to the script directory
cd "$(dirname "$0")"

# Check if PID file exists
if [ ! -f cache_processor.pid ]; then
    echo "Cache processor is not running (PID file not found)"
    exit 0
fi

# Get PID from file
PID=$(cat cache_processor.pid)

# Check if process is still running
if ! ps -p $PID > /dev/null; then
    echo "Cache processor is not running (PID $PID not found)"
    rm cache_processor.pid
    exit 0
fi

# Kill the process
echo "Stopping cache processor (PID $PID)..."
kill $PID

# Wait for process to terminate
for i in {1..10}; do
    if ! ps -p $PID > /dev/null; then
        echo "Cache processor stopped"
        rm cache_processor.pid
        exit 0
    fi
    sleep 1
done

# Force kill if needed
echo "Cache processor did not terminate gracefully, force killing..."
kill -9 $PID
rm cache_processor.pid
echo "Cache processor stopped"
EOF
    chmod +x stop_cache_processor.sh
    echo "Created stop_cache_processor.sh script"
fi

echo "Cache processor is now running in the background"
echo "To stop it, run ./stop_cache_processor.sh"