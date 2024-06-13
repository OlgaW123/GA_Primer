#!/bin/bash

# File to log resource usage
LOG_FILE="resource_usage.txt"

# Clear the log file
echo "Time (s), CPU (%), Memory (MB)" > $LOG_FILE

# Function to log resource usage
log_resources() {
    start_time=$(date +%s)
    while true; do
        # Get the current time in seconds
        current_time=$(date +%s)
        elapsed_time=$((current_time - start_time))

        # Get the CPU and Memory usage
        cpu_usage=$(top -b -n1 | grep "Cpu(s)" | awk '{print $2 + $4}')
        memory_usage=$(ps --no-headers -Ao rss | awk '{sum+=$1} END {print sum / 1024}')  # Convert KB to MB

        # Write to the log file
        echo "$elapsed_time, $cpu_usage, $memory_usage" >> $LOG_FILE

        # Sleep for 1 second
        sleep 1
    done
}

# Start logging resources in the background
log_resources &
LOGGER_PID=$!

# Run the Python script
python3 cpu.py

# Stop logging after the script finishes
kill $LOGGER_PID
