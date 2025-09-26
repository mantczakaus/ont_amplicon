#!/usr/bin/bash

# Expand relative paths in CSV input to absolute paths

FILE=$1

if [ -z "$FILE" ]; then
  echo "Usage: $0 <path to CSV file>"
  exit 1
fi

sed -i 's/uploads\//\/mnt\/data\/user-data\/uploads\//' "$FILE"
