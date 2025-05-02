#!/bin/bash
# generate_sbom.sh - Generate Software Bill of Materials (SBOM) for Docker images
# 
# This script generates an SBOM for Docker images using Trivy
# and can store it in various formats and locations.
#
# Usage: ./generate_sbom.sh [options] IMAGE_NAME
#
# Options:
#   --format FORMAT       Output format (cyclonedx, spdx, json) [default: cyclonedx]
#   --output FILE         Output file name [default: sbom.{format}]
#   --store-dir DIR       Directory to store SBOM files [default: ./sbom]
#   --version VERSION     Version tag for the SBOM file [default: latest]
#   --help                Show this help message and exit

set -e

# Default values
FORMAT="cyclonedx"
OUTPUT=""
STORE_DIR="./sbom"
VERSION="latest"

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --format)
      FORMAT="$2"
      shift 2
      ;;
    --output)
      OUTPUT="$2"
      shift 2
      ;;
    --store-dir)
      STORE_DIR="$2"
      shift 2
      ;;
    --version)
      VERSION="$2"
      shift 2
      ;;
    --help)
      echo "Usage: ./generate_sbom.sh [options] IMAGE_NAME"
      echo ""
      echo "Options:"
      echo "  --format FORMAT       Output format (cyclonedx, spdx, json) [default: cyclonedx]"
      echo "  --output FILE         Output file name [default: sbom.{format}]"
      echo "  --store-dir DIR       Directory to store SBOM files [default: ./sbom]"
      echo "  --version VERSION     Version tag for the SBOM file [default: latest]"
      echo "  --help                Show this help message and exit"
      exit 0
      ;;
    *)
      IMAGE_NAME="$1"
      shift
      ;;
  esac
done

# Check if image name is provided
if [ -z "$IMAGE_NAME" ]; then
  echo "Error: IMAGE_NAME is required"
  echo "Run './generate_sbom.sh --help' for usage information"
  exit 1
fi

# Set default output file if not provided
if [ -z "$OUTPUT" ]; then
  case "$FORMAT" in
    cyclonedx)
      OUTPUT="sbom.json"
      ;;
    spdx)
      OUTPUT="sbom.spdx.json"
      ;;
    json)
      OUTPUT="sbom.json"
      ;;
    *)
      OUTPUT="sbom.json"
      ;;
  esac
fi

echo "Generating SBOM for Docker image: $IMAGE_NAME"
echo "Output format: $FORMAT"
echo "Output file: $OUTPUT"
echo "Version: $VERSION"

# Check if Trivy is installed
if ! command -v trivy &> /dev/null; then
  echo "Error: Trivy is not installed"
  echo "Please install Trivy: https://aquasecurity.github.io/trivy/latest/getting-started/installation/"
  exit 1
fi

# Create storage directory if it doesn't exist
mkdir -p "$STORE_DIR"

# Generate SBOM
trivy image --format "$FORMAT" --output "$OUTPUT" "$IMAGE_NAME"

# Create a versioned copy in the storage directory
TIMESTAMP=$(date -u +"%Y%m%d%H%M%S")
FILENAME=$(basename "$OUTPUT")
EXTENSION="${FILENAME##*.}"
VERSIONED_FILENAME="${IMAGE_NAME//\//_}_${VERSION}_${TIMESTAMP}.${EXTENSION}"
VERSIONED_FILENAME=$(echo "$VERSIONED_FILENAME" | tr ':' '_')

cp "$OUTPUT" "$STORE_DIR/$VERSIONED_FILENAME"

# Create a latest symlink or copy
LATEST_FILENAME="${IMAGE_NAME//\//_}_latest.${EXTENSION}"
LATEST_FILENAME=$(echo "$LATEST_FILENAME" | tr ':' '_')

if [ -f "$STORE_DIR/$LATEST_FILENAME" ]; then
  rm "$STORE_DIR/$LATEST_FILENAME"
fi

cp "$OUTPUT" "$STORE_DIR/$LATEST_FILENAME"

echo "SBOM generated successfully."
echo "Stored in: $STORE_DIR/$VERSIONED_FILENAME"
echo "Latest version: $STORE_DIR/$LATEST_FILENAME"

# Generate metadata file
METADATA_FILE="$STORE_DIR/${VERSIONED_FILENAME%.${EXTENSION}}.metadata.json"
cat > "$METADATA_FILE" << EOF
{
  "image": "$IMAGE_NAME",
  "version": "$VERSION",
  "timestamp": "$(date -u +"%Y-%m-%dT%H:%M:%SZ")",
  "format": "$FORMAT",
  "sbom_file": "$VERSIONED_FILENAME"
}
EOF

echo "Metadata stored in: $METADATA_FILE"