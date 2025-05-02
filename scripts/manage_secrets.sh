#!/bin/bash
# CryoProtect v2 Secret Management Script
# This script helps create, update, and manage Docker secrets for CryoProtect v2

set -e

# Default values
SECRET_PREFIX="cryoprotect"
SECRET_DIR="./.secrets"
MODE="help"
SECRET_NAME=""
SECRET_VALUE=""
ENVIRONMENT=""
DEPLOYMENT_COLOR=""
LIST_SECRETS=false
ROTATE_SECRETS=false
EXTERNAL_PROVIDER=""
DOCKER_COMPOSE_FILE="docker-compose.yml"

# ANSI color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to display help
show_help() {
  echo -e "${BLUE}CryoProtect v2 Secret Management Script${NC}"
  echo ""
  echo "Usage: $0 [options]"
  echo ""
  echo "Options:"
  echo "  -m, --mode MODE           Operation mode: create, update, delete, dev, swarm, k8s"
  echo "  -n, --name NAME           Secret name (e.g., SUPABASE_URL)"
  echo "  -v, --value VALUE         Secret value (or - to read from stdin)"
  echo "  -e, --environment ENV     Environment prefix: staging, production"
  echo "  -c, --color COLOR         Deployment color: blue, green"
  echo "  -p, --prefix PREFIX       Secret name prefix (default: cryoprotect)"
  echo "  -d, --dir DIRECTORY       Directory for dev secrets (default: ./.secrets)"
  echo "  -f, --file FILE           Docker compose file (default: docker-compose.yml)"
  echo "  -l, --list                List all secrets"
  echo "  -r, --rotate              Rotate secrets (update timestamp)"
  echo "  -x, --external PROVIDER   Configure external provider: aws, vault, azure, gcp"
  echo "  -h, --help                Show this help message"
  echo ""
  echo "Examples:"
  echo "  # Create a Docker Swarm secret"
  echo "  $0 --mode swarm --name SUPABASE_URL --value \"https://example.supabase.co\""
  echo ""
  echo "  # Create a secret for staging environment"
  echo "  $0 --mode swarm --name SUPABASE_URL --value \"https://staging.supabase.co\" --environment staging"
  echo ""
  echo "  # Create a secret for blue deployment"
  echo "  $0 --mode swarm --name SUPABASE_URL --value \"https://blue.supabase.co\" --color blue"
  echo ""
  echo "  # Set up development secrets"
  echo "  $0 --mode dev --name SUPABASE_URL --value \"https://dev.supabase.co\""
  echo ""
  echo "  # List all secrets"
  echo "  $0 --list"
  echo ""
  echo "  # Rotate secrets (update timestamp)"
  echo "  $0 --rotate"
  echo ""
  echo "  # Configure AWS Secrets Manager"
  echo "  $0 --external aws --name AWS_ACCESS_KEY_ID --value \"your-access-key\""
  echo ""
}

# Function to validate input
validate_input() {
  if [[ "$MODE" != "help" && "$MODE" != "list" && "$MODE" != "rotate" && -z "$SECRET_NAME" ]]; then
    echo -e "${RED}Error: Secret name is required${NC}"
    exit 1
  fi
  
  if [[ "$MODE" == "create" || "$MODE" == "update" || "$MODE" == "swarm" || "$MODE" == "dev" || "$MODE" == "k8s" ]]; then
    if [[ -z "$SECRET_VALUE" && "$SECRET_VALUE" != "-" ]]; then
      echo -e "${RED}Error: Secret value is required${NC}"
      exit 1
    fi
  fi
}

# Function to create a full secret name with prefix, environment, and color
get_full_secret_name() {
  local base_name="$1"
  local full_name="${SECRET_PREFIX}_${base_name}"
  
  if [[ -n "$ENVIRONMENT" ]]; then
    full_name="${SECRET_PREFIX}_${ENVIRONMENT}_${base_name}"
  fi
  
  if [[ -n "$DEPLOYMENT_COLOR" ]]; then
    full_name="${SECRET_PREFIX}_${DEPLOYMENT_COLOR}_${base_name}"
  fi
  
  echo "$full_name"
}

# Function to create a Docker Swarm secret
create_swarm_secret() {
  local name="$1"
  local value="$2"
  local full_name=$(get_full_secret_name "$name")
  
  # Check if secret already exists
  if docker secret inspect "$full_name" &>/dev/null; then
    echo -e "${YELLOW}Secret $full_name already exists. Use update mode to change it.${NC}"
    return 1
  fi
  
  # Create the secret
  if [[ "$value" == "-" ]]; then
    echo -e "${BLUE}Enter secret value for $name (input will be hidden):${NC}"
    read -s value
    echo "$value" | docker secret create "$full_name" -
  else
    echo "$value" | docker secret create "$full_name" -
  fi
  
  echo -e "${GREEN}Created Docker Swarm secret: $full_name${NC}"
}

# Function to update a Docker Swarm secret
update_swarm_secret() {
  local name="$1"
  local value="$2"
  local full_name=$(get_full_secret_name "$name")
  
  # Check if secret exists
  if ! docker secret inspect "$full_name" &>/dev/null; then
    echo -e "${YELLOW}Secret $full_name does not exist. Creating it instead.${NC}"
    create_swarm_secret "$name" "$value"
    return
  fi
  
  # Docker Swarm secrets cannot be updated directly, so we need to delete and recreate
  docker secret rm "$full_name"
  
  # Create the secret
  if [[ "$value" == "-" ]]; then
    echo -e "${BLUE}Enter new secret value for $name (input will be hidden):${NC}"
    read -s value
    echo "$value" | docker secret create "$full_name" -
  else
    echo "$value" | docker secret create "$full_name" -
  fi
  
  echo -e "${GREEN}Updated Docker Swarm secret: $full_name${NC}"
}

# Function to delete a Docker Swarm secret
delete_swarm_secret() {
  local name="$1"
  local full_name=$(get_full_secret_name "$name")
  
  # Check if secret exists
  if ! docker secret inspect "$full_name" &>/dev/null; then
    echo -e "${YELLOW}Secret $full_name does not exist.${NC}"
    return 1
  fi
  
  # Delete the secret
  docker secret rm "$full_name"
  
  echo -e "${GREEN}Deleted Docker Swarm secret: $full_name${NC}"
}

# Function to create a development secret
create_dev_secret() {
  local name="$1"
  local value="$2"
  
  # Create the secrets directory if it doesn't exist
  mkdir -p "$SECRET_DIR"
  
  # Add to .gitignore if it's not already there
  if [[ -f ".gitignore" ]]; then
    if ! grep -q "^$SECRET_DIR" .gitignore; then
      echo "$SECRET_DIR" >> .gitignore
      echo -e "${BLUE}Added $SECRET_DIR to .gitignore${NC}"
    fi
  else
    echo "$SECRET_DIR" > .gitignore
    echo -e "${BLUE}Created .gitignore with $SECRET_DIR${NC}"
  fi
  
  # Create the secret file
  local secret_file="$SECRET_DIR/$name"
  
  if [[ -f "$secret_file" ]]; then
    echo -e "${YELLOW}Secret file $secret_file already exists. Overwriting.${NC}"
  fi
  
  if [[ "$value" == "-" ]]; then
    echo -e "${BLUE}Enter secret value for $name (input will be hidden):${NC}"
    read -s value
    echo -n "$value" > "$secret_file"
  else
    echo -n "$value" > "$secret_file"
  fi
  
  # Set proper permissions
  chmod 600 "$secret_file"
  
  echo -e "${GREEN}Created development secret: $secret_file${NC}"
  
  # Update docker-compose.yml if it exists
  if [[ -f "$DOCKER_COMPOSE_FILE" ]]; then
    echo -e "${BLUE}To use this secret in development, add the following to your $DOCKER_COMPOSE_FILE:${NC}"
    echo ""
    echo "secrets:"
    echo "  $name:"
    echo "    file: $secret_file"
    echo ""
    echo "services:"
    echo "  cryoprotect-dev:"
    echo "    secrets:"
    echo "      - source: $name"
    echo "        target: $name"
    echo "        mode: 0400"
  fi
}

# Function to create a Kubernetes secret
create_k8s_secret() {
  local name="$1"
  local value="$2"
  local namespace="${K8S_NAMESPACE:-default}"
  
  # Check if kubectl is available
  if ! command -v kubectl &>/dev/null; then
    echo -e "${RED}Error: kubectl is not installed or not in PATH${NC}"
    exit 1
  fi
  
  # Check if secret already exists
  if kubectl get secret "$SECRET_PREFIX" -n "$namespace" &>/dev/null; then
    # Secret exists, update it
    if [[ "$value" == "-" ]]; then
      echo -e "${BLUE}Enter secret value for $name (input will be hidden):${NC}"
      read -s value
    fi
    
    kubectl patch secret "$SECRET_PREFIX" -n "$namespace" -p "{\"data\":{\"$name\":\"$(echo -n "$value" | base64 -w 0)\"}}"
    echo -e "${GREEN}Updated Kubernetes secret: $name in $SECRET_PREFIX${NC}"
  else
    # Secret doesn't exist, create it
    if [[ "$value" == "-" ]]; then
      echo -e "${BLUE}Enter secret value for $name (input will be hidden):${NC}"
      read -s value
    fi
    
    kubectl create secret generic "$SECRET_PREFIX" -n "$namespace" --from-literal="$name=$value"
    echo -e "${GREEN}Created Kubernetes secret: $name in $SECRET_PREFIX${NC}"
  fi
  
  echo -e "${BLUE}To use this secret in Kubernetes, set the following environment variables:${NC}"
  echo "K8S_SECRET_NAMESPACE=$namespace"
  echo "K8S_SECRET_NAME=$SECRET_PREFIX"
}

# Function to list secrets
list_secrets() {
  echo -e "${BLUE}Listing secrets:${NC}"
  
  # Check for Docker Swarm secrets
  if command -v docker &>/dev/null; then
    echo -e "${BLUE}Docker Swarm secrets:${NC}"
    docker secret ls | grep "$SECRET_PREFIX" || echo "No Docker Swarm secrets found"
    echo ""
  fi
  
  # Check for development secrets
  if [[ -d "$SECRET_DIR" ]]; then
    echo -e "${BLUE}Development secrets:${NC}"
    ls -la "$SECRET_DIR" || echo "No development secrets found"
    echo ""
  fi
  
  # Check for Kubernetes secrets
  if command -v kubectl &>/dev/null; then
    echo -e "${BLUE}Kubernetes secrets:${NC}"
    kubectl get secret "$SECRET_PREFIX" --all-namespaces -o json 2>/dev/null | \
      jq -r '.items[] | "Namespace: \(.metadata.namespace), Name: \(.metadata.name), Keys: \((.data | keys) | join(", "))"' || \
      echo "No Kubernetes secrets found or jq not installed"
    echo ""
  fi
}

# Function to rotate secrets
rotate_secrets() {
  local timestamp=$(date +%s)
  
  echo -e "${BLUE}Rotating secrets (updating timestamp to $timestamp)${NC}"
  
  # Update Docker Swarm secret
  if command -v docker &>/dev/null; then
    local rotation_secret="${SECRET_PREFIX}_SECRET_ROTATION_TIMESTAMP"
    
    # Check if secret already exists
    if docker secret inspect "$rotation_secret" &>/dev/null; then
      docker secret rm "$rotation_secret"
    fi
    
    # Create the secret
    echo "$timestamp" | docker secret create "$rotation_secret" -
    echo -e "${GREEN}Updated Docker Swarm secret rotation timestamp: $rotation_secret${NC}"
  fi
  
  # Update development secret
  if [[ -d "$SECRET_DIR" ]]; then
    local rotation_file="$SECRET_DIR/SECRET_ROTATION_TIMESTAMP"
    echo "$timestamp" > "$rotation_file"
    chmod 600 "$rotation_file"
    echo -e "${GREEN}Updated development secret rotation timestamp: $rotation_file${NC}"
  fi
  
  # Update Kubernetes secret
  if command -v kubectl &>/dev/null; then
    local namespace="${K8S_NAMESPACE:-default}"
    
    if kubectl get secret "$SECRET_PREFIX" -n "$namespace" &>/dev/null; then
      kubectl patch secret "$SECRET_PREFIX" -n "$namespace" -p "{\"data\":{\"SECRET_ROTATION_TIMESTAMP\":\"$(echo -n "$timestamp" | base64 -w 0)\"}}"
      echo -e "${GREEN}Updated Kubernetes secret rotation timestamp in $SECRET_PREFIX${NC}"
    else
      kubectl create secret generic "$SECRET_PREFIX" -n "$namespace" --from-literal="SECRET_ROTATION_TIMESTAMP=$timestamp"
      echo -e "${GREEN}Created Kubernetes secret rotation timestamp in $SECRET_PREFIX${NC}"
    fi
  fi
}

# Function to configure external provider
configure_external_provider() {
  local provider="$1"
  local name="$2"
  local value="$3"
  
  case "$provider" in
    aws)
      echo -e "${BLUE}Configuring AWS Secrets Manager integration${NC}"
      if [[ -z "$name" || -z "$value" ]]; then
        echo -e "${YELLOW}For AWS integration, you need to create the following secrets:${NC}"
        echo "  - AWS_ACCESS_KEY_ID"
        echo "  - AWS_SECRET_ACCESS_KEY"
        echo "  - AWS_REGION (optional)"
        echo ""
        echo -e "${YELLOW}And set the following environment variable:${NC}"
        echo "  - AWS_SECRET_PREFIX (e.g., /cryoprotect/)"
        return
      fi
      
      # Create the secret based on the current mode
      if [[ "$MODE" == "swarm" ]]; then
        create_swarm_secret "$name" "$value"
      elif [[ "$MODE" == "dev" ]]; then
        create_dev_secret "$name" "$value"
      elif [[ "$MODE" == "k8s" ]]; then
        create_k8s_secret "$name" "$value"
      else
        echo -e "${RED}Error: Please specify a mode (swarm, dev, k8s) when configuring external provider${NC}"
        exit 1
      fi
      
      echo -e "${BLUE}To use AWS Secrets Manager, set EXTERNAL_SECRET_PROVIDER=aws in your environment${NC}"
      ;;
      
    vault)
      echo -e "${BLUE}Configuring HashiCorp Vault integration${NC}"
      if [[ -z "$name" || -z "$value" ]]; then
        echo -e "${YELLOW}For Vault integration, you need to create the following secrets:${NC}"
        echo "  - VAULT_TOKEN"
        echo "  - VAULT_ADDR (optional)"
        echo ""
        echo -e "${YELLOW}And set the following environment variable:${NC}"
        echo "  - VAULT_SECRET_PATH (e.g., secret/cryoprotect)"
        return
      fi
      
      # Create the secret based on the current mode
      if [[ "$MODE" == "swarm" ]]; then
        create_swarm_secret "$name" "$value"
      elif [[ "$MODE" == "dev" ]]; then
        create_dev_secret "$name" "$value"
      elif [[ "$MODE" == "k8s" ]]; then
        create_k8s_secret "$name" "$value"
      else
        echo -e "${RED}Error: Please specify a mode (swarm, dev, k8s) when configuring external provider${NC}"
        exit 1
      fi
      
      echo -e "${BLUE}To use HashiCorp Vault, set EXTERNAL_SECRET_PROVIDER=vault in your environment${NC}"
      ;;
      
    azure)
      echo -e "${BLUE}Configuring Azure Key Vault integration${NC}"
      if [[ -z "$name" || -z "$value" ]]; then
        echo -e "${YELLOW}For Azure integration, you need to create the following secrets:${NC}"
        echo "  - AZURE_CLIENT_ID"
        echo "  - AZURE_CLIENT_SECRET"
        echo "  - AZURE_TENANT_ID"
        echo ""
        echo -e "${YELLOW}And set the following environment variable:${NC}"
        echo "  - AZURE_KEYVAULT_NAME"
        return
      fi
      
      # Create the secret based on the current mode
      if [[ "$MODE" == "swarm" ]]; then
        create_swarm_secret "$name" "$value"
      elif [[ "$MODE" == "dev" ]]; then
        create_dev_secret "$name" "$value"
      elif [[ "$MODE" == "k8s" ]]; then
        create_k8s_secret "$name" "$value"
      else
        echo -e "${RED}Error: Please specify a mode (swarm, dev, k8s) when configuring external provider${NC}"
        exit 1
      fi
      
      echo -e "${BLUE}To use Azure Key Vault, set EXTERNAL_SECRET_PROVIDER=azure in your environment${NC}"
      ;;
      
    gcp)
      echo -e "${BLUE}Configuring Google Cloud Secret Manager integration${NC}"
      if [[ -z "$name" || -z "$value" ]]; then
        echo -e "${YELLOW}For GCP integration, you need to create the following secret:${NC}"
        echo "  - GOOGLE_APPLICATION_CREDENTIALS (service account key file content)"
        echo ""
        echo -e "${YELLOW}And set the following environment variables:${NC}"
        echo "  - GCP_PROJECT_ID"
        echo "  - GCP_SECRET_PREFIX"
        return
      fi
      
      # Create the secret based on the current mode
      if [[ "$MODE" == "swarm" ]]; then
        create_swarm_secret "$name" "$value"
      elif [[ "$MODE" == "dev" ]]; then
        create_dev_secret "$name" "$value"
      elif [[ "$MODE" == "k8s" ]]; then
        create_k8s_secret "$name" "$value"
      else
        echo -e "${RED}Error: Please specify a mode (swarm, dev, k8s) when configuring external provider${NC}"
        exit 1
      fi
      
      echo -e "${BLUE}To use Google Cloud Secret Manager, set EXTERNAL_SECRET_PROVIDER=gcp in your environment${NC}"
      ;;
      
    *)
      echo -e "${RED}Error: Unknown external provider: $provider${NC}"
      echo -e "${YELLOW}Supported providers: aws, vault, azure, gcp${NC}"
      exit 1
      ;;
  esac
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    -m|--mode)
      MODE="$2"
      shift 2
      ;;
    -n|--name)
      SECRET_NAME="$2"
      shift 2
      ;;
    -v|--value)
      SECRET_VALUE="$2"
      shift 2
      ;;
    -e|--environment)
      ENVIRONMENT="$2"
      shift 2
      ;;
    -c|--color)
      DEPLOYMENT_COLOR="$2"
      shift 2
      ;;
    -p|--prefix)
      SECRET_PREFIX="$2"
      shift 2
      ;;
    -d|--dir)
      SECRET_DIR="$2"
      shift 2
      ;;
    -f|--file)
      DOCKER_COMPOSE_FILE="$2"
      shift 2
      ;;
    -l|--list)
      LIST_SECRETS=true
      shift
      ;;
    -r|--rotate)
      ROTATE_SECRETS=true
      shift
      ;;
    -x|--external)
      EXTERNAL_PROVIDER="$2"
      shift 2
      ;;
    -h|--help)
      show_help
      exit 0
      ;;
    *)
      echo -e "${RED}Error: Unknown option: $1${NC}"
      show_help
      exit 1
      ;;
  esac
done

# Handle list mode
if [[ "$LIST_SECRETS" == true ]]; then
  list_secrets
  exit 0
fi

# Handle rotate mode
if [[ "$ROTATE_SECRETS" == true ]]; then
  rotate_secrets
  exit 0
fi

# Handle external provider configuration
if [[ -n "$EXTERNAL_PROVIDER" ]]; then
  configure_external_provider "$EXTERNAL_PROVIDER" "$SECRET_NAME" "$SECRET_VALUE"
  exit 0
fi

# Handle help mode
if [[ "$MODE" == "help" ]]; then
  show_help
  exit 0
fi

# Validate input
validate_input

# Execute the requested operation
case "$MODE" in
  create|swarm)
    create_swarm_secret "$SECRET_NAME" "$SECRET_VALUE"
    ;;
  update)
    update_swarm_secret "$SECRET_NAME" "$SECRET_VALUE"
    ;;
  delete)
    delete_swarm_secret "$SECRET_NAME"
    ;;
  dev)
    create_dev_secret "$SECRET_NAME" "$SECRET_VALUE"
    ;;
  k8s)
    create_k8s_secret "$SECRET_NAME" "$SECRET_VALUE"
    ;;
  *)
    echo -e "${RED}Error: Unknown mode: $MODE${NC}"
    show_help
    exit 1
    ;;
esac