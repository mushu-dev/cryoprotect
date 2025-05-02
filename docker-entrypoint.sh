#!/bin/bash
set -eo pipefail

# Function to log messages with timestamp and log level
log_message() {
  local level="${2:-INFO}"
  echo "[$(date -Iseconds)] [${level}] $1"
}

# Function to log debug messages (only if DEBUG=1)
log_debug() {
  if [ "${DEBUG:-0}" = "1" ]; then
    log_message "$1" "DEBUG"
  fi
}

# Function to log error messages
log_error() {
  log_message "$1" "ERROR"
}

# Function to log warning messages
log_warning() {
  log_message "$1" "WARNING"
}

# Start timing for performance measurement
START_TIME=$(date +%s.%N)

# Function to load secrets into environment variables
load_secrets() {
  log_message "Loading secrets..."
  
  # Define required secrets - application will fail if these are missing
  REQUIRED_SECRETS=(
    "SUPABASE_URL"
    "SUPABASE_KEY"
    "SECRET_KEY"
  )
  
  # Define optional secrets - application will continue if these are missing
  OPTIONAL_SECRETS=(
    "SUPABASE_SERVICE_KEY"
    "SUPABASE_SERVICE_ROLE_KEY"
    "SUPABASE_DB_PASSWORD"
    "DATABASE_URL"
    "REDIS_URL"
    "NOTIFY_SMTP_PASSWORD"
    "MAIL_PASSWORD"
    "API_KEY"
    "JWT_SECRET"
    "ENCRYPTION_KEY"
    "BACKUP_ENCRYPTION_KEY"
    "MONITORING_API_KEY"
    "HEALTH_CHECK_TOKEN"
  )
  
  # Combine all secrets for processing
  ALL_SECRETS=("${REQUIRED_SECRETS[@]}" "${OPTIONAL_SECRETS[@]}")
  MISSING_REQUIRED=()

  # Load each secret if it exists
  for SECRET in "${ALL_SECRETS[@]}"; do
    SECRET_FILE="/run/secrets/${SECRET}"
    if [ -f "$SECRET_FILE" ]; then
      # Read the secret value and set it as an environment variable
      # Use -n to preserve newlines for multi-line secrets like certificates
      export "$SECRET"="$(cat -n "$SECRET_FILE" | sed 's/^[[:space:]]*[[:digit:]]\+[[:space:]]//g' | tr -d '\r')"
      log_message "Loaded secret: $SECRET"
    elif [[ " ${REQUIRED_SECRETS[@]} " =~ " ${SECRET} " ]]; then
      # Track missing required secrets
      MISSING_REQUIRED+=("$SECRET")
    fi
  done

  # Check for missing required secrets
  if [ ${#MISSING_REQUIRED[@]} -gt 0 ]; then
    log_message "ERROR: Missing required secrets: ${MISSING_REQUIRED[*]}"
    if [ "${STRICT_SECRET_MODE:-true}" = "true" ]; then
      log_message "Exiting due to missing required secrets. Set STRICT_SECRET_MODE=false to continue anyway."
      exit 1
    else
      log_message "Warning: Continuing despite missing required secrets (STRICT_SECRET_MODE=false)"
    fi
  fi

  # Handle environment-specific secrets with deployment color awareness
  ENV_PREFIXES=("STAGING_" "PRODUCTION_")
  
  # If DEPLOYMENT_COLOR is set (blue/green), add color-specific prefix
  if [ -n "$DEPLOYMENT_COLOR" ]; then
    ENV_PREFIXES+=("${DEPLOYMENT_COLOR^^}_")  # Convert to uppercase
  fi
  
  for PREFIX in "${ENV_PREFIXES[@]}"; do
    for SECRET in "${ALL_SECRETS[@]}"; do
      PREFIXED_SECRET="${PREFIX}${SECRET}"
      SECRET_FILE="/run/secrets/${PREFIXED_SECRET}"
      if [ -f "$SECRET_FILE" ]; then
        export "$SECRET"="$(cat -n "$SECRET_FILE" | sed 's/^[[:space:]]*[[:digit:]]\+[[:space:]]//g' | tr -d '\r')"
        log_message "Loaded environment-specific secret: $PREFIXED_SECRET (overriding $SECRET)"
      fi
    done
  done

  # Handle SSH keys for deployment
  SSH_SECRETS=("STAGING_SSH_KEY" "PRODUCTION_SSH_KEY")
  for SSH_SECRET in "${SSH_SECRETS[@]}"; do
    SECRET_FILE="/run/secrets/${SSH_SECRET}"
    if [ -f "$SECRET_FILE" ]; then
      # For SSH keys, we need to set proper permissions
      SSH_DIR="/home/appuser/.ssh"
      if [ ! -d "$SSH_DIR" ]; then
        mkdir -p "$SSH_DIR"
        chmod 700 "$SSH_DIR"
      fi
      
      # Write the key to a file with proper permissions
      KEY_FILE="$SSH_DIR/id_rsa"
      cat "$SECRET_FILE" > "$KEY_FILE"
      chmod 600 "$KEY_FILE"
      
      # Generate public key if it doesn't exist
      if [ ! -f "${KEY_FILE}.pub" ] && command -v ssh-keygen > /dev/null; then
        ssh-keygen -y -f "$KEY_FILE" > "${KEY_FILE}.pub"
        chmod 644 "${KEY_FILE}.pub"
        log_message "Generated public SSH key from private key"
      fi
      
      log_message "Loaded SSH key: $SSH_SECRET to $KEY_FILE"
    fi
  done
  
  # Handle external secret providers
  handle_external_secret_providers
  
  # Handle Kubernetes secrets if running in Kubernetes
  if [ -d "/var/run/secrets/kubernetes.io" ] || [ "${KUBERNETES_SERVICE_HOST:-}" != "" ]; then
    log_message "Kubernetes environment detected, checking for additional secrets"
    handle_kubernetes_secrets
  fi
  
  # Verify secret rotation timestamp if available
  check_secret_rotation
  
  log_message "Secret loading complete"
}

# Function to handle external secret providers
handle_external_secret_providers() {
  if [ -n "$EXTERNAL_SECRET_PROVIDER" ]; then
    log_message "External secret provider configured: $EXTERNAL_SECRET_PROVIDER"
    case "$EXTERNAL_SECRET_PROVIDER" in
      "aws")
        handle_aws_secrets
        ;;
      "vault")
        handle_vault_secrets
        ;;
      "azure")
        handle_azure_secrets
        ;;
      "gcp")
        handle_gcp_secrets
        ;;
      *)
        log_message "Warning: Unknown external secret provider: $EXTERNAL_SECRET_PROVIDER"
        ;;
    esac
  fi
}

# Function to handle AWS Secrets Manager
handle_aws_secrets() {
  if [ -f "/run/secrets/AWS_ACCESS_KEY_ID" ] && [ -f "/run/secrets/AWS_SECRET_ACCESS_KEY" ]; then
    export AWS_ACCESS_KEY_ID="$(cat /run/secrets/AWS_ACCESS_KEY_ID)"
    export AWS_SECRET_ACCESS_KEY="$(cat /run/secrets/AWS_SECRET_ACCESS_KEY)"
    log_message "AWS credentials loaded for secrets management"
    
    # If AWS region is specified, use it
    if [ -f "/run/secrets/AWS_REGION" ]; then
      export AWS_REGION="$(cat /run/secrets/AWS_REGION)"
    elif [ -n "$AWS_REGION" ]; then
      log_message "Using AWS_REGION from environment: $AWS_REGION"
    else
      export AWS_REGION="us-east-1"
      log_message "AWS_REGION not specified, defaulting to us-east-1"
    fi
    
    # If AWS CLI is available and AWS_SECRET_PREFIX is set, fetch secrets
    if command -v aws > /dev/null && [ -n "$AWS_SECRET_PREFIX" ]; then
      log_message "Fetching secrets from AWS Secrets Manager with prefix: $AWS_SECRET_PREFIX"
      
      # List secrets with the specified prefix
      SECRET_LIST=$(aws secretsmanager list-secrets --filter Key=name,Values=$AWS_SECRET_PREFIX --query "SecretList[].Name" --output text 2>/dev/null)
      
      if [ $? -eq 0 ] && [ -n "$SECRET_LIST" ]; then
        for SECRET_NAME in $SECRET_LIST; do
          # Get the secret value
          SECRET_VALUE=$(aws secretsmanager get-secret-value --secret-id "$SECRET_NAME" --query SecretString --output text 2>/dev/null)
          
          if [ $? -eq 0 ] && [ -n "$SECRET_VALUE" ]; then
            # Extract the environment variable name from the secret name
            ENV_NAME=$(echo "$SECRET_NAME" | sed "s/$AWS_SECRET_PREFIX//")
            
            # Set the environment variable
            export "$ENV_NAME"="$SECRET_VALUE"
            log_message "Loaded AWS secret: $SECRET_NAME as $ENV_NAME"
          else
            log_message "Failed to fetch AWS secret: $SECRET_NAME"
          fi
        done
      else
        log_message "No AWS secrets found with prefix: $AWS_SECRET_PREFIX"
      fi
    fi
  else
    log_message "AWS credentials not found, skipping AWS Secrets Manager integration"
  fi
}

# Function to handle HashiCorp Vault
handle_vault_secrets() {
  if [ -f "/run/secrets/VAULT_TOKEN" ]; then
    export VAULT_TOKEN="$(cat /run/secrets/VAULT_TOKEN)"
    log_message "HashiCorp Vault token loaded for secrets management"
    
    # If Vault address is specified, use it
    if [ -f "/run/secrets/VAULT_ADDR" ]; then
      export VAULT_ADDR="$(cat /run/secrets/VAULT_ADDR)"
    elif [ -n "$VAULT_ADDR" ]; then
      log_message "Using VAULT_ADDR from environment: $VAULT_ADDR"
    else
      log_message "VAULT_ADDR not specified, Vault integration may not work"
    fi
    
    # If Vault path is specified and vault CLI is available, fetch secrets
    if command -v vault > /dev/null && [ -n "$VAULT_SECRET_PATH" ]; then
      log_message "Fetching secrets from HashiCorp Vault path: $VAULT_SECRET_PATH"
      
      # Get secrets from Vault
      VAULT_SECRETS=$(vault kv get -format=json "$VAULT_SECRET_PATH" 2>/dev/null)
      
      if [ $? -eq 0 ] && [ -n "$VAULT_SECRETS" ]; then
        # Parse JSON and set environment variables
        # This requires jq, but we'll check for it
        if command -v jq > /dev/null; then
          KEYS=$(echo "$VAULT_SECRETS" | jq -r '.data.data | keys[]')
          
          for KEY in $KEYS; do
            VALUE=$(echo "$VAULT_SECRETS" | jq -r ".data.data.\"$KEY\"")
            export "$KEY"="$VALUE"
            log_message "Loaded Vault secret: $KEY"
          done
        else
          log_message "jq not available, cannot parse Vault secrets"
        fi
      else
        log_message "Failed to fetch secrets from Vault path: $VAULT_SECRET_PATH"
      fi
    fi
  else
    log_message "Vault token not found, skipping HashiCorp Vault integration"
  fi
}

# Function to handle Azure Key Vault
handle_azure_secrets() {
  if [ -f "/run/secrets/AZURE_CLIENT_ID" ] && [ -f "/run/secrets/AZURE_CLIENT_SECRET" ] && [ -f "/run/secrets/AZURE_TENANT_ID" ]; then
    export AZURE_CLIENT_ID="$(cat /run/secrets/AZURE_CLIENT_ID)"
    export AZURE_CLIENT_SECRET="$(cat /run/secrets/AZURE_CLIENT_SECRET)"
    export AZURE_TENANT_ID="$(cat /run/secrets/AZURE_TENANT_ID)"
    log_message "Azure credentials loaded for secrets management"
    
    # If Azure Key Vault name is specified and az CLI is available, fetch secrets
    if command -v az > /dev/null && [ -n "$AZURE_KEYVAULT_NAME" ]; then
      log_message "Fetching secrets from Azure Key Vault: $AZURE_KEYVAULT_NAME"
      
      # Login to Azure
      az login --service-principal -u "$AZURE_CLIENT_ID" -p "$AZURE_CLIENT_SECRET" --tenant "$AZURE_TENANT_ID" > /dev/null
      
      if [ $? -eq 0 ]; then
        # List secrets in the Key Vault
        SECRET_LIST=$(az keyvault secret list --vault-name "$AZURE_KEYVAULT_NAME" --query "[].name" -o tsv 2>/dev/null)
        
        if [ $? -eq 0 ] && [ -n "$SECRET_LIST" ]; then
          for SECRET_NAME in $SECRET_LIST; do
            # Get the secret value
            SECRET_VALUE=$(az keyvault secret show --vault-name "$AZURE_KEYVAULT_NAME" --name "$SECRET_NAME" --query "value" -o tsv 2>/dev/null)
            
            if [ $? -eq 0 ] && [ -n "$SECRET_VALUE" ]; then
              # Set the environment variable
              export "$SECRET_NAME"="$SECRET_VALUE"
              log_message "Loaded Azure Key Vault secret: $SECRET_NAME"
            else
              log_message "Failed to fetch Azure Key Vault secret: $SECRET_NAME"
            fi
          done
        else
          log_message "No secrets found in Azure Key Vault: $AZURE_KEYVAULT_NAME"
        fi
        
        # Logout from Azure
        az logout > /dev/null
      else
        log_message "Failed to login to Azure, skipping Azure Key Vault integration"
      fi
    fi
  else
    log_message "Azure credentials not found, skipping Azure Key Vault integration"
  fi
}

# Function to handle Google Cloud Secret Manager
handle_gcp_secrets() {
  if [ -f "/run/secrets/GOOGLE_APPLICATION_CREDENTIALS" ]; then
    export GOOGLE_APPLICATION_CREDENTIALS="/run/secrets/GOOGLE_APPLICATION_CREDENTIALS"
    log_message "GCP credentials loaded for secrets management"
    
    # If GCP project ID is specified and gcloud CLI is available, fetch secrets
    if command -v gcloud > /dev/null && [ -n "$GCP_PROJECT_ID" ] && [ -n "$GCP_SECRET_PREFIX" ]; then
      log_message "Fetching secrets from GCP Secret Manager in project: $GCP_PROJECT_ID with prefix: $GCP_SECRET_PREFIX"
      
      # Activate service account
      gcloud auth activate-service-account --key-file="$GOOGLE_APPLICATION_CREDENTIALS" > /dev/null
      
      if [ $? -eq 0 ]; then
        # List secrets with the specified prefix
        SECRET_LIST=$(gcloud secrets list --project="$GCP_PROJECT_ID" --filter="name:$GCP_SECRET_PREFIX" --format="value(name)" 2>/dev/null)
        
        if [ $? -eq 0 ] && [ -n "$SECRET_LIST" ]; then
          for SECRET_NAME in $SECRET_LIST; do
            # Get the secret value
            SECRET_VALUE=$(gcloud secrets versions access latest --secret="$SECRET_NAME" --project="$GCP_PROJECT_ID" 2>/dev/null)
            
            if [ $? -eq 0 ] && [ -n "$SECRET_VALUE" ]; then
              # Extract the environment variable name from the secret name
              ENV_NAME=$(echo "$SECRET_NAME" | sed "s/$GCP_SECRET_PREFIX//")
              
              # Set the environment variable
              export "$ENV_NAME"="$SECRET_VALUE"
              log_message "Loaded GCP secret: $SECRET_NAME as $ENV_NAME"
            else
              log_message "Failed to fetch GCP secret: $SECRET_NAME"
            fi
          done
        else
          log_message "No GCP secrets found with prefix: $GCP_SECRET_PREFIX"
        fi
      else
        log_message "Failed to activate GCP service account, skipping GCP Secret Manager integration"
      fi
    fi
  else
    log_message "GCP credentials not found, skipping GCP Secret Manager integration"
  fi
}

# Function to handle Kubernetes secrets
handle_kubernetes_secrets() {
  # Check for Kubernetes service account token
  if [ -f "/var/run/secrets/kubernetes.io/serviceaccount/token" ]; then
    log_message "Kubernetes service account token found"
    
    # If K8S_SECRET_NAMESPACE and K8S_SECRET_NAME are specified, fetch secrets
    if [ -n "$K8S_SECRET_NAMESPACE" ] && [ -n "$K8S_SECRET_NAME" ]; then
      log_message "Fetching Kubernetes secret: $K8S_SECRET_NAME in namespace: $K8S_SECRET_NAMESPACE"
      
      # Get Kubernetes API server
      KUBERNETES_API_SERVER="https://${KUBERNETES_SERVICE_HOST}:${KUBERNETES_SERVICE_PORT}"
      
      # Get service account token
      TOKEN=$(cat /var/run/secrets/kubernetes.io/serviceaccount/token)
      
      # Get CA certificate
      CA_CERT="/var/run/secrets/kubernetes.io/serviceaccount/ca.crt"
      
      # Fetch the secret using curl
      SECRET_DATA=$(curl -s --cacert "$CA_CERT" -H "Authorization: Bearer $TOKEN" \
        "$KUBERNETES_API_SERVER/api/v1/namespaces/$K8S_SECRET_NAMESPACE/secrets/$K8S_SECRET_NAME" 2>/dev/null)
      
      if [ $? -eq 0 ] && [ -n "$SECRET_DATA" ]; then
        # Parse JSON and set environment variables
        # This requires jq, but we'll check for it
        if command -v jq > /dev/null; then
          # Extract data field from the secret
          DATA=$(echo "$SECRET_DATA" | jq -r '.data')
          
          if [ "$DATA" != "null" ]; then
            # Get all keys in the data field
            KEYS=$(echo "$DATA" | jq -r 'keys[]')
            
            for KEY in $KEYS; do
              # Get base64-encoded value
              ENCODED_VALUE=$(echo "$DATA" | jq -r ".[\"$KEY\"]")
              
              # Decode the value
              VALUE=$(echo "$ENCODED_VALUE" | base64 -d)
              
              # Set the environment variable
              export "$KEY"="$VALUE"
              log_message "Loaded Kubernetes secret: $KEY"
            done
          else
            log_message "No data found in Kubernetes secret: $K8S_SECRET_NAME"
          fi
        else
          log_message "jq not available, cannot parse Kubernetes secret"
        fi
      else
        log_message "Failed to fetch Kubernetes secret: $K8S_SECRET_NAME"
      fi
    fi
  fi
}

# Function to check secret rotation
check_secret_rotation() {
  # Check if secret rotation timestamp file exists
  ROTATION_FILE="/run/secrets/SECRET_ROTATION_TIMESTAMP"
  if [ -f "$ROTATION_FILE" ]; then
    ROTATION_TIMESTAMP=$(cat "$ROTATION_FILE")
    CURRENT_TIMESTAMP=$(date +%s)
    
    log_message "Secret rotation timestamp: $ROTATION_TIMESTAMP"
    
    # Calculate age of secrets in days
    SECRET_AGE_SECONDS=$((CURRENT_TIMESTAMP - ROTATION_TIMESTAMP))
    SECRET_AGE_DAYS=$((SECRET_AGE_SECONDS / 86400))
    
    # Log warning if secrets are older than the threshold
    if [ -n "$SECRET_ROTATION_THRESHOLD_DAYS" ] && [ "$SECRET_AGE_DAYS" -gt "$SECRET_ROTATION_THRESHOLD_DAYS" ]; then
      log_message "WARNING: Secrets are $SECRET_AGE_DAYS days old, exceeding threshold of $SECRET_ROTATION_THRESHOLD_DAYS days"
      
      # If strict rotation is enabled, exit
      if [ "${STRICT_ROTATION_MODE:-false}" = "true" ]; then
        log_message "Exiting due to expired secrets (STRICT_ROTATION_MODE=true)"
        exit 1
      fi
    fi
  fi
}

# Function to set up health check endpoints
setup_health_check() {
  log_message "Setting up health check endpoints..."
  
  # Create health check directory if it doesn't exist
  HEALTH_DIR="/app/health"
  if [ ! -d "$HEALTH_DIR" ]; then
    mkdir -p "$HEALTH_DIR"
  fi
  
  # Create main health check file
  cat > "$HEALTH_DIR/status.json" << EOF
{
  "status": "healthy",
  "version": "${APP_VERSION:-unknown}",
  "deployment": "${DEPLOYMENT_COLOR:-production}",
  "timestamp": "$(date -Iseconds)",
  "hostname": "$(hostname)",
  "uptime": "$(uptime -p)"
}
EOF

  # Create liveness probe file
  cat > "$HEALTH_DIR/liveness.json" << EOF
{
  "status": "alive",
  "timestamp": "$(date -Iseconds)",
  "hostname": "$(hostname)"
}
EOF

  # Create readiness probe file
  cat > "$HEALTH_DIR/readiness.json" << EOF
{
  "status": "ready",
  "timestamp": "$(date -Iseconds)",
  "hostname": "$(hostname)",
  "deployment": "${DEPLOYMENT_COLOR:-production}"
}
EOF

  # Create startup probe file
  cat > "$HEALTH_DIR/startup.json" << EOF
{
  "status": "started",
  "timestamp": "$(date -Iseconds)",
  "hostname": "$(hostname)",
  "deployment": "${DEPLOYMENT_COLOR:-production}",
  "version": "${APP_VERSION:-unknown}"
}
EOF
  
  log_message "Health check endpoints configured"
}

# Load secrets into environment variables
log_debug "Starting secrets loading process"
load_secrets
log_debug "Secrets loading completed"

# Set up health check endpoint
log_debug "Setting up health check endpoints"
setup_health_check
log_debug "Health check endpoints configured"

# Print deployment information
log_message "Starting CryoProtect v2 in ${FLASK_ENV:-production} mode"
if [ -n "$DEPLOYMENT_COLOR" ]; then
  log_message "Deployment color: $DEPLOYMENT_COLOR"
fi

# Calculate startup time so far
CURRENT_TIME=$(date +%s.%N)
STARTUP_TIME=$(echo "$CURRENT_TIME - $START_TIME" | bc)
log_message "Entrypoint initialization completed in ${STARTUP_TIME} seconds"

# Optimize Python for production if in production mode
if [ "${FLASK_ENV:-production}" = "production" ]; then
  log_debug "Applying production optimizations"
  # Set Python optimization flags
  export PYTHONOPTIMIZE=2
  # Disable bytecode writing to improve startup time
  export PYTHONDONTWRITEBYTECODE=1
  # Enable fault handler for better crash diagnostics
  export PYTHONFAULTHANDLER=1
fi

# Check if we're running in a container orchestration environment
if [ -n "${KUBERNETES_SERVICE_HOST:-}" ] || [ -n "${DOCKER_SWARM:-}" ]; then
  log_message "Running in container orchestration environment"
  # Apply container-specific optimizations
  export PYTHONUNBUFFERED=1
fi

# Verify required environment variables are set
verify_required_env() {
  local missing=()
  for var in "$@"; do
    if [ -z "${!var}" ]; then
      missing+=("$var")
    fi
  done
  
  if [ ${#missing[@]} -gt 0 ]; then
    log_message "ERROR: Missing required environment variables: ${missing[*]}"
    if [ "${STRICT_ENV_MODE:-true}" = "true" ]; then
      log_message "Exiting due to missing required environment variables. Set STRICT_ENV_MODE=false to continue anyway."
      return 1
    else
      log_message "Warning: Continuing despite missing required environment variables (STRICT_ENV_MODE=false)"
      return 0
    fi
  fi
  
  return 0
}

# Verify required environment variables
REQUIRED_ENV_VARS=("SUPABASE_URL" "SUPABASE_KEY" "SECRET_KEY")
if ! verify_required_env "${REQUIRED_ENV_VARS[@]}"; then
  log_error "Required environment variables missing, exiting"
  exit 1
fi

# Execute the command passed to the container
log_message "Executing: $*"

# Calculate total startup time
CURRENT_TIME=$(date +%s.%N)
TOTAL_STARTUP_TIME=$(echo "$CURRENT_TIME - $START_TIME" | bc)
log_message "Total entrypoint startup time: ${TOTAL_STARTUP_TIME} seconds"

# Execute the command with exec to replace the shell process
exec "$@"