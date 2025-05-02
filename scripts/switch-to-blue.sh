#!/bin/bash
set -e

# Switch traffic to blue environment
# This script updates the NGINX configuration to route traffic to the blue environment

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
NGINX_CONF_DIR="$PROJECT_DIR/nginx/conf.d"
ACTIVE_CONF="$NGINX_CONF_DIR/active.conf"

# Check if blue environment is healthy
if ! docker ps --filter "name=cryoprotect-blue" --filter "health=healthy" --quiet | grep -q .; then
  echo "Error: Blue environment is not healthy. Cannot switch traffic."
  exit 1
fi

# Create new active.conf pointing to blue
cat > "$ACTIVE_CONF" << EOF
# Active environment is blue
# Last updated: $(date)

server {
    listen 80;
    server_name localhost;

    location / {
        proxy_pass http://blue;
        proxy_set_header Host \$host;
        proxy_set_header X-Real-IP \$remote_addr;
        proxy_set_header X-Forwarded-For \$proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto \$scheme;
    }

    # Health check endpoint
    location /health {
        proxy_pass http://blue/health;
        proxy_set_header Host \$host;
        proxy_set_header X-Real-IP \$remote_addr;
        access_log off;
    }
}
EOF

# Reload NGINX configuration
echo "Reloading NGINX configuration..."
docker exec "$(docker ps --filter "name=nginx" --quiet)" nginx -s reload || {
  echo "Error: Failed to reload NGINX configuration."
  exit 1
}

echo "Traffic successfully switched to blue environment!"