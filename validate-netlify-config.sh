#!/bin/bash
# Validate and fix Netlify configuration

set -e
echo "🔍 Validating Netlify configuration..."

# Define colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Define paths
ROOT_DIR="$(pwd)"
FRONTEND_DIR="$ROOT_DIR/frontend"
ROOT_NETLIFY_TOML="$ROOT_DIR/netlify.toml"
FRONTEND_NETLIFY_TOML="$FRONTEND_DIR/netlify.toml"

# Functions
validate_file() {
  local file="$1"
  if [[ -f "$file" ]]; then
    echo -e "${GREEN}✅ $file exists${NC}"
    return 0
  else
    echo -e "${RED}❌ $file does not exist${NC}"
    return 1
  fi
}

check_env_var() {
  local file="$1"
  local var_name="$2"
  
  if grep -q "$var_name" "$file"; then
    echo -e "${GREEN}✅ $var_name is set in $file${NC}"
    return 0
  else
    echo -e "${YELLOW}⚠️ $var_name is missing in $file${NC}"
    return 1
  fi
}

# Step 1: Verify netlify.toml files
echo -e "\n${YELLOW}Checking netlify.toml files...${NC}"
validate_file "$ROOT_NETLIFY_TOML"
validate_file "$FRONTEND_NETLIFY_TOML"

# Step 2: Check which netlify.toml is being used by Netlify CLI
echo -e "\n${YELLOW}Checking which netlify.toml is used by Netlify CLI...${NC}"
NETLIFY_SITE_INFO=$(netlify status)
USED_TOML=$(echo "$NETLIFY_SITE_INFO" | grep "Netlify TOML" | awk '{print $NF}')

echo -e "Netlify is using: ${GREEN}$USED_TOML${NC}"

# Step 3: Check for required environment variables in the used netlify.toml
echo -e "\n${YELLOW}Checking environment variables in the used netlify.toml...${NC}"
USED_TOML_PATH=""

if [[ "$USED_TOML" == *"frontend/netlify.toml"* ]]; then
  USED_TOML_PATH="$FRONTEND_NETLIFY_TOML"
else
  USED_TOML_PATH="$ROOT_NETLIFY_TOML"
fi

REQUIRED_VARS=(
  "NEXT_PUBLIC_NETLIFY"
  "NEXT_PUBLIC_API_URL"
  "NEXT_PUBLIC_ENVIRONMENT"
)

MISSING_VARS=()
for var in "${REQUIRED_VARS[@]}"; do
  if ! check_env_var "$USED_TOML_PATH" "$var"; then
    MISSING_VARS+=("$var")
  fi
done

# Step 4: Check for API redirects in the used netlify.toml
echo -e "\n${YELLOW}Checking API redirects in the used netlify.toml...${NC}"
if grep -q "/api/v1/health/connectivity" "$USED_TOML_PATH"; then
  echo -e "${GREEN}✅ API connectivity endpoint redirect found${NC}"
else
  echo -e "${YELLOW}⚠️ API connectivity endpoint redirect not found${NC}"
fi

if grep -q "404.html" "$USED_TOML_PATH"; then
  echo -e "${GREEN}✅ 404 tracking redirect found${NC}"
else
  echo -e "${YELLOW}⚠️ 404 tracking redirect not found${NC}"
fi

# Step 5: Verify analytics components
echo -e "\n${YELLOW}Verifying analytics components...${NC}"
ANALYTICS_FILES=(
  "$FRONTEND_DIR/src/app/netlify-analytics.js"
  "$FRONTEND_DIR/src/hooks/useAnalytics.ts"
  "$FRONTEND_DIR/src/components/analytics/AnalyticsProvider.tsx"
  "$FRONTEND_DIR/src/components/analytics/AnalyticsConsent.tsx"
)

for file in "${ANALYTICS_FILES[@]}"; do
  validate_file "$file"
done

# Step 6: Check integration in layout.tsx
echo -e "\n${YELLOW}Checking integration in layout.tsx...${NC}"
LAYOUT_FILE="$FRONTEND_DIR/src/app/layout.tsx"

validate_file "$LAYOUT_FILE"
if grep -q "NetlifyAnalytics" "$LAYOUT_FILE"; then
  echo -e "${GREEN}✅ NetlifyAnalytics component is imported in layout.tsx${NC}"
else
  echo -e "${RED}❌ NetlifyAnalytics component is not imported in layout.tsx${NC}"
fi

if grep -q "NEXT_PUBLIC_NETLIFY" "$LAYOUT_FILE"; then
  echo -e "${GREEN}✅ Conditional Netlify analytics rendering found in layout.tsx${NC}"
else
  echo -e "${RED}❌ Conditional Netlify analytics rendering not found in layout.tsx${NC}"
fi

# Step 7: Validate CSP in netlify.toml
echo -e "\n${YELLOW}Validating Content-Security-Policy...${NC}"
if grep -q "plausible.io" "$USED_TOML_PATH"; then
  echo -e "${GREEN}✅ Plausible domain is included in Content-Security-Policy${NC}"
else
  echo -e "${YELLOW}⚠️ Plausible domain is missing in Content-Security-Policy${NC}"
fi

# Step 8: Final report and suggestions
echo -e "\n${YELLOW}==== Validation Summary ====${NC}"

if [[ ${#MISSING_VARS[@]} -eq 0 ]]; then
  echo -e "${GREEN}✅ All required environment variables are set${NC}"
else
  echo -e "${YELLOW}⚠️ Missing environment variables: ${MISSING_VARS[*]}${NC}"
  echo -e "Consider adding these to the Netlify dashboard or to $USED_TOML_PATH"
fi

# Check Netlify CLI deploy status
echo -e "\n${YELLOW}Testing Netlify deployment status...${NC}"
if netlify sites:list | grep -q "cryoprotect"; then
  echo -e "${GREEN}✅ Netlify site 'cryoprotect' exists${NC}"
  echo -e "You can deploy using: ${GREEN}netlify deploy${NC}"
else
  echo -e "${YELLOW}⚠️ Netlify site may not be properly set up${NC}"
  echo -e "Run: ${GREEN}netlify init${NC} to set up your site"
fi

echo -e "\n${GREEN}Validation complete!${NC}"
echo -e "To deploy your site with analytics, run: ${GREEN}netlify deploy --prod${NC}"
echo -e "To purchase Netlify Analytics: Go to https://app.netlify.com/projects/cryoprotect/analytics"