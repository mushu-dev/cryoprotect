#!/bin/bash
# test_connectivity.sh
# Comprehensive connectivity test script for CryoProtect
# Tests both IPv4 and IPv6 connectivity to Supabase and provides detailed diagnostics

set -e

# ANSI color codes for better readability
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
MAGENTA='\033[0;35m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m' # No Color

# Print header
echo -e "${BOLD}${BLUE}=============================================${NC}"
echo -e "${BOLD}${BLUE}  CryoProtect Dual-Stack Connectivity Test   ${NC}"
echo -e "${BOLD}${BLUE}=============================================${NC}"

# Check for required arguments
if [ -z "$1" ]; then
    echo -e "${YELLOW}Usage: $0 <supabase_url> [database_name]${NC}"
    echo -e "${YELLOW}Example: $0 https://yourproject.supabase.co postgres${NC}"
    exit 1
fi

SUPABASE_URL="$1"
DB_NAME="${2:-postgres}"

# Extract hostname from Supabase URL
HOSTNAME=$(echo "$SUPABASE_URL" | sed -E 's|^https?://||' | sed -E 's|/.*$||' | sed -E 's|:[0-9]+$||')
DB_PORT=5432

echo -e "${CYAN}Testing connectivity to:${NC} $HOSTNAME"
echo -e "${CYAN}Database:${NC} $DB_NAME"
echo

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Check for required tools
echo -e "${BOLD}Checking for required tools...${NC}"
MISSING_TOOLS=0

for tool in nc host dig nslookup curl openssl ping iputils-ping psql; do
    if ! command_exists "$tool"; then
        echo -e "${RED}✗ $tool not found${NC}"
        MISSING_TOOLS=$((MISSING_TOOLS+1))
    else
        echo -e "${GREEN}✓ $tool found${NC}"
    fi
done

if [ $MISSING_TOOLS -gt 0 ]; then
    echo -e "${YELLOW}Some tools are missing. Install them with:${NC}"
    echo -e "  sudo dnf install -y netcat hostname bind-utils curl openssl iputils postgresql"
    echo
fi

# Function to print section header
section() {
    echo
    echo -e "${BOLD}${MAGENTA}● $1${NC}"
    echo -e "${MAGENTA}$(printf '%.0s─' $(seq 1 50))${NC}"
}

# 1. System network stack information
section "System Network Configuration"

echo -e "${BOLD}Kernel IP configuration:${NC}"
sysctl -a 2>/dev/null | grep -E 'ipv6|ip_forward' | sort | head -10

echo -e "\n${BOLD}Network interfaces:${NC}"
ip -br addr show

echo -e "\n${BOLD}Default routes:${NC}"
ip route show default
ip -6 route show default 2>/dev/null || echo -e "${YELLOW}No IPv6 default route${NC}"

# 2. DNS resolution tests
section "DNS Resolution Tests"

echo -e "${BOLD}System DNS configuration:${NC}"
cat /etc/resolv.conf | grep -v '^#' | grep .

echo -e "\n${BOLD}IPv4 address records (A):${NC}"
host -t A "$HOSTNAME" || echo -e "${YELLOW}No IPv4 records found${NC}"

echo -e "\n${BOLD}IPv6 address records (AAAA):${NC}"
host -t AAAA "$HOSTNAME" || echo -e "${YELLOW}No IPv6 records found${NC}"

echo -e "\n${BOLD}All IP addresses for hostname:${NC}"
getent ahosts "$HOSTNAME"

# 3. Basic connectivity tests
section "Basic Connectivity Tests"

echo -e "${BOLD}ICMP ping test (IPv4):${NC}"
ping -c 3 -4 "$HOSTNAME" || echo -e "${YELLOW}IPv4 ping failed or blocked${NC}"

echo -e "\n${BOLD}ICMP ping test (IPv6):${NC}"
ping -c 3 -6 "$HOSTNAME" 2>/dev/null || echo -e "${YELLOW}IPv6 ping failed, blocked, or not available${NC}"

# 4. Port connectivity tests
section "Port Connectivity Tests"

echo -e "${BOLD}Testing PostgreSQL connectivity (IPv4):${NC}"
nc -zv -4 "$HOSTNAME" $DB_PORT -w 5 || echo -e "${RED}IPv4 connection to port $DB_PORT failed${NC}"

echo -e "\n${BOLD}Testing PostgreSQL connectivity (IPv6):${NC}"
nc -zv -6 "$HOSTNAME" $DB_PORT -w 5 2>/dev/null || echo -e "${YELLOW}IPv6 connection to port $DB_PORT failed or not available${NC}"

echo -e "\n${BOLD}Testing HTTPS connectivity (IPv4):${NC}"
nc -zv -4 "$HOSTNAME" 443 -w 5 || echo -e "${RED}IPv4 connection to port 443 failed${NC}"

echo -e "\n${BOLD}Testing HTTPS connectivity (IPv6):${NC}"
nc -zv -6 "$HOSTNAME" 443 -w 5 2>/dev/null || echo -e "${YELLOW}IPv6 connection to port 443 failed or not available${NC}"

# 5. SSL/TLS verification
section "SSL/TLS Verification"

echo -e "${BOLD}Checking HTTPS certificate:${NC}"
echo | openssl s_client -connect "$HOSTNAME:443" 2>/dev/null | openssl x509 -noout -dates -subject -issuer || echo -e "${RED}SSL certificate verification failed${NC}"

echo -e "\n${BOLD}Checking PostgreSQL SSL:${NC}"
echo | openssl s_client -connect "$HOSTNAME:$DB_PORT" 2>/dev/null | grep "Verify return code" || echo -e "${YELLOW}PostgreSQL SSL verification failed or not available${NC}"

# 6. Traceroute (if available)
section "Network Path Analysis"

if command_exists traceroute; then
    echo -e "${BOLD}Traceroute to host (IPv4):${NC}"
    traceroute -n -4 "$HOSTNAME" || echo -e "${YELLOW}IPv4 traceroute failed${NC}"
    
    if command_exists traceroute6; then
        echo -e "\n${BOLD}Traceroute to host (IPv6):${NC}"
        traceroute6 -n "$HOSTNAME" 2>/dev/null || echo -e "${YELLOW}IPv6 traceroute failed or not available${NC}"
    fi
else
    echo -e "${YELLOW}Traceroute not installed. Install with: sudo dnf install -y traceroute${NC}"
fi

# 7. Database connectivity tests (if PostgreSQL client is available)
section "Database Connection Tests"

if command_exists psql; then
    echo -e "${BOLD}Testing PostgreSQL connection string construction:${NC}"
    IPV4_ADDR=$(getent ahostsv4 "$HOSTNAME" | grep STREAM | head -1 | awk '{print $1}')
    IPV6_ADDR=$(getent ahostsv6 "$HOSTNAME" | grep STREAM | head -1 | awk '{print $1}')
    
    if [ -n "$IPV4_ADDR" ]; then
        echo -e "${GREEN}✓ IPv4 address resolved:${NC} $IPV4_ADDR"
        echo -e "  Connection string would be: postgresql://postgres:password@$IPV4_ADDR:$DB_PORT/$DB_NAME?sslmode=require"
    else
        echo -e "${RED}✗ Failed to resolve IPv4 address${NC}"
    fi
    
    if [ -n "$IPV6_ADDR" ]; then
        echo -e "${GREEN}✓ IPv6 address resolved:${NC} $IPV6_ADDR"
        echo -e "  Connection string would be: postgresql://postgres:password@[$IPV6_ADDR]:$DB_PORT/$DB_NAME?sslmode=require"
    else
        echo -e "${YELLOW}✗ Failed to resolve IPv6 address${NC}"
    fi
    
    echo -e "\n${BOLD}To test actual database connection, use:${NC}"
    echo "PGPASSWORD=your_password psql -h $HOSTNAME -p $DB_PORT -U postgres -d $DB_NAME"
    
    if [ -n "$IPV4_ADDR" ]; then
        echo "# or for explicit IPv4:"
        echo "PGPASSWORD=your_password psql -h $IPV4_ADDR -p $DB_PORT -U postgres -d $DB_NAME"
    fi
    
    if [ -n "$IPV6_ADDR" ]; then
        echo "# or for explicit IPv6:"
        echo "PGPASSWORD=your_password psql -h $IPV6_ADDR -p $DB_PORT -U postgres -d $DB_NAME"
    fi
else
    echo -e "${YELLOW}PostgreSQL client not installed. Install with: sudo dnf install -y postgresql${NC}"
fi

# 8. API connectivity test
section "API Connectivity Test"

echo -e "${BOLD}Testing Supabase REST API (IPv4):${NC}"
curl -s -4 -o /dev/null -w "%{http_code}\n" "$SUPABASE_URL/rest/v1/" -H "apikey: REPLACE_WITH_ANON_KEY" || echo -e "${RED}IPv4 API connection failed${NC}"

echo -e "\n${BOLD}Testing Supabase REST API (IPv6):${NC}"
curl -s -6 -o /dev/null -w "%{http_code}\n" "$SUPABASE_URL/rest/v1/" -H "apikey: REPLACE_WITH_ANON_KEY" 2>/dev/null || echo -e "${YELLOW}IPv6 API connection failed or not available${NC}"

# 9. Summary and recommendations
section "Connectivity Summary"

IPV4_WORKS=0
IPV6_WORKS=0

# Check IPv4 connectivity
if ping -c 1 -4 "$HOSTNAME" > /dev/null 2>&1 && nc -z -4 "$HOSTNAME" $DB_PORT -w 5 > /dev/null 2>&1; then
    IPV4_WORKS=1
fi

# Check IPv6 connectivity
if ping -c 1 -6 "$HOSTNAME" > /dev/null 2>&1 && nc -z -6 "$HOSTNAME" $DB_PORT -w 5 > /dev/null 2>&1; then
    IPV6_WORKS=1
fi

if [ $IPV4_WORKS -eq 1 ] && [ $IPV6_WORKS -eq 1 ]; then
    echo -e "${GREEN}✓ Both IPv4 and IPv6 connectivity available${NC}"
    echo -e "${BOLD}Recommendation:${NC} Use dual-stack connection strategy with fallback"
elif [ $IPV4_WORKS -eq 1 ]; then
    echo -e "${YELLOW}⚠ Only IPv4 connectivity available${NC}"
    echo -e "${BOLD}Recommendation:${NC} Use IPv4-only connection configuration"
elif [ $IPV6_WORKS -eq 1 ]; then
    echo -e "${YELLOW}⚠ Only IPv6 connectivity available${NC}"
    echo -e "${BOLD}Recommendation:${NC} Use IPv6-only connection configuration"
else
    echo -e "${RED}✗ No connectivity available to the database${NC}"
    echo -e "${BOLD}Recommendation:${NC} Check your network configuration and firewall rules"
fi

echo
echo -e "${BOLD}${BLUE}=============================================${NC}"
echo -e "${BOLD}${BLUE}       Connectivity Test Completed           ${NC}"
echo -e "${BOLD}${BLUE}=============================================${NC}"