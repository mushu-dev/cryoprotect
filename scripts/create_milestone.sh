#!/bin/bash
# Script to create milestones and assign issues to them

create_milestone() {
  local title="$1"
  local description="$2"
  local due_date="$3"  # Format: YYYY-MM-DD
  
  echo "Creating milestone: $title"
  
  # Create the milestone
  response=$(gh api repos/blueprint-house/CryoProtect/milestones \
    --method POST \
    -f title="$title" \
    -f description="$description" \
    -f due_on="${due_date}T23:59:59Z" \
    --jq '.number')
    
  echo "Created milestone #$response"
  return $response
}

# Create standard project milestones
create_phase_milestones() {
  echo "Creating standard project phase milestones..."
  
  # Calculate dates (30 days apart)
  local today=$(date +%Y-%m-%d)
  local phase1_date=$(date -d "$today + 30 days" +%Y-%m-%d)
  local phase2_date=$(date -d "$today + 60 days" +%Y-%m-%d)
  local phase3_date=$(date -d "$today + 90 days" +%Y-%m-%d)
  local phase4_date=$(date -d "$today + 120 days" +%Y-%m-%d)
  
  # Create the milestones
  create_milestone "Phase 1: Database Architecture & Foundation" \
    "Database schema, connections, and core database functionality." \
    "$phase1_date"
  local phase1=$?
  
  create_milestone "Phase 2: API Implementation & Core Features" \
    "API endpoints, core application logic, and data processing." \
    "$phase2_date"
  local phase2=$?
  
  create_milestone "Phase 3: Deployment & Integration" \
    "Production deployment, integration with external systems." \
    "$phase3_date"
  local phase3=$?
  
  create_milestone "Phase 4: Documentation & Testing" \
    "Comprehensive documentation, testing, and quality assurance." \
    "$phase4_date"
  local phase4=$?
  
  echo "Created standard milestones:"
  echo "Phase 1: #$phase1"
  echo "Phase 2: #$phase2"
  echo "Phase 3: #$phase3"
  echo "Phase 4: #$phase4"
}

# Assign issues to milestones
assign_issues() {
  local milestone_number="$1"
  shift
  local issues=("$@")
  
  echo "Assigning issues to milestone #$milestone_number..."
  
  for issue in "${issues[@]}"; do
    echo "Assigning issue #$issue to milestone #$milestone_number"
    gh issue edit "$issue" --milestone "$milestone_number"
  done
  
  echo "All issues assigned to milestone #$milestone_number"
}

# Main function
main() {
  echo "GitHub Milestone Manager"
  echo "======================="
  
  echo -e "\nOptions:"
  echo "1. Create standard project milestones"
  echo "2. Create custom milestone"
  echo "3. Assign issues to milestone"
  echo "4. Exit"
  
  read -p "Enter your choice (1-4): " choice
  
  case $choice in
    1)
      create_phase_milestones
      ;;
    2)
      read -p "Enter milestone title: " title
      read -p "Enter milestone description: " description
      read -p "Enter due date (YYYY-MM-DD): " due_date
      create_milestone "$title" "$description" "$due_date"
      ;;
    3)
      read -p "Enter milestone number: " milestone_number
      read -p "Enter issue numbers (space-separated): " -a issues
      assign_issues "$milestone_number" "${issues[@]}"
      ;;
    4)
      echo "Exiting."
      exit 0
      ;;
    *)
      echo "Invalid choice. Exiting."
      exit 1
      ;;
  esac
}

# Run main function
main