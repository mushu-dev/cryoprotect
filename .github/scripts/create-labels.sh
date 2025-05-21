#!/bin/bash

# Create labels script
# Creates all labels needed for the GitHub repository

echo "Creating labels for the GitHub repository..."

# Type labels
gh label create "type:epic" --color "6E49CB" --description "Large initiative containing multiple issues"
gh label create "type:feature" --color "b60205" --description "New feature implementation"
gh label create "type:bugfix" --color "88E3AE" --description "Bug fix"
gh label create "type:refactor" --color "F7C4A7" --description "Code refactoring"
gh label create "type:documentation" --color "B01778" --description "Documentation updates"
gh label create "type:validation" --color "6BA5AF" --description "Validation and testing tasks"

# Status labels
gh label create "status:planning" --color "FEF2C0" --description "In planning/design phase"
gh label create "status:ready" --color "0E8A16" --description "Ready for implementation"
gh label create "status:in-progress" --color "BFDADC" --description "Currently being worked on"
gh label create "status:needs-review" --color "008F91" --description "Implementation complete, needs review"
gh label create "status:completed" --color "806362" --description "Fully completed and validated"
gh label create "status:blocked" --color "D93F0B" --description "Blocked by dependencies or issues"

# Area labels
gh label create "area:database" --color "0052CC" --description "Database-related issues"
gh label create "area:api" --color "52CC00" --description "API-related issues"
gh label create "area:ui" --color "CC5200" --description "User interface issues"
gh label create "area:auth" --color "5200CC" --description "Authentication and security issues"
gh label create "area:testing" --color "00CC52" --description "Testing and QA issues"
gh label create "area:chembl" --color "CC0052" --description "ChEMBL integration issues"
gh label create "area:docs" --color "00CCCC" --description "Documentation issues"
gh label create "area:devops" --color "CCCC00" --description "Infrastructure and DevOps issues"

# Priority labels
gh label create "priority:high" --color "B60205" --description "High priority"
gh label create "priority:medium" --color "F2CB05" --description "Medium priority" 
gh label create "priority:low" --color "4CBB17" --description "Low priority"

# Special labels
gh label create "stale" --color "D3D3D3" --description "No activity for 30+ days"
gh label create "duplicate" --color "EEEEEE" --description "Duplicate of another issue"
gh label create "ai-generated" --color "FF00FF" --description "AI-generated test issue"

echo "Labels created successfully!"