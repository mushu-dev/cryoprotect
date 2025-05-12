#!/bin/bash
# Script to consolidate GitHub issues for project organization

# Set terminal colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}===== GitHub Issue Consolidation Script =====${NC}"
echo "This script will organize GitHub issues for Cursor migration"

# Create milestone for project phases
echo -e "\n${YELLOW}Creating milestones for project phases...${NC}"
gh api repos/:owner/:repo/milestones -X POST -F title="Phase 1: Technical Foundation" -F state="open" -F description="Database Architecture & Authentication System"
gh api repos/:owner/:repo/milestones -X POST -F title="Phase 2: Feature Completion" -F state="open" -F description="API Layer, Core Functionality, User Interface"
gh api repos/:owner/:repo/milestones -X POST -F title="Phase 3: Production Readiness" -F state="open" -F description="Deployment, Monitoring, Security"
gh api repos/:owner/:repo/milestones -X POST -F title="Phase 4: Documentation" -F state="open" -F description="Documentation and Knowledge Transfer"
gh api repos/:owner/:repo/milestones -X POST -F title="ChEMBL Integration" -F state="open" -F description="ChEMBL Database Integration"

# Get milestone IDs
PHASE1=$(gh api repos/:owner/:repo/milestones --jq '.[] | select(.title=="Phase 1: Technical Foundation") | .number')
PHASE2=$(gh api repos/:owner/:repo/milestones --jq '.[] | select(.title=="Phase 2: Feature Completion") | .number')
PHASE3=$(gh api repos/:owner/:repo/milestones --jq '.[] | select(.title=="Phase 3: Production Readiness") | .number')
PHASE4=$(gh api repos/:owner/:repo/milestones --jq '.[] | select(.title=="Phase 4: Documentation") | .number')
CHEMBL=$(gh api repos/:owner/:repo/milestones --jq '.[] | select(.title=="ChEMBL Integration") | .number')

# Create new comprehensive labels
echo -e "\n${YELLOW}Creating more descriptive labels...${NC}"
gh label create "priority:high" --color "#FF0000" --description "High priority task"
gh label create "priority:medium" --color "#FFA500" --description "Medium priority task"
gh label create "priority:low" --color "#FFFF00" --description "Low priority task"
gh label create "component:database" --color "#0075CA" --description "Related to database functionality"
gh label create "component:api" --color "#A2EEEF" --description "Related to API functionality"
gh label create "component:ui" --color "#7057FF" --description "Related to UI functionality"
gh label create "component:chembl" --color "#D876E3" --description "Related to ChEMBL integration"
gh label create "component:authentication" --color "#008672" --description "Related to authentication system"
gh label create "component:testing" --color "#FBCA04" --description "Related to testing infrastructure"
gh label create "cursor-migration" --color "#1D76DB" --description "Part of Cursor migration effort"

# Consolidate issues by category
echo -e "\n${YELLOW}Consolidating ChEMBL integration issues...${NC}"
# Update ChEMBL integration issues
for id in 227 228 234 196 201 225 226; do
  gh issue edit $id --add-label "component:chembl,cursor-migration" --milestone $CHEMBL
done

# Assign implementation and validation pairs
echo -e "\n${YELLOW}Linking implementation issues with their validation counterparts...${NC}"
# Link issues 230 and 236
gh issue edit 230 --add-label "component:database" --milestone $PHASE1
gh issue edit 236 --add-label "component:testing" --milestone $PHASE1
gh issue comment 230 -b "This issue is linked to validation issue #236. Both issues will be migrated to Cursor for continued development."
gh issue comment 236 -b "This issue is validating the implementation in #230. Both issues will be migrated to Cursor for continued development."

# Link issues 229 and 235
gh issue edit 229 --add-label "component:database" --milestone $PHASE1
gh issue edit 235 --add-label "component:testing" --milestone $PHASE1
gh issue comment 229 -b "This issue is linked to validation issue #235. Both issues will be migrated to Cursor for continued development."
gh issue comment 235 -b "This issue is validating the implementation in #229. Both issues will be migrated to Cursor for continued development."

# Link issues 231 and 239
gh issue edit 231 --add-label "component:api" --milestone $PHASE2
gh issue edit 239 --add-label "component:testing" --milestone $PHASE2
gh issue comment 231 -b "This issue is linked to validation issue #239. Both issues will be migrated to Cursor for continued development."
gh issue comment 239 -b "This issue is validating the implementation in #231. Both issues will be migrated to Cursor for continued development."

# Link issues 237 and 238
gh issue edit 237 --add-label "component:api,type:bugfix" --milestone $PHASE2
gh issue edit 238 --add-label "component:testing" --milestone $PHASE2
gh issue comment 237 -b "This issue is linked to validation issue #238. Both issues will be migrated to Cursor for continued development."
gh issue comment 238 -b "This issue is validating the implementation in #237. Both issues will be migrated to Cursor for continued development."

# Consolidate Windows/Linux related tasks
echo -e "\n${YELLOW}Organizing Windows/Linux compatibility issues...${NC}"
for id in 209 241 242; do
  gh issue edit $id --add-label "component:testing,cursor-migration" --milestone $PHASE3
done
gh issue comment 209 -b "This issue is a parent for Windows compatibility issues #241 and #242. All will be migrated to Cursor."

# Create a meta-issue for Cursor migration
echo -e "\n${YELLOW}Creating a main Cursor migration tracking issue...${NC}"
gh issue create --title "Cursor IDE Migration Framework" \
  --body "# Cursor IDE Migration Framework

This issue tracks our migration from RooRoo/Claude Code to Cursor IDE.

## Migration Phases
1. **Issue Organization**: Consolidate GitHub issues and apply appropriate labels
2. **Project Structure Review**: Ensure project works correctly in Cursor
3. **Tool Update**: Update scripts and tools to work with Cursor 
4. **Documentation**: Update documentation to reflect the new workflow

## Related Issues
- Windows/Linux Compatibility: #209, #241, #242
- ChEMBL Integration: #227, #228, #234, #196, #201, #225, #226

## Current Status
* [x] Created milestone structure
* [x] Created component-based labels
* [x] Consolidated related issues
* [ ] Updated documentation
* [ ] Tested project in Cursor IDE
* [ ] Completed transition

## Next Steps
1. Review open PRs and ensure they're properly labeled
2. Update README with Cursor-specific instructions
3. Complete any remaining RooRoo tasks before full migration" \
  --label "priority:high,cursor-migration" \
  --milestone $PHASE1

echo -e "\n${GREEN}Issue consolidation complete! Issues are now organized for Cursor migration.${NC}"