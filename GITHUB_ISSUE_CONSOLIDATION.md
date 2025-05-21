# GitHub Issue Consolidation Guide

## Purpose
This document provides guidelines for consolidating duplicate or fragmented issues in the CryoProtect repository, with a focus on maintaining clear organization and traceability.

## When to Consolidate Issues
Issues should be consolidated when:
1. Multiple issues describe the same work
2. Several small issues are better managed as parts of a larger issue
3. Related features are scattered across numerous issues
4. A clear parent-child relationship exists between issues

## Consolidation Process

### 1. Identify Related Issues
- Use search filters to find issues with similar titles or topics
- Look for issues with overlapping descriptions
- Check for issues that reference the same components or features
- Prioritize consolidating active/open issues first

### 2. Select a Master Issue
- Choose the most comprehensive or clearly defined issue as the master
- Prefer issues that already have:
  - Detailed descriptions
  - Clear acceptance criteria
  - Proper labels and milestone assignments
- If no suitable issue exists, create a new master issue

### 3. Enhance the Master Issue
- Update the description to include all relevant information
- Add a "Consolidated Issues" section listing all related issues
- Ensure the description includes:
  - Comprehensive problem/feature overview
  - Complete implementation details
  - All necessary acceptance criteria
  - Combined dependencies
  - Merged technical approach details

### 4. Reference the Master Issue
- Add comments to all related issues pointing to the master issue
- Use consistent wording: "This issue has been consolidated into the master issue #XXX for [feature]. Please refer to that issue for the complete implementation plan and tracking."
- Do not close the related issues if they contain unique information or ongoing discussion

### 5. Label Management
- Ensure the master issue has appropriate:
  - Area label(s)
  - Type label
  - Status label
  - Priority label
- Consider adding a custom "consolidated" label to identify issues that have undergone this process

## Example Consolidation

### ChEMBL/PubChem Import Consolidation
The first major consolidation was the ChEMBL and PubChem data import pipeline:

1. **Master Issue**: #243 "Implement unified ChEMBL and PubChem data import pipeline"
2. **Consolidated Issues**:
   - #232: Unified Import Script Development
   - #227: ChEMBL Integration: Overall Plan
   - #196: ChEMBL Cryoprotectant Data Population: Overall Plan
   - #211: Fix ChEMBL Population Script and Perform Full Database Population
   - #201: ChEMBL Database Full Population with Updated Connection
   - #148: Orchestrate Complete ChEMBL Integration Process
   - #178: Main PubChem Population Script

3. **Enhancement**: The master issue was updated with:
   - Comprehensive implementation plan
   - Combined objectives
   - Merged technical details
   - Complete acceptance criteria
   - Reference to all consolidated issues

4. **References**: Comments were added to all related issues pointing to the master issue

## Best Practices
1. **Maintain Traceability**: Always reference the consolidated issues
2. **Preserve Context**: Don't delete information; merge it into the master issue
3. **Clear Communication**: Keep stakeholders informed about consolidation
4. **Systematic Approach**: Tackle one category of issues at a time
5. **Milestone Alignment**: Ensure the master issue is assigned to the appropriate milestone

## Next Consolidation Targets
1. **RLS Implementation Issues**: Consolidate the various RLS implementation tasks
2. **Testing Framework Issues**: Combine the scattered testing framework tasks
3. **Documentation Issues**: Merge related documentation tasks</content>