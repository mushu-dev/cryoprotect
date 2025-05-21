# CryoProtect v2 Implementation Guide

This document provides guidance on using the phase-based execution model for completing the CryoProtect v2 project efficiently.

## Phase-Based Implementation Process

1. **Primary Chat (Current Chat)**
   - For progress updates and phase handoffs
   - High-level oversight of project status
   - Decision-making on when to move to next phase
   - Repository of completed phases and lessons learned

2. **Phase Chat**
   - Dedicated to completing a single cohesive phase
   - Complete autonomy within the phase boundaries
   - Maintains all context about the phase
   - Reports back when phase is complete

## How to Use This Approach

### Starting a New Phase

1. **Create a Phase Chat**
   - Use the template in PHASE_PM_PROMPT.md
   - Customize for the specific phase (2.1, 2.2, etc.)
   - Include all relevant technical context

2. **Set Clear Boundaries**
   - Define exactly what's in scope for the phase
   - List specific files to modify
   - Establish concrete acceptance criteria

3. **Delegate Complete Phase**
   - Give the Phase PM autonomy to organize work
   - Provide access to necessary resources
   - Set expectations for reporting

### During Phase Execution

1. **Minimize Intervention**
   - Let the Phase PM work autonomously
   - Avoid splitting attention across multiple phases
   - Trust the process within phase boundaries

2. **Monitor Progress**
   - Request periodic updates if phase is lengthy
   - Focus on blockers or integration issues
   - Provide guidance if requested

### Phase Completion

1. **Review Deliverables**
   - Verify that all acceptance criteria are met
   - Review code changes and documentation
   - Check test results

2. **Capture Lessons**
   - Document what worked well or could be improved
   - Note any technical debt for future phases
   - Update implementation approach if needed

3. **Plan Next Phase**
   - Decide if ready to move to the next phase
   - Prepare the prompt for the next phase
   - Ensure smooth transition between phases

## Current Phase Plan

| Phase | Description | Status | Dependencies |
|-------|-------------|--------|-------------|
| 2.1 | API Layer Completion | Not Started | - |
| 2.2 | Core Functionality | Not Started | Phase 2.1 |
| 2.3 | User Interface | Not Started | Phase 2.2 |
| 3.1 | Deployment Infrastructure | Not Started | Phase 2.3 |
| 3.2 | Monitoring and Maintenance | Not Started | Phase 3.1 |
| 3.3 | Security | Not Started | Phase 3.2 |
| 4.1 | Documentation | Not Started | Phase 3.3 |
| 4.2 | Knowledge Transfer | Not Started | Phase 4.1 |

## Phase Completion Checklist

For each phase, verify:

- [ ] All implementation tasks completed
- [ ] Tests passing
- [ ] Documentation updated
- [ ] No regressions in functionality
- [ ] Acceptance criteria met
- [ ] Phase summary report provided

## Next Steps

1. Begin with Phase 2.1: API Layer Completion
2. Use the PHASE_PM_PROMPT.md template
3. Complete the phase end-to-end
4. Report back for phase handoff
5. Proceed to Phase 2.2