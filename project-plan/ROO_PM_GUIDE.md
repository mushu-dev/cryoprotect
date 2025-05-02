# Roo Code Project Management Guide for CryoProtect v2

## Overview
This guide provides instructions for using Roo Code to manage the implementation of CryoProtect v2 through an agentic approach. It explains how to interpret and use the project plan documents, assign tasks to specialized agents, and manage the implementation process effectively.

## Project Structure

### Agentic Plan Documents
The project implementation is structured using agentic plan documents located in the `/project-plan/` directory:

- **Overall Plans**
  - `PROJECT_PLAN_README.md` - High-level project overview
  - `ROO_PM_GUIDE.md` - This guide for Roo PM usage

- **Phase Plans**
  - `PHASE_2_1_API_DOCUMENTATION_AGENTIC.md`
  - `PHASE_2_2_PROTOCOL_DESIGNER_AGENTIC.md`
  - `PHASE_2_2_EXPORT_SHARING_AGENTIC.md`
  - (Additional phase plans will be added as needed)

### Plan Structure
Each agentic plan document follows a consistent structure:

1. **Objective** - Overall goal of the phase
2. **Progress Status** - Current completion status
3. **Key Files** - Important files involved
4. **Current Focus** - Specific area of focus
5. **Micro-Tasks** - Individual, self-contained tasks for agents
6. **Implementation Instructions** - Guidelines for Roo PM
7. **Testing Protocol** - Testing approach
8. **Specialized Agent Selection** - Agent type recommendations
9. **Next Phase Planning** - Future directions

## Using Roo Code with Agentic Plans

### Task Assignment Process

1. **Review the Current Phase**
   ```
   I want to implement the following plan: @project-plan/PROJECT_PLAN_README.md. Please analyze our current progress and identify which phase we should focus on next.
   ```

2. **Select a Specific Micro-Task**
   ```
   I want to implement the following plan: @project-plan/PHASE_2_1_API_DOCUMENTATION_AGENTIC.md. Please analyze Task 3A: Molecule Schema Documentation and create a detailed implementation plan for a specialized Schema Agent to execute.
   ```

3. **Assign to Specialized Agent**
   ```
   As a Schema Specialist Agent, implement the Molecule Schema Documentation micro-task as described in PHASE_2_1_API_DOCUMENTATION_AGENTIC.md. Focus specifically on creating a well-documented MoleculeSchema with proper validation and examples.
   ```

4. **Review and Integrate**
   ```
   Please review the implementation of the Molecule Schema Documentation micro-task against the acceptance criteria in PHASE_2_1_API_DOCUMENTATION_AGENTIC.md and verify it meets all requirements.
   ```

### Multi-Agent Coordination

For complex tasks requiring multiple specialized agents:

```
I want to implement Task 2 from @project-plan/PHASE_2_2_PROTOCOL_DESIGNER_AGENTIC.md using multiple specialized agents. Please:

1. Assign Task 2A to a Validation Logic Specialist
2. Assign Task 2B to a Scientific Domain Specialist
3. Coordinate their work to ensure the scientific validation builds on the core validation logic
4. Create a plan for integrating their outputs
```

### Progress Tracking

Track progress on micro-tasks and update the status in the plan documents:

```
Please update the Progress Status section in @project-plan/PHASE_2_1_API_DOCUMENTATION_AGENTIC.md to reflect the completion of Task 3A (Molecule Schema Documentation) and Task 3B (Mixture Schema Documentation).
```

## Specialized Agent Types

Each micro-task is designed for a specific type of specialized agent:

- **Schema Specialist** - For data modeling and schema creation
- **API Documentation Specialist** - For documenting API endpoints
- **Scientific Domain Specialist** - For scientific aspects of cryopreservation
- **Security Specialist** - For authentication and data security
- **Frontend Specialist** - For UI components and interactions
- **Database Specialist** - For database schema and operations
- **Integration Specialist** - For connecting different components

## Best Practices

1. **Single Responsibility**: Each agent should focus on one micro-task at a time
2. **Clear Boundaries**: Define clear inputs and outputs for each agent
3. **Incremental Progress**: Complete and verify one micro-task before moving to the next
4. **Documentation**: Maintain documentation alongside code changes
5. **Testing**: Ensure each micro-task includes appropriate tests
6. **Knowledge Transfer**: Document agent decisions and reasoning for future reference

## Example Roo PM Workflow

1. PM reviews project status and selects Task 3A from PHASE_2_1_API_DOCUMENTATION_AGENTIC.md
2. PM assigns task to Schema Specialist Agent
3. Schema Specialist completes MoleculeSchema implementation
4. PM reviews implementation against acceptance criteria
5. PM updates progress status and selects next task
6. Process repeats for subsequent tasks

## Handling Dependencies

When tasks depend on each other:

```
I notice that Task 3C depends on Task 3A and Task 3B in @project-plan/PHASE_2_1_API_DOCUMENTATION_AGENTIC.md. Please:

1. Verify that Tasks 3A and 3B are complete and meet acceptance criteria
2. If complete, proceed with Task 3C implementation
3. If incomplete, prioritize completing the dependencies first
```

## Conclusion

This agentic approach to project management with Roo Code enables efficient implementation of the CryoProtect v2 project by breaking down complex tasks into manageable micro-tasks that can be handled by specialized agents. By following this guide, you'll be able to coordinate multiple agents, track progress, and ensure successful project completion.