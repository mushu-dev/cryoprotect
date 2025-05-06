# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [v0.1.2] - 2025-04-28
### Changed
- **Workflow Coordinator Dispatch Logic Overhauled (Hybrid Model):** The `workflow-coordinator` (now v5 implicitly) uses a fundamentally different dispatch mechanism:
    - **Type-Based Mapping:** Now uses a fixed mapping (defined in its instructions) to determine the target mode based on the task `type`.
    - **Hybrid Delegation:**
        - Specialist Tasks (`tech-design`, `ui-design`, `ux-design`, `validation`, `test-execution`, `documentation-*`): Delegated to the corresponding custom specialist modes (`solution-architect`, `ux-specialist`, `guardian-validator`, `docu-crafter`) using `<new_task>` and passing the `taskStateFile`. Relies on these specialists to update the overview status directly upon completion/error.
        - Coding/Fixing Tasks (`feature`, `refactor`, `chore`, `bugfix`): Explicitly delegated to the built-in `code` and `debug` modes using `<new_task>` **without** passing the `taskStateFile`.
    - **Completion Handling for Built-in Modes:** Shifted from primary reliance on inference to **explicitly requiring user confirmation** via `<ask_followup_question>` before updating the status of tasks handled by built-in modes to `Implemented`.
- **Strategic Planner Task Definition:**
    - Now assigns tasks using a **predefined, fixed list of `type` values**.
    - No longer defines or plans an `intendedMode`; the `type` alone dictates the Coordinator's dispatch action.
- **Specialist Agent Instructions Refined:** Instructions for `solution-architect`, `ux-specialist`, `guardian-validator`, and `docu-crafter` updated to:
    - Explicitly reference the specific `type`(s) of tasks they handle.
    - Reinforce their responsibility to directly update **only** the `project_overview.json` status for their assigned task upon completion ('Done', 'Validated', 'Failed') or error ('Error', 'Blocked-Debug').
- **State Schema (`project_overview.json`):**
    - The list of allowed values for the `type` field is now explicitly defined and fixed within the Planner/Coordinator instructions.
    - The `Implemented` status remains necessary for the multi-step flow involving built-in modes and user confirmation/testing.

### Removed
- **Coordinator Logic:** Removed the primary reliance on *inferring* completion status for tasks delegated to built-in modes. Replaced with an explicit user confirmation loop.


## [v0.1.1] - 2025-04-28

### Changed
- **Workflow Coordinator Delegation Refined:** Updated `customInstructions` for `workflow-coordinator` (v3 implicitly) to refine task delegation logic:
    - `tech-design` tasks are now directed to `solution-architect`.
    - `ui-design` and `ux-design` tasks are now directed to `ux-specialist`.
    - Both `test_execution` and `validation` tasks are now explicitly directed to `guardian-validator`.

## [v0.1.0] - 2025-04-28

### BREAKING CHANGE
- **Major Architectural Shift:** Replaced the single `master-orchestrator` agent with a two-part system (`strategic-planner` and `workflow-coordinator`) to separate planning and execution management.
    - **`strategic-planner`:** Now handles initial goal interpretation, high-level planning, and initial state file setup (`project_overview.json`, `.state/tasks/`). It then hands off execution management.
    - **`workflow-coordinator`:** Now manages the execution phase, monitoring `project_overview.json`, delegating tasks, interpreting outcomes, handling user test decisions, and managing errors.
- **Delegation Model Changed:** The `workflow-coordinator` now delegates implementation/fixing tasks (`feature`, `bugfix`, `refactor`, `chore`) primarily to the **built-in `code` and `debug` modes** instead of the previous dedicated `apex-implementer` agent. Specialist tasks (design, validation, docs) are still delegated to specialist agents.
- **Removed Agent:** Removed the `apex-implementer` specialist agent. Its responsibilities are now handled by the built-in `code` and `debug` modes under the direction of the `workflow-coordinator`.
- **State Management for Built-in Modes:** The `workflow-coordinator` is now responsible for inferring the completion status of tasks delegated to built-in modes (which cannot update state files directly) and updating `project_overview.json` accordingly.

### Added
- **New Agent:** `strategic-planner` - Focuses solely on initial planning and state setup.
- **New Agent:** `workflow-coordinator` - Focuses solely on execution management, task delegation (including to built-in modes), and state updates during the workflow.
- **Handoff Mechanism:** Introduced an explicit handoff where `strategic-planner` uses `<switch_mode>` to pass control to `workflow-coordinator`.

### Changed
- **Error Handling:** Protocols adapted to handle errors originating from both specialist agents (diagnosed via their state files) and built-in modes (often requiring more user interaction or diagnostics via `debug` mode).
- **Test Execution Flow:** The `workflow-coordinator` now manages the user decision point for test execution after implementation tasks (handled by `code`/`debug`) are inferred to be complete.
- **Specialist Agent Interaction:** `solution-architect`, `ux-specialist`, `guardian-validator`, `docu-crafter` retain their core roles but now interact with the `workflow-coordinator`.

### Removed
- **Agent:** `master-orchestrator` - Responsibilities split between the new `strategic-planner` and `workflow-coordinator`.
- **Agent:** `apex-implementer` - Responsibilities transferred to built-in `code`/`debug` modes managed by the `workflow-coordinator`.

## [v0.0.7] - 2025-04-27

### Changed
- **Agent Directives Updated (Split State):** Updated `customInstructions` across multiple agents to implement a split state management strategy using `project_overview.json` for high-level status and `.state/tasks/{taskId}.json` for detailed task information. If previsouly work with `project_state.json`, we suggest deleting it and using the new `project_overview.json` and `.state/tasks/{taskId}.json` files.
- **State Management Strategy:** Refined the "CRITICAL JSON EDITING STRATEGY" to apply specifically to either the own task file or the overview file depending on the agent's role and the information being updated.
- **Core Optimization Strategy:** Implemented a core strategy for agents to optimize state I/O. This strategy prioritizes the task payload, reads the agent's own task file only once initially, performs selective external reads (overview, docs, specs, design) only when necessary, consolidates multiple updates to the agent's own task file into a single edit operation per logical step or completion, and minimizes updates to the project overview file by only updating the summary status for the agent's task ID upon reaching a final state (Done, Implemented, Error, Failed, Validated).

## [v0.0.6] - 2025-04-26

### Fixed
- **Apex Implementer Behavior:** Fixed an issue where the Apex Implementer mode would repeatedly read files unnecessarily.
- **Safe JSON Editing Strategy:** Added explicit "CRITICAL JSON EDITING STRATEGY" instructions to all agents (`Master Orchestrator`, `Solution Architect`, `UX Specialist`, `Apex Implementer`, `Guardian Validator`, `DocuCrafter`) that modify `project_state.json`. This strategy mandates reading the relevant object, reconstructing the full object with changes, mentally validating syntax, and replacing the entire object via `edit` to minimize syntax errors caused by incremental patching. Includes error handling instructions if edits fail.

### Changed
- **Agent Directives Updated (v5/v6/v7/v8/v9/v10):** Updated `customInstructions` across multiple agents to incorporate the Safe JSON Editing Strategy alongside previous context-aware and workflow changes.

## [v0.0.5] - 2025-04-26

### Added
- **New Task Status:** Added `'Implemented'` status to `project_state.json` task schema to signify code/unit tests written but tests not yet executed, creating a user decision point.

### Changed
- **Workflow Modification:** Introduced an explicit user decision point managed by the Master Orchestrator after a task reaches the `'Implemented'` status.
- **Master Orchestrator Workflow:** Orchestrator now detects `'Implemented'` status, **prompts the user** to decide on test execution (run tests, skip to validation, defer), and delegates the next step (a specific test execution task or a validation task) based on user input. Updated `customInstructions` to v6.
- **Apex Implementer Behavior:** Implementer now writes unit tests but **stops before executing them**, updating task status to `'Implemented'` upon successful code/test writing. Updated `customInstructions` to v8.
- **Test Execution Delegation:** Orchestrator now delegates test execution (when requested by user) as a distinct task (e.g., `type: 'chore'`), often targeting `Guardian Validator` or potentially `Apex Implementer`.
- **Guardian Validator Capability:** Validator instructions updated (v6) to handle specific "run tests" tasks in addition to full validation tasks, reporting pass/fail accordingly.

## [v0.0.4] - 2025-04-26

### Added
- **Task `type` Field:** Introduced a mandatory `type` field (`feature`, `bugfix`, `refactor`, `chore`, `init`, `design`, `validation`, `documentation`) to the `project_state.json` task schema to enable more adaptive workflows.
- **Enhanced Task References:** Added optional `references.reports` and `references.sourceCode` to the task schema.
- **Acceptance Criteria:** Added optional `acceptanceCriteria` array to the task schema.

### Changed
- **Agent Directives Updated (v4/v6):** Updated `customInstructions` for Master Orchestrator (v4), Solution Architect (v4), Apex Implementer (v6), and Guardian Validator (v4) to recognize and adapt behavior based on the new task `type`.
- **Adaptive Orchestration:** Master Orchestrator now delegates tasks based on their `type`, enabling different workflows for features, bug fixes, refactoring, etc.
- **Adaptive Design/Analysis:** Solution Architect now handles analysis for bug fixes/refactors and creates type-specific plans/specs.
- **Adaptive Implementation:** Apex Implementer now adjusts its process based on task type (reading specs, bug reports, or refactoring plans) and emphasizes type-specific testing.
- **Adaptive Validation:** Guardian Validator now tailors its validation strategy based on task type (checking feature criteria, confirming bug fixes, ensuring refactor regressions).
- **Task Status Clarification:** Added `Validated` and `Failed` to the list of possible task statuses in schema documentation for clarity.

## [v0.0.3] - 2025-04-26

### Changed
- **Agent Directives Updated (v3/v4/v5):** Enhanced custom instructions across all agents for better coordination, stricter state management, and refined error handling centered around `project_state.json`.
- **Stricter Schema Enforcement:** Explicitly defined and mandated adherence to the `project_state.json` schema structure for all state interactions.
- **Master Orchestrator:** Added schema reference, enabled concurrent implementation task delegation, and improved reactive error/debug monitoring.
- **Apex Implementer:** Made post-implementation code refinement optional/opportunistic and added structured error details to the state log on failure.
- **Guardian Validator:** **Shifted detailed failure reports** from the state file log to dedicated files within a new `.reports/` directory, linking them from the state log.
- **General:** Minor clarifications on state interaction (`read`/`edit`), task references, and debugging protocols across agents.

## [v0.0.2] - 2025-04-25

### Changed
- **Agent Directive Overhaul (v2/v3/v4):** Updated agent instructions in `.roomodes` to enhance coordination and robustness.
- **Centralized State:** Mandated `project_state.json` for all workflow state, task tracking, and status updates.
- **Standardized Protocols:** Clarified task handoff (Architect defines, Orchestrator delegates based on state), status reporting (agents update `project_state.json`), logging (within `project_state.json`), and `new_task` payload structure.
- **Refined Error Handling:** Restricted Orchestrator's direct fixes and improved the interactive debugging protocol using `Blocked-Debug` status.
- **Increased Clarity:** Improved consistency and definitions across agent roles.

## [v0.0.1] - Initial Release
- Initial project setup and definition of core agent roles.