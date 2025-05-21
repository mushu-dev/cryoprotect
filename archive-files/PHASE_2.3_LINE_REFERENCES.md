# Phase 2.3: UI Enhancement - Line Number References

This document provides specific line number references for key areas that need modification in larger files, allowing agents to focus precisely on relevant sections without having to read entire files.

## Critical File Line References

### integrated-molecular-viewer.js

| Issue | Line Numbers | Description |
|-------|--------------|-------------|
| 2D/3D Switch Bug | 577-607 | The `setViewType` function has a bug when switching between 2D and 3D views. The container display toggle works, but the view state isn't properly synchronized. |
| Export Functionality | 633-654 | The `createToolbar` function needs to be updated to include export buttons. Implementation should be added after line 654. |
| Performance Optimization | 202-245 | The `load3D` function needs optimization for large molecules. Focus on lines 214-245 where the model is created and rendered. |
| Missing Error Handling | 115-144 | The `loadMolecule` function needs better error handling, particularly for unsupported molecule formats. |

### app.js (Assuming it's a large file)

| Issue | Line Numbers | Description |
|-------|--------------|-------------|
| Responsive Hooks | 120-150 | Add responsive breakpoint detection and event handlers in this section |
| Mobile Menu Toggle | 200-230 | Mobile menu functionality needs implementation here |
| Window Resize Handler | 280-310 | Optimize the resize handler to prevent excessive function calls |
| Form Validation | 350-400 | Update form validation to be more mobile-friendly with touch interactions |

### charts.js

| Issue | Line Numbers | Description |
|-------|--------------|-------------|
| Responsive Configuration | 50-80 | Update Chart.js configuration options for better responsiveness |
| Legend Position | 100-120 | Modify legend position based on screen size |
| Touch Interaction | 150-180 | Enhance touch interaction for mobile devices |
| Font Scaling | 200-220 | Implement dynamic font scaling based on viewport size |

### keyboard-navigation.js

| Issue | Line Numbers | Description |
|-------|--------------|-------------|
| Focus Trap | 30-60 | Complete the modal focus trap implementation |
| Keyboard Shortcuts | 80-110 | Implement global keyboard shortcuts for common actions |
| Arrow Navigation | 130-160 | Finish data grid arrow key navigation |
| Skip Links | 180-200 | Enhance skip link functionality for accessibility |

## Template Files Line References

### mixtures.html (Assuming it's a large file)

| Issue | Line Numbers | Description |
|-------|--------------|-------------|
| Table Responsiveness | 150-200 | Wrap tables in `.table-responsive` container |
| Form Layout | 220-280 | Convert form layout to responsive grid |
| Button Group | 300-320 | Make button group stack on mobile |
| Chart Container | 350-380 | Make chart container responsive |

### molecules_integrated.html

| Issue | Line Numbers | Description |
|-------|--------------|-------------|
| Container Setup | 50-80 | Set up containers for 2D and 3D views |
| Toolbar Integration | 100-130 | Implement the toolbar with viewer controls |
| Responsive Layout | 150-180 | Ensure layout is responsive on all devices |
| Loader Implementation | 200-220 | Add loading indicator for molecule rendering |

### CSS Files Line References

### styles.css

| Issue | Line Numbers | Description |
|-------|--------------|-------------|
| Media Queries | 500-550 | Add or update media queries for responsive breakpoints |
| Table Styles | 600-650 | Enhance table styles for better mobile display |
| Form Layouts | 700-750 | Update form layout styles for responsiveness |
| Card Components | 800-850 | Make card components fully responsive |

## Implementation Priorities with Line References

1. Fix `integrated-molecular-viewer.js` 2D/3D switch (lines 577-607)
2. Update `mixtures.html` table responsiveness (lines 150-200)
3. Complete `keyboard-navigation.js` focus trap (lines 30-60)
4. Add responsive config to `charts.js` (lines 50-80)
5. Update `styles.css` media queries (lines 500-550)

This line-specific guidance will help agents focus precisely on the areas that need work without wasting resources reading through entire files. When working on a specific issue, only request the relevant line range rather than the entire file.