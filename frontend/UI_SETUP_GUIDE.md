# CryoProtect UI Setup Guide

This document provides a comprehensive guide to set up and use the modern UI components in the CryoProtect application. The UI system is designed to provide a consistent, accessible, and visually appealing interface across the application.

## Table of Contents

1. [Getting Started](#getting-started)
2. [Component Library](#component-library)
3. [Layout System](#layout-system)
4. [Theme Customization](#theme-customization)
5. [Best Practices](#best-practices)
6. [Troubleshooting](#troubleshooting)

## Getting Started

### Prerequisites

Make sure you have the following:

- Node.js (v12+)
- npm or yarn

### Installation

The UI components are already included in the project. If you need to install or update dependencies, run:

```bash
# Using npm
npm install

# Using yarn
yarn
```

### Running the Component Library Documentation

To see all available UI components and their usage:

```bash
npm run dev
```

Then navigate to `/components` in your browser.

## Component Library

The UI system includes the following components:

### Basic Components

- **Button**: For interactive actions
- **Card**: Container for content
- **Badge**: Small status indicators
- **Alert**: Information and notification messages
- **Skeleton**: Loading states
- **Avatar**: User or entity representation
- **Separator**: Visual dividers

### Form Components

- **Input**: Text input fields
- **Textarea**: Multi-line text input
- **Select**: Dropdown selection
- **Switch**: Toggle controls
- **Label**: Form field labels

### Layout Components

- **Sheet**: Slide-in panels
- **Dialog**: Modal dialogs
- **Tabs**: Content organization
- **ToastProvider**: Notification system

### Navigation

- **Navigation Header**: Top navigation with responsive mobile view
- **Dashboard Layout**: Full application layout with sidebar

## Layout System

### Dashboard Layout

The dashboard layout provides a consistent structure for pages with:

- Sidebar navigation
- Top header bar
- Content area
- Responsive behavior for mobile devices

Usage:

```jsx
import DashboardLayout from '../components/layouts/dashboard-layout';

export default function YourPage() {
  return (
    <DashboardLayout>
      {/* Your page content */}
    </DashboardLayout>
  );
}
```

For pages that should automatically use the dashboard layout, add them to the `dashboardPages` array in `_app.js`.

## Theme Customization

Colors and styling are controlled through Tailwind CSS variables defined in `globals.css`. Key variables include:

- Primary colors (`--primary`)
- Secondary colors (`--secondary`)
- Background colors (`--background`)
- Text colors (`--foreground`)
- UI element colors (card, border, etc.)

To adjust the theme, modify these variables in `globals.css`.

## Best Practices

### Component Usage

1. Use the appropriate component for each UI element
2. Maintain consistent spacing using the Tailwind spacing system
3. Use the card component for grouping related content
4. Keep the UI clean and minimal, focusing on functionality

### Responsive Design

1. Use the responsive utilities provided by Tailwind (`sm:`, `md:`, `lg:` prefixes)
2. Test on various screen sizes
3. Ensure content is accessible on mobile devices

### Accessibility

1. Include proper ARIA labels
2. Ensure sufficient color contrast
3. Test keyboard navigation
4. Support screen readers

## Troubleshooting

### Common Issues

**Styling Inconsistencies**
- Check that you're using the proper component and utility classes
- Verify that global CSS is being applied

**Component Errors**
- Check import paths
- Verify that all required props are provided
- Ensure proper nesting of related components

**Layout Problems**
- Use browser developer tools to inspect element positioning
- Check responsive breakpoints
- Verify CSS variables are properly defined

### Getting Help

If you encounter issues not covered in this guide, check:

1. The component documentation page (`/components`)
2. The source code for each component for detailed implementation
3. Tailwind CSS documentation for utility classes

For specific project-related questions, refer to the README or contact the development team.