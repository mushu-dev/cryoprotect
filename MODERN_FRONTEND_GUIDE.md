# Modern Frontend Guide for CryoProtect

This guide details the modern frontend implementation for CryoProtect using Next.js and other cutting-edge technologies.

## Overview

The new CryoProtect frontend is built with Next.js and deployed on Vercel, providing a fast, responsive, and accessible user interface for molecular data analysis. This document outlines the architecture, technology choices, and implementation details.

## Technology Stack

- **Framework**: Next.js 14+ with App Router
- **Language**: TypeScript 5.0+
- **Styling**: Tailwind CSS with Shadcn UI components
- **State Management**:
  - React Query for server state
  - Zustand for client state
- **Visualization Libraries**:
  - 3DMol.js for molecular visualization
  - Plotly.js for data visualization
  - D3.js for custom visualizations
- **Form Management**: React Hook Form with Zod validation
- **Authentication**: NextAuth.js (supports Supabase integration)
- **Testing**:
  - Vitest for unit and component testing
  - Playwright for end-to-end testing

## Features

- **Responsive Design**: Works on desktop, tablet, and mobile devices
- **Dark/Light Mode**: Complete theme support for user preference
- **Accessible UI**: ARIA compliance and keyboard navigation
- **Interactive Visualizations**:
  - 3D molecular structure rendering
  - Interactive charts for property analysis
  - Composition visualization for mixtures
- **Advanced Data Filtering**: Complex scientific data exploration
- **Real-time Data Updates**: Automatic refreshing of data
- **Authentication**: Secure user login and registration
- **Cross-browser Compatibility**: Works in all modern browsers

## Directory Structure

The frontend follows a feature-based architecture:

```
frontend/
├── src/
│   ├── app/               # Next.js App Router
│   ├── components/        # Shared UI components
│   ├── features/          # Feature modules
│   │   ├── molecules/     # Molecule-related features
│   │   ├── mixtures/      # Mixture-related features
│   │   ├── dashboard/     # Dashboard features
│   │   └── auth/          # Authentication features
│   ├── hooks/             # Shared custom hooks
│   ├── lib/               # Library code
│   ├── styles/            # Global styles
│   ├── types/             # TypeScript type definitions
│   └── utils/             # Utility functions
├── public/                # Static assets
```

## Feature Modules

Each feature module includes:

- **Components**: UI components specific to the feature
- **Hooks**: Custom hooks for feature-specific logic
- **Utils**: Utility functions for the feature
- **Types**: TypeScript types for the feature
- **API**: API integration code

## Getting Started

### Setup Development Environment

1. Install dependencies:
   ```bash
   cd frontend
   npm install
   ```

2. Start the development server:
   ```bash
   npm run dev
   ```

3. Open [http://localhost:3000](http://localhost:3000) in your browser.

### Building for Production

```bash
npm run build
```

### Running Tests

```bash
# Unit tests
npm run test

# End-to-end tests
npm run test:e2e
```

## Deployment on Vercel

The frontend is designed for seamless deployment on Vercel:

1. Connect your GitHub repository to Vercel
2. Configure environment variables:
   - `NEXT_PUBLIC_API_URL`
   - `NEXT_PUBLIC_SUPABASE_URL` (if applicable)
   - `NEXT_PUBLIC_SUPABASE_ANON_KEY` (if applicable)
3. Deploy with the default Next.js settings

## API Integration

The frontend integrates with the CryoProtect API through a set of typed client functions in the `src/lib/api.ts` file. These functions use the Axios library for HTTP requests and React Query for caching and state management.

Example API integration:

```typescript
// src/features/molecules/api.ts
import axios from 'axios';
import { Molecule } from '@/types';

const API_URL = process.env.NEXT_PUBLIC_API_URL;

export async function getMolecules(params?: {
  limit?: number;
  offset?: number;
  sort?: string;
  order?: 'asc' | 'desc';
}) {
  const response = await axios.get(`${API_URL}/molecules`, { params });
  return response.data.data;
}

export async function getMolecule(id: string) {
  const response = await axios.get(`${API_URL}/molecules/${id}`);
  return response.data.data;
}
```

## UI Components

The frontend uses a combination of custom components and Shadcn UI components. Key components include:

- **MoleculeViewer3D**: 3D visualization of molecules using 3DMol.js
- **MixtureCompositionChart**: Visualization of mixture compositions using Plotly.js
- **PropertyExplorer**: Advanced filtering and visualization of molecular properties
- **ConsolidatedMoleculeManager**: Interface for managing consolidated molecules

## Theming

The application uses Tailwind CSS with a custom theme configuration that supports both light and dark modes. The theme is defined in `tailwind.config.js` and includes:

- **Colors**: Primary, secondary, accent, and other color variables
- **Typography**: Font families, sizes, and weights
- **Spacing**: Margins, paddings, and gaps
- **Breakpoints**: Responsive design breakpoints

## Performance Optimizations

Several performance optimizations are implemented:

1. **Code Splitting**: Dynamic imports for large components
2. **Image Optimization**: Next.js Image component for optimized images
3. **Memoization**: React.memo, useMemo, and useCallback for expensive operations
4. **Virtualization**: For large lists of molecules
5. **Incremental Static Regeneration**: For static pages with dynamic data
6. **API Response Caching**: With React Query's caching strategies

## Accessibility

Accessibility features include:

1. **Keyboard Navigation**: Full keyboard support for all interactive elements
2. **Screen Reader Support**: ARIA labels and roles
3. **Focus Management**: Visible focus indicators
4. **Color Contrast**: WCAG 2.1 compliant color contrast
5. **Reduced Motion**: Support for reduced motion preferences

## Future Enhancements

Planned enhancements for the frontend include:

1. **Advanced Search**: Complex query builder for molecular searches
2. **Collaborative Features**: Real-time collaboration on mixtures and experiments
3. **Machine Learning Integration**: Visual interface for ML-based property predictions
4. **Offline Support**: Progressive Web App features for offline use
5. **Export Formats**: Additional export formats for visualization and data

## Conclusion

The modern frontend for CryoProtect provides a robust, scalable, and user-friendly interface for molecular data analysis. By leveraging modern web technologies and best practices, it offers a significant improvement over traditional interfaces for scientific applications.