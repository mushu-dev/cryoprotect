# CryoProtect Frontend

A modern web interface for the CryoProtect application, providing tools for analyzing and optimizing cryoprotectant molecules and mixtures.

## Features

- 3D Molecule Visualization with 3DMol.js
- Mixture Composition Analysis with Plotly.js
- Interactive Property Explorer
- User Authentication and Profile Management
- Responsive Design with Dark/Light Mode

## Technology Stack

- **Framework**: Next.js 14 with App Router
- **Language**: TypeScript
- **Styling**: Tailwind CSS with Shadcn UI
- **State Management**:
  - React Query for server state
  - Zustand for client state
- **Authentication**: NextAuth.js
- **Visualization**:
  - 3DMol.js for molecular visualization
  - Plotly.js for data visualization
- **Form Handling**: React Hook Form with Zod validation
- **Deployment**: Vercel

## Getting Started

### Prerequisites

- Node.js 18.x or higher
- npm or yarn

### Installation

1. Clone the repository
```bash
git clone https://github.com/yourusername/cryoprotect-app.git
cd cryoprotect-app/frontend
```

2. Install dependencies
```bash
npm install
# or
yarn install
```

3. Create a `.env.local` file in the project root with the following variables:
```
# API URL for development
NEXT_PUBLIC_API_URL=http://localhost:5000

# Default API version
NEXT_PUBLIC_API_VERSION=v1

# Enable/disable mock data in development
NEXT_PUBLIC_USE_MOCK_DATA=true

# Authentication settings
NEXTAUTH_URL=http://localhost:3000
NEXTAUTH_SECRET=your-secret-key-goes-here
```

4. Start the development server
```bash
npm run dev
# or
yarn dev
```

5. Open [http://localhost:3000](http://localhost:3000) in your browser to see the application.

## Project Structure

```
/frontend
├── public/              # Static assets
├── src/                 # Source code
│   ├── app/             # App router pages
│   │   ├── api/         # API routes
│   │   ├── auth/        # Authentication pages
│   │   ├── mixtures/    # Mixture related pages
│   │   ├── molecules/   # Molecule related pages
│   │   ├── profile/     # User profile page
│   │   ├── settings/    # User settings page
│   │   ├── page.tsx     # Home page
│   │   └── ...          # Other pages
│   ├── components/      # Shared components
│   │   ├── ui/          # UI components (Shadcn UI)
│   │   └── ...          # Other components
│   ├── features/        # Feature-based organization
│   │   ├── auth/        # Authentication features
│   │   ├── dashboard/   # Dashboard features
│   │   ├── mixtures/    # Mixture features
│   │   ├── molecules/   # Molecule features
│   │   └── properties/  # Property explorer features
│   ├── lib/             # Utility functions and libraries
│   └── styles/          # Global styles
├── .env.local           # Environment variables
├── next.config.js       # Next.js configuration
├── package.json         # Dependencies and scripts
├── tailwind.config.js   # Tailwind CSS configuration
└── tsconfig.json        # TypeScript configuration
```

## Key Features and Usage

### Molecule Viewer

The molecule viewer uses 3DMol.js to render interactive 3D visualizations of molecular structures. Features include:

- Multiple visualization styles (stick, sphere, cartoon, surface)
- Rotation and zooming controls
- Optional animation (spinning)
- Molecular information display

### Mixture Composition

The mixture composition component uses Plotly.js to visualize the components of cryoprotectant mixtures:

- Pie charts for composition visualization
- Concentration and ratio analysis
- Interactive tooltips with detailed information
- Component selection for detailed view

### Property Explorer

The property explorer allows users to analyze and visualize molecular properties:

- Interactive charts for property visualization
- Property comparison across multiple molecules
- Correlation analysis between different properties
- Data export functionality

## Authentication and User Management

The application includes a complete authentication system with:

- User registration and login
- Profile management
- Password reset
- Role-based access control
- Protected routes

## Development

### Adding New Components

1. Create a new component in the appropriate directory
2. Import and use Shadcn UI components for consistent styling
3. Use the following pattern for client components:

```tsx
'use client'

import { useState } from 'react'
import { ComponentName } from '@/components/ui/component-name'

export function YourComponent() {
  // Component logic here
  
  return (
    // JSX here
  )
}
```

### API Integration

The application integrates with the CryoProtect API using React Query:

```tsx
import { useQuery } from '@tanstack/react-query'
import { moleculeService } from '@/features/molecules/services/molecule-service'

export function useMolecules(params) {
  return useQuery({
    queryKey: ['molecules', params],
    queryFn: () => moleculeService.getMolecules(params),
    staleTime: 5 * 60 * 1000, // 5 minutes
  })
}
```

### Deployment

The application is configured for deployment on Vercel. To deploy:

1. Push your changes to your repository
2. Connect your repository to Vercel
3. Configure the environment variables
4. Deploy!

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgements

- [Next.js](https://nextjs.org/)
- [Tailwind CSS](https://tailwindcss.com/)
- [Shadcn UI](https://ui.shadcn.com/)
- [3DMol.js](https://3dmol.csb.pitt.edu/)
- [Plotly.js](https://plotly.com/javascript/)
- [React Query](https://tanstack.com/query/latest)
- [NextAuth.js](https://next-auth.js.org/)