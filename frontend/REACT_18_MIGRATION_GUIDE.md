# React 18 Migration Guide for CryoProtect Frontend

This guide outlines the steps required to migrate the CryoProtect frontend from React 17 to React 18 and update Next.js from version 12.3.4 to 14.x.

## Prerequisites

Before starting the migration:

1. Create a new branch for the migration
2. Make sure all tests pass on the current version
3. Have a good understanding of the codebase structure

## Step 1: Update Package Dependencies

Update package.json with the latest React and Next.js versions:

```bash
# Save the updated package.json
cp package.json.updated package.json

# Install dependencies
npm install
```

## Step 2: Update React Root Rendering

React 18 introduces a new root API for rendering. Update each entry point:

### For Next.js App Router (New)

Create `src/app/layout.tsx` if not already present:

```tsx
import React from 'react';
import { Metadata } from 'next';
import '../styles/globals.css';
import { ThemeProvider } from '../components/theme-provider';
import NavigationHeader from '../components/navigation-header-updated';
import Footer from '../components/footer';
import { ConvexClientProvider } from '../convex/ConvexClientProvider';

export const metadata: Metadata = {
  title: 'CryoProtect',
  description: 'A platform for cryoprotectant analysis',
};

export default function RootLayout({
  children,
}: {
  children: React.ReactNode;
}) {
  return (
    <html lang="en" suppressHydrationWarning>
      <body>
        <ConvexClientProvider>
          <ThemeProvider 
            attribute="class" 
            defaultTheme="system" 
            enableSystem
            disableTransitionOnChange
          >
            <div className="flex min-h-screen flex-col">
              <NavigationHeader />
              <main className="flex-grow">
                {children}
              </main>
              <Footer />
            </div>
          </ThemeProvider>
        </ConvexClientProvider>
      </body>
    </html>
  );
}
```

### For Pages Router (Existing)

Update `src/pages/_app.js` to use the new React 18 API:

```tsx
import { useState, useEffect } from 'react'
import { ThemeProvider } from 'next-themes'
import '../styles/globals.css'
import NavigationHeader from '../components/navigation-header-updated'
import Footer from '../components/footer'
import { ConvexClientProvider } from '../convex/ConvexClientProvider'

function MyApp({ Component, pageProps }) {
  // Add this if you're using client-side features that need to be protected from SSR
  const [mounted, setMounted] = useState(false)
  
  useEffect(() => {
    setMounted(true)
  }, [])

  if (!mounted) {
    return null
  }

  return (
    <ConvexClientProvider>
      <ThemeProvider attribute="class" defaultTheme="system" enableSystem>
        <div className="flex min-h-screen flex-col">
          <NavigationHeader />
          <main className="flex-grow">
            <Component {...pageProps} />
          </main>
          <Footer />
        </div>
      </ThemeProvider>
    </ConvexClientProvider>
  )
}

export default MyApp
```

## Step 3: Update Next.js Configuration

Update your Next.js configuration to work with the new version:

```js
// next.config.js
/** @type {import('next').NextConfig} */
const nextConfig = {
  reactStrictMode: true,
  swcMinify: true,
  images: {
    unoptimized: true,
  },
  eslint: {
    ignoreDuringBuilds: true
  },
  typescript: {
    ignoreBuildErrors: true
  },
  // For transitioning, support both pages and app router
  experimental: {
    appDir: true,
  },
  // Static generation configuration
  output: 'export',
  trailingSlash: true,
};

module.exports = nextConfig;
```

## Step 4: Update React API Usage

### Update Event Pooling

React 17 removed the event pooling feature, so you can remove any event.persist() calls:

```tsx
// Before
function handleChange(event) {
  event.persist();
  setFormState(formState => ({
    ...formState,
    [event.target.name]: event.target.value
  }));
}

// After
function handleChange(event) {
  setFormState(formState => ({
    ...formState,
    [event.target.name]: event.target.value
  }));
}
```

### Update useEffect Cleanup

Make sure all useEffect hooks with cleanup functions are properly implemented:

```tsx
// Correct implementation
useEffect(() => {
  const subscription = subscribeToData(props.id);
  
  // Return cleanup function
  return () => {
    subscription.unsubscribe();
  };
}, [props.id]);
```

### Automatic Batching

React 18 implements automatic batching for all state updates. Check for places where state updates might behave differently:

```tsx
// Before (in React 17, these would trigger 2 renders)
function handleClick() {
  setCount(c => c + 1);
  setFlag(f => !f);
}

// After (in React 18, these will be batched into a single render)
function handleClick() {
  setCount(c => c + 1);
  setFlag(f => !f);
}
```

## Step 5: Update Libraries and Components

### Update Convex Integration

Update the Convex client provider:

```tsx
// src/convex/ConvexClientProvider.tsx
import { ReactNode, useState, useEffect } from 'react';
import { ConvexReactClient } from 'convex/react';
import { ConvexProviderWithAuth0 } from 'convex/react-auth0';
import { Auth0Provider } from '@auth0/auth0-react';

const convex = new ConvexReactClient(process.env.NEXT_PUBLIC_CONVEX_URL || '');

export function ConvexClientProvider({ children }: { children: ReactNode }) {
  const [mounted, setMounted] = useState(false);
  
  useEffect(() => {
    setMounted(true);
    return () => setMounted(false);
  }, []);
  
  // Don't render anything on the server
  if (!mounted && typeof window !== 'undefined') {
    return null;
  }
  
  return (
    <Auth0Provider
      domain={process.env.NEXT_PUBLIC_AUTH0_DOMAIN || ''}
      clientId={process.env.NEXT_PUBLIC_AUTH0_CLIENT_ID || ''}
      authorizationParams={{
        redirect_uri: typeof window !== 'undefined' ? window.location.origin : '',
      }}
    >
      <ConvexProviderWithAuth0 client={convex}>
        {children}
      </ConvexProviderWithAuth0>
    </Auth0Provider>
  );
}
```

### Update Third-Party Libraries

Some libraries might need updates to work with React 18:

```bash
# Update form libraries
npm install react-hook-form@latest @hookform/resolvers@latest

# Update visualization libraries
npm install @react-three/fiber@latest @react-three/drei@latest three@latest
```

## Step 6: Handle React 18 Strict Mode

React 18's Strict Mode now double-invokes effect functions to help find side effects. Update any components that might be affected:

```tsx
function ExampleComponent() {
  const [data, setData] = useState(null);
  
  useEffect(() => {
    // This will run twice in development mode with React 18 Strict Mode
    let isMounted = true;
    fetchData().then(result => {
      // Check if component is still mounted before updating state
      if (isMounted) {
        setData(result);
      }
    });
    
    return () => {
      isMounted = false;
    };
  }, []);
  
  return <div>{/* component content */}</div>;
}
```

## Step 7: Add TypeScript Types

For better type safety, add proper TypeScript types:

```tsx
// Example of adding types to component props
interface ButtonProps extends React.ButtonHTMLAttributes<HTMLButtonElement> {
  variant?: 'default' | 'destructive' | 'outline' | 'secondary' | 'ghost' | 'link';
  size?: 'default' | 'sm' | 'lg' | 'icon';
  asChild?: boolean;
}

// Example of adding types to state variables
const [isLoading, setIsLoading] = useState<boolean>(false);
const [data, setData] = useState<ApiResponse | null>(null);
```

## Step 8: Use New React 18 Features

Take advantage of new React 18 features:

### Transitions API

```tsx
import { useTransition, useState } from 'react';

function TabContainer() {
  const [isPending, startTransition] = useTransition();
  const [tab, setTab] = useState('home');
  
  function selectTab(nextTab) {
    startTransition(() => {
      setTab(nextTab);
    });
  }
  
  return (
    <div>
      <TabButton
        isActive={tab === 'home'}
        onClick={() => selectTab('home')}
      >
        Home
      </TabButton>
      <TabButton
        isActive={tab === 'about'}
        onClick={() => selectTab('about')}
      >
        About
      </TabButton>
      {isPending ? (
        <Spinner />
      ) : (
        <TabPanel tab={tab} />
      )}
    </div>
  );
}
```

### Suspense for Data Fetching

```tsx
import { Suspense } from 'react';

function MoleculeViewer({ id }) {
  return (
    <div>
      <h2>Molecule Structure</h2>
      <Suspense fallback={<Spinner />}>
        <MoleculeData id={id} />
      </Suspense>
    </div>
  );
}
```

## Step 9: Testing

1. Run the existing test suite to check for regressions:

```bash
npm run test:all
```

2. Create new tests for React 18-specific features:

```tsx
// tests/react18/transitions.test.js
import { render, screen, fireEvent } from '@testing-library/react';
import { TabContainer } from '../components/TabContainer';

test('should show spinner during tab transition', async () => {
  render(<TabContainer />);
  
  // Click on About tab
  fireEvent.click(screen.getByText('About'));
  
  // Should show spinner briefly
  expect(screen.getByTestId('spinner')).toBeInTheDocument();
  
  // Should eventually show About content
  await screen.findByText('About Content');
});
```

## Step 10: Progressive Rollout

1. First deploy to a staging environment
2. Monitor for unexpected behavior
3. Fix any issues that arise
4. Roll out to production

## Common Issues and Solutions

### SSR Hydration Issues

React 18 is stricter about hydration mismatches. Fix by ensuring server and client render the same content:

```tsx
// Problem: Date is rendered differently on server vs client
function Example() {
  return <div>{new Date().toString()}</div>
}

// Solution: Use useEffect to update client-only content
function Example() {
  const [dateString, setDateString] = useState('');
  
  useEffect(() => {
    setDateString(new Date().toString());
  }, []);
  
  return <div>{dateString}</div>
}
```

### State Updates in useEffect

Double-invocation of effects can cause unexpected state updates:

```tsx
// Problem: Counter increases twice in dev mode
function Example() {
  const [count, setCount] = useState(0);
  
  useEffect(() => {
    setCount(count + 1);
  }, []);
  
  return <div>{count}</div>
}

// Solution: Use functional updates
function Example() {
  const [count, setCount] = useState(0);
  
  useEffect(() => {
    setCount(c => c + 1);
  }, []);
  
  return <div>{count}</div>
}
```

### Concurrent Rendering Issues

Concurrent rendering can reveal race conditions in your code:

```tsx
// Problem: Race condition with fetching data
useEffect(() => {
  fetchData(props.id).then(data => {
    setData(data);
  });
}, [props.id]);

// Solution: Use cleanup function and track component mount state
useEffect(() => {
  let isMounted = true;
  fetchData(props.id).then(data => {
    if (isMounted) {
      setData(data);
    }
  });
  
  return () => {
    isMounted = false;
  };
}, [props.id]);
```

## Conclusion

Migrating to React 18 provides significant benefits in terms of performance, features, and developer experience. By following this guide and carefully testing each step, you can successfully upgrade the CryoProtect frontend to the latest React and Next.js versions.