// Import jest-dom's custom assertions
require('@testing-library/jest-dom');
const React = require('react');
const { configure } = require('@testing-library/react');

// React 18 testing library configuration
configure({ 
  testIdAttribute: 'data-testid',
  asyncUtilTimeout: 5000 
});

// React 18 - Mock createRoot from react-dom/client
jest.mock('react-dom/client', () => ({
  createRoot: jest.fn((container) => ({
    render: jest.fn((element) => {
      // Simulate React 18 rendering
      return element;
    }),
    unmount: jest.fn()
  }))
}));

// Mock useId hook which is new in React 18
jest.mock('react', () => ({
  ...jest.requireActual('react'),
  useId: () => 'test-id',
}));

// Setup React 18 concurrent mode globals
global.IS_REACT_ACT_ENVIRONMENT = true;

// Mock Next.js components and hooks - App Router style (Next.js 13+)
jest.mock('next/navigation', () => ({
  useRouter() {
    return {
      push: jest.fn(),
      replace: jest.fn(),
      prefetch: jest.fn(),
      back: jest.fn(),
      forward: jest.fn(),
      refresh: jest.fn(),
      pathname: '',
      query: {},
    };
  },
  usePathname() {
    return '';
  },
  useSearchParams() {
    return new URLSearchParams();
  },
  useParams() {
    return {};
  },
}));

// Keep legacy router mock for backwards compatibility
jest.mock('next/router', () => ({
  useRouter() {
    return {
      push: jest.fn(),
      replace: jest.fn(),
      prefetch: jest.fn(),
      back: jest.fn(),
      query: {},
      pathname: '',
      asPath: '',
      events: {
        on: jest.fn(),
        off: jest.fn(),
      },
    };
  },
}));

// Mock next/image component
jest.mock('next/image', () => ({
  __esModule: true,
  default: (props) => {
    // eslint-disable-next-line jsx-a11y/alt-text
    return <img {...props} data-testid="next-image" />; 
  },
}));

// Mock next/link component
jest.mock('next/link', () => ({
  __esModule: true,
  default: ({ children, href, ...props }) => (
    <a href={href} {...props} data-testid="next-link">
      {children}
    </a>
  ),
}));

// Mock next/head component
jest.mock('next/head', () => ({
  __esModule: true,
  default: ({ children }) => <>{children}</>,
}));

// Mock for components that use custom UI elements
jest.mock('@/components/ui/select', () => ({
  Select: (props) => <div className="Select-mock" data-testid="select" {...props} />,
  SelectContent: (props) => <div className="SelectContent-mock" data-testid="select-content" {...props} />,
  SelectItem: (props) => <div className="SelectItem-mock" data-testid="select-item" {...props} />,
  SelectTrigger: (props) => <div className="SelectTrigger-mock" data-testid="select-trigger" {...props} />,
  SelectValue: (props) => <div className="SelectValue-mock" data-testid="select-value" {...props} />,
}));

jest.mock('@/components/ui/button', () => ({
  Button: ({ isLoading, loadingText, children, ...props }) => (
    <button className="Button-mock" data-testid="button" {...props}>
      {isLoading ? loadingText || 'Loading...' : children}
    </button>
  ),
}));

jest.mock('@/components/ui/card', () => ({
  Card: (props) => <div className="Card-mock" data-testid="card" {...props} />,
  CardHeader: (props) => <div className="CardHeader-mock" data-testid="card-header" {...props} />,
  CardTitle: (props) => <div className="CardTitle-mock" data-testid="card-title" {...props} />,
  CardDescription: (props) => <div className="CardDescription-mock" data-testid="card-description" {...props} />,
  CardContent: (props) => <div className="CardContent-mock" data-testid="card-content" {...props} />,
  CardFooter: (props) => <div className="CardFooter-mock" data-testid="card-footer" {...props} />,
}));

// Mock Lucide React icons
const mockIcon = (name) => (props) => <div className={`${name}-mock`} data-testid={`icon-${name.toLowerCase()}`} {...props} />;

jest.mock('lucide-react', () => {
  const icons = [
    'RotateCcw', 'RotateCw', 'ZoomIn', 'ZoomOut', 'Download', 'Maximize', 'Layers',
    'Check', 'X', 'ChevronLeft', 'ChevronRight', 'ChevronUp', 'ChevronDown',
    'Search', 'Plus', 'Minus', 'Settings', 'User', 'Menu', 'Edit', 'Trash',
    'AlertCircle', 'Info', 'Moon', 'Sun', 'Bell', 'Calendar', 'Clock', 'Home',
    'Filter', 'SortAsc', 'SortDesc', 'Eye', 'EyeOff', 'Heart', 'Mail',
    'Database', 'Beaker', 'Flask', 'TestTube', 'Microscope', 'BarChart', 'ServerCog'
  ];
  
  const iconMocks = {};
  icons.forEach(icon => {
    iconMocks[icon] = mockIcon(icon);
  });
  
  return iconMocks;
});

// Conditionally silent console errors during tests based on environment
// In CI we want to see errors
if (process.env.NODE_ENV === 'test' && !process.env.CI) {
  console.error = jest.fn();
  console.warn = jest.fn();
}