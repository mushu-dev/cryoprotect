import React from 'react';
import { render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { ThemeProvider } from './theme-provider';
import { setup, actAsync } from '../../tests/utils/test-utils';

// Create a test component to use inside the ThemeProvider
const TestComponent = () => {
  const [theme, setTheme] = React.useState('light');
  
  return (
    <div>
      <span data-testid="current-theme">Current theme: {theme}</span>
      <button 
        data-testid="toggle-theme" 
        onClick={() => setTheme(theme === 'light' ? 'dark' : 'light')}
      >
        Toggle Theme
      </button>
    </div>
  );
};

describe('ThemeProvider Component', () => {
  beforeEach(() => {
    // Reset DOM attributes before each test
    document.documentElement.classList.remove('dark');
    document.documentElement.removeAttribute('data-theme');
    document.documentElement.removeAttribute('style');
    
    // Reset localStorage mock
    localStorage.clear();
  });
  
  test('renders children without crashing', () => {
    render(
      <ThemeProvider>
        <div data-testid="test-child">Test Child</div>
      </ThemeProvider>
    );
    
    expect(screen.getByTestId('test-child')).toBeInTheDocument();
  });
  
  test('provides theme context to children', async () => {
    // Setup with enableSystem=false to ensure deterministic theme
    const { user } = setup(
      <ThemeProvider attribute="class" defaultTheme="light" enableSystem={false}>
        <TestComponent />
      </ThemeProvider>
    );
    
    const themeText = screen.getByTestId('current-theme');
    expect(themeText).toHaveTextContent('Current theme: light');
    
    // Simulate toggling the theme using userEvent
    const toggleButton = screen.getByTestId('toggle-theme');
    await user.click(toggleButton);
    
    // After toggling, the theme should be "dark"
    expect(themeText).toHaveTextContent('Current theme: dark');
    
    // Check that "dark" class was added to the document element (React 18 feature)
    expect(document.documentElement.classList.contains('dark')).toBe(true);
  });
  
  test('respects defaultTheme prop', () => {
    render(
      <ThemeProvider attribute="class" defaultTheme="dark" enableSystem={false}>
        <TestComponent />
      </ThemeProvider>
    );
    
    expect(screen.getByTestId('current-theme')).toHaveTextContent('Current theme: dark');
    expect(document.documentElement.classList.contains('dark')).toBe(true);
  });
  
  test('handles forced color scheme', async () => {
    // This tests the React 18 feature for handling forced color schemes
    await actAsync(async () => {
      render(
        <ThemeProvider 
          attribute="class" 
          defaultTheme="light" 
          enableSystem={false}
          forcedTheme="dark"
        >
          <TestComponent />
        </ThemeProvider>
      );
    });
    
    // Theme should be forced to dark regardless of the default
    expect(document.documentElement.classList.contains('dark')).toBe(true);
    
    // Even after clicking, the theme should remain dark due to forcedTheme
    const toggleButton = screen.getByTestId('toggle-theme');
    await actAsync(async () => {
      await userEvent.setup().click(toggleButton);
    });
    
    expect(document.documentElement.classList.contains('dark')).toBe(true);
  });
  
  test('uses system preference when enableSystem is true', () => {
    // Mock window.matchMedia to simulate system dark mode preference
    Object.defineProperty(window, 'matchMedia', {
      writable: true,
      value: jest.fn().mockImplementation(query => ({
        matches: query === '(prefers-color-scheme: dark)' ? true : false,
        media: query,
        onchange: null,
        addEventListener: jest.fn(),
        removeEventListener: jest.fn(),
        dispatchEvent: jest.fn(),
      })),
    });
    
    render(
      <ThemeProvider attribute="class" enableSystem={true}>
        <TestComponent />
      </ThemeProvider>
    );
    
    // Should detect system preference for dark theme
    expect(document.documentElement.classList.contains('dark')).toBe(true);
  });
  
  test('persists theme selection in localStorage', async () => {
    // Spy on localStorage
    const setItemSpy = jest.spyOn(Storage.prototype, 'setItem');
    
    const { user } = setup(
      <ThemeProvider attribute="class" defaultTheme="light" enableSystem={false}>
        <TestComponent />
      </ThemeProvider>
    );
    
    // Toggle to dark theme
    const toggleButton = screen.getByTestId('toggle-theme');
    await user.click(toggleButton);
    
    // Check localStorage was called with the correct values
    expect(setItemSpy).toHaveBeenCalledWith('theme', 'dark');
    
    // Cleanup
    setItemSpy.mockRestore();
  });
});