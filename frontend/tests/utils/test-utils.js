import React from 'react';
import { render } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { ThemeProvider } from 'next-themes';

/**
 * Custom render function for components with theme provider
 * @param {React.ReactElement} ui - The component to render
 * @param {Object} options - Additional render options
 * @returns {Object} Rendered component with testing utilities
 */
export function renderWithTheme(ui, { theme = 'light', ...options } = {}) {
  const Wrapper = ({ children }) => (
    <ThemeProvider attribute="class" defaultTheme={theme} enableSystem={false}>
      {children}
    </ThemeProvider>
  );
  
  return render(ui, { wrapper: Wrapper, ...options });
}

/**
 * Setup a component for testing with user events
 * @param {React.ReactElement} ui - The component to render
 * @param {Object} options - Additional render options
 * @returns {Object} Rendered component with testing utilities including user events
 */
export function setup(ui, options = {}) {
  const user = userEvent.setup();
  const utils = renderWithTheme(ui, options);
  
  return {
    ...utils,
    user
  };
}

/**
 * Mocks Next.js router in tests
 * @param {Object} overrides - Router properties to override
 * @returns {Object} Mocked router object
 */
export function mockRouter(overrides = {}) {
  return {
    push: jest.fn(),
    replace: jest.fn(),
    prefetch: jest.fn(),
    back: jest.fn(),
    forward: jest.fn(),
    refresh: jest.fn(),
    pathname: '',
    query: {},
    asPath: '',
    events: {
      on: jest.fn(),
      off: jest.fn(),
    },
    ...overrides
  };
}

/**
 * Async wait utility - useful for waiting for UI updates
 * @param {number} ms - Milliseconds to wait
 * @returns {Promise} Promise that resolves after the specified time
 */
export const waitFor = (ms = 0) => new Promise(resolve => setTimeout(resolve, ms));

/**
 * React 18 compliant act wrapper for asynchronous actions
 * This helps with asynchronous state updates in React 18's concurrent mode
 * @param {Function} callback - Function to execute inside act
 * @returns {Promise} Result of the callback execution
 */
export async function actAsync(callback) {
  const { act } = require('react-dom/test-utils');
  return await act(async () => {
    return await callback();
  });
}