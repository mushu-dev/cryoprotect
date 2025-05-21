import React from 'react';
import { render, screen, fireEvent } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { Button } from './button';
import { setup } from '../../../tests/utils/test-utils';

// Note: We're using the actual Button component here, not a mock,
// to properly test its React 18 features like useTransition and useId

describe('Button Component', () => {
  test('renders with default props', () => {
    render(<Button>Click me</Button>);
    const button = screen.getByRole('button', { name: /click me/i });
    expect(button).toBeInTheDocument();
    expect(button).toHaveClass('bg-primary');
    expect(button).not.toBeDisabled();
  });

  test('applies different variants correctly', () => {
    const { rerender } = render(<Button variant="secondary">Secondary</Button>);
    expect(screen.getByRole('button')).toHaveClass('bg-secondary');

    rerender(<Button variant="destructive">Destructive</Button>);
    expect(screen.getByRole('button')).toHaveClass('bg-destructive');

    rerender(<Button variant="outline">Outline</Button>);
    expect(screen.getByRole('button')).toHaveClass('border');

    rerender(<Button variant="ghost">Ghost</Button>);
    expect(screen.getByRole('button')).toHaveClass('hover:bg-accent');

    rerender(<Button variant="link">Link</Button>);
    expect(screen.getByRole('button')).toHaveClass('text-primary');
  });

  test('applies different sizes correctly', () => {
    const { rerender } = render(<Button size="default">Default</Button>);
    expect(screen.getByRole('button')).toHaveClass('h-10');

    rerender(<Button size="sm">Small</Button>);
    expect(screen.getByRole('button')).toHaveClass('h-8');

    rerender(<Button size="lg">Large</Button>);
    expect(screen.getByRole('button')).toHaveClass('h-11');

    rerender(<Button size="icon">Icon</Button>);
    expect(screen.getByRole('button')).toHaveClass('h-10 w-10');
  });

  test('handles click events', async () => {
    // Using userEvent instead of fireEvent for more realistic interactions
    const handleClick = jest.fn();
    const { user } = setup(<Button onClick={handleClick}>Click me</Button>);
    
    const button = screen.getByRole('button', { name: /click me/i });
    await user.click(button);
    
    expect(handleClick).toHaveBeenCalledTimes(1);
  });

  test('handles disabled state', async () => {
    const handleClick = jest.fn();
    const { user } = setup(
      <Button disabled onClick={handleClick}>
        Disabled button
      </Button>
    );
    
    const button = screen.getByRole('button', { name: /disabled button/i });
    expect(button).toBeDisabled();
    
    await user.click(button);
    expect(handleClick).not.toHaveBeenCalled();
  });

  test('shows loading state when isLoading is true', () => {
    render(
      <Button isLoading loadingText="Loading...">
        Submit
      </Button>
    );
    
    expect(screen.getByRole('button')).toHaveTextContent('Loading...');
    expect(screen.getByRole('button')).toHaveClass('cursor-wait');
    
    // Should have loading spinner
    expect(screen.getByRole('button').querySelector('.animate-spin')).toBeInTheDocument();
  });

  test('applies custom classNames', () => {
    render(<Button className="test-custom-class">Custom</Button>);
    expect(screen.getByRole('button')).toHaveClass('test-custom-class');
  });

  test('passes additional props to button element', () => {
    render(<Button data-testid="test-button" type="submit">Submit</Button>);
    const button = screen.getByTestId('test-button');
    expect(button).toHaveAttribute('type', 'submit');
  });
});