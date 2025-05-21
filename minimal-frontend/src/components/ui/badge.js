/**
 * Badge Component
 * Small status indicator or label
 */
import React from 'react';

export function Badge({
  children,
  className = '',
  variant = 'default',
  ...props
}) {
  // Base classes
  const baseClasses = 'inline-flex items-center rounded-full px-2.5 py-0.5 text-xs font-semibold transition-colors';
  
  // Variant classes
  const variantClasses = {
    default: 'bg-primary text-primary-foreground',
    secondary: 'bg-secondary text-secondary-foreground',
    destructive: 'bg-destructive text-destructive-foreground',
    outline: 'text-foreground border border-input bg-background',
  };

  const classes = `${baseClasses} ${variantClasses[variant] || variantClasses.default} ${className}`;

  return (
    <span className={classes} {...props}>
      {children}
    </span>
  );
}