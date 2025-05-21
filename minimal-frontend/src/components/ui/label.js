/**
 * Label Component
 * Accessible label for form controls
 */
import React from 'react';

export function Label({
  className = '',
  htmlFor,
  children,
  ...props
}) {
  return (
    <label
      htmlFor={htmlFor}
      className={`text-sm font-medium leading-none peer-disabled:cursor-not-allowed peer-disabled:opacity-70 ${className}`}
      {...props}
    >
      {children}
    </label>
  );
}