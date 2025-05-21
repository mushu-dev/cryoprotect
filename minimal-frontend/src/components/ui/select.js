/**
 * Select Component
 * Simplified select components for the minimal frontend
 */
import React from 'react';

export function Select({ children, defaultValue, onValueChange, className = '', ...props }) {
  return (
    <select
      className={`h-10 w-full rounded-md border border-input bg-background px-3 py-2 text-sm ring-offset-background focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring focus-visible:ring-offset-2 ${className}`}
      defaultValue={defaultValue}
      onChange={(e) => onValueChange && onValueChange(e.target.value)}
      {...props}
    >
      {children}
    </select>
  );
}

export function SelectTrigger({ children, className = '', ...props }) {
  return (
    <div
      className={`flex h-10 w-full items-center justify-between rounded-md border border-input bg-background px-3 py-2 text-sm ring-offset-background ${className}`}
      {...props}
    >
      {children}
    </div>
  );
}

export function SelectValue({ children, placeholder, ...props }) {
  return <span {...props}>{children || placeholder}</span>;
}

export function SelectContent({ children, className = '', ...props }) {
  return (
    <div
      className={`relative z-50 min-w-[8rem] overflow-hidden rounded-md border bg-popover text-popover-foreground shadow-md animate-in fade-in-80 ${className}`}
      {...props}
    >
      <div className="p-1">{children}</div>
    </div>
  );
}

export function SelectItem({ children, value, className = '', ...props }) {
  return (
    <option value={value} className={`${className}`} {...props}>
      {children}
    </option>
  );
}