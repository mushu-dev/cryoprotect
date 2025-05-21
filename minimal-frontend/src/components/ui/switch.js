/**
 * Switch Component
 * Toggle switch for boolean values
 */
import React from 'react';

export function Switch({
  className = '',
  id,
  checked = false,
  onCheckedChange,
  ...props
}) {
  return (
    <input
      type="checkbox"
      id={id}
      role="switch"
      checked={checked}
      onChange={(e) => onCheckedChange && onCheckedChange(e.target.checked)}
      className={`appearance-none w-11 h-6 rounded-full bg-muted transition-colors 
        relative before:content-[''] before:inline-block before:size-5 before:rounded-full 
        before:bg-background before:shadow-sm before:transform before:translate-x-0.5 
        before:transition-transform checked:bg-primary checked:before:translate-x-5 
        focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring 
        focus-visible:ring-offset-2 ${className}`}
      {...props}
    />
  );
}