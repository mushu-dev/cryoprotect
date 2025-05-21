/**
 * Pagination Component
 * Reusable pagination component for lists and tables
 */
import React from 'react';

export function Pagination({ className = '', ...props }) {
  return (
    <nav
      role="navigation"
      aria-label="pagination"
      className={`mx-auto flex w-full justify-center ${className}`}
      {...props}
    />
  );
}

export function PaginationContent({ className = '', ...props }) {
  return (
    <ul className={`flex flex-row items-center gap-1 ${className}`} {...props} />
  );
}

export function PaginationItem({ className = '', ...props }) {
  return <li className={className} {...props} />;
}

export function PaginationLink({
  className = '',
  isActive = false,
  children,
  onClick,
  ...props
}) {
  return (
    <button
      onClick={onClick}
      aria-current={isActive ? 'page' : undefined}
      className={`flex h-10 w-10 items-center justify-center rounded-md border ${
        isActive
          ? 'bg-primary text-primary-foreground'
          : 'bg-background text-foreground hover:bg-muted hover:text-muted-foreground'
      } ${className}`}
      {...props}
    >
      {children}
    </button>
  );
}

export function PaginationPrevious({ className = '', onClick, ...props }) {
  return (
    <PaginationLink
      onClick={onClick}
      aria-label="Go to previous page"
      className={`gap-1 pl-2.5 ${className}`}
      {...props}
    >
      <span>Previous</span>
    </PaginationLink>
  );
}

export function PaginationNext({ className = '', onClick, ...props }) {
  return (
    <PaginationLink
      onClick={onClick}
      aria-label="Go to next page"
      className={`gap-1 pr-2.5 ${className}`}
      {...props}
    >
      <span>Next</span>
    </PaginationLink>
  );
}

export function PaginationEllipsis({ className = '', ...props }) {
  return (
    <span
      aria-hidden
      className={`flex h-10 w-10 items-center justify-center ${className}`}
      {...props}
    >
      <span>...</span>
    </span>
  );
}