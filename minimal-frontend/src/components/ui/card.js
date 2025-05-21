/**
 * Card Component
 * Simple card component for displaying content in a contained card format
 */
import React from 'react';

export function Card({ children, className = '', ...props }) {
  return (
    <div 
      className={`bg-white rounded-lg border shadow-sm overflow-hidden ${className}`} 
      {...props}
    >
      {children}
    </div>
  );
}

export function CardHeader({ children, className = '', ...props }) {
  return (
    <div 
      className={`p-4 border-b ${className}`} 
      {...props}
    >
      {children}
    </div>
  );
}

export function CardTitle({ children, className = '', ...props }) {
  return (
    <h3 
      className={`font-semibold text-lg ${className}`} 
      {...props}
    >
      {children}
    </h3>
  );
}

export function CardContent({ children, className = '', ...props }) {
  return (
    <div 
      className={`p-4 ${className}`} 
      {...props}
    >
      {children}
    </div>
  );
}

export function CardFooter({ children, className = '', ...props }) {
  return (
    <div 
      className={`p-4 border-t ${className}`} 
      {...props}
    >
      {children}
    </div>
  );
}