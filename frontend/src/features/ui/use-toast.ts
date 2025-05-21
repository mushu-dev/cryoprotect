/**
 * Hook for displaying toast notifications
 * 
 * Provides a standardized way to show notifications to users.
 * This is a simple implementation that can be replaced with any
 * toast notification library.
 */

import { useCallback } from 'react';

// Types for toast notifications
export type ToastVariant = 'default' | 'destructive' | 'success' | 'warning' | 'info';

export interface ToastOptions {
  title: string;
  description?: string;
  variant?: ToastVariant;
  duration?: number;
  action?: {
    label: string;
    onClick: () => void;
  };
}

/**
 * Hook for displaying toast notifications
 */
export function useToast() {
  /**
   * Show a toast notification
   */
  const toast = useCallback((options: ToastOptions) => {
    // This is a simple implementation that can be replaced with
    // any toast notification library like react-hot-toast,
    // react-toastify, or shadcn/ui toast, etc.
    
    // Default values
    const variant = options.variant || 'default';
    const duration = options.duration || 5000;
    
    // Log to console in development
    if (process.env.NODE_ENV === 'development') {
      console.info(`[TOAST/${variant}] ${options.title}`, options.description);
    }
    
    // TODO: Replace this with your actual toast implementation
    // For now, we'll display toasts using a simple DOM-based approach
    
    // Create toast container if it doesn't exist
    let container = document.getElementById('toast-container');
    if (!container) {
      container = document.createElement('div');
      container.id = 'toast-container';
      container.style.position = 'fixed';
      container.style.top = '1rem';
      container.style.right = '1rem';
      container.style.zIndex = '9999';
      container.style.display = 'flex';
      container.style.flexDirection = 'column';
      container.style.gap = '0.5rem';
      document.body.appendChild(container);
    }
    
    // Create toast element
    const toast = document.createElement('div');
    toast.style.backgroundColor = getVariantColor(variant);
    toast.style.color = variant === 'default' ? '#000' : '#fff';
    toast.style.padding = '1rem';
    toast.style.borderRadius = '0.375rem';
    toast.style.boxShadow = '0 2px 8px rgba(0, 0, 0, 0.15)';
    toast.style.minWidth = '300px';
    toast.style.maxWidth = '500px';
    toast.style.animation = 'toast-enter 0.3s ease-out';
    
    // Add title
    const title = document.createElement('div');
    title.style.fontWeight = 'bold';
    title.style.marginBottom = options.description ? '0.25rem' : '0';
    title.textContent = options.title;
    toast.appendChild(title);
    
    // Add description if provided
    if (options.description) {
      const description = document.createElement('div');
      description.style.fontSize = '0.875rem';
      description.textContent = options.description;
      toast.appendChild(description);
    }
    
    // Add action button if provided
    if (options.action) {
      const button = document.createElement('button');
      button.textContent = options.action.label;
      button.style.marginTop = '0.5rem';
      button.style.padding = '0.25rem 0.5rem';
      button.style.backgroundColor = 'rgba(255, 255, 255, 0.2)';
      button.style.border = 'none';
      button.style.borderRadius = '0.25rem';
      button.style.cursor = 'pointer';
      button.onclick = options.action.onClick;
      toast.appendChild(button);
    }
    
    // Add to container
    container.appendChild(toast);
    
    // Remove after duration
    setTimeout(() => {
      toast.style.animation = 'toast-exit 0.3s ease-in';
      setTimeout(() => {
        container?.removeChild(toast);
        
        // Remove container if empty
        if (container?.childElementCount === 0) {
          document.body.removeChild(container);
        }
      }, 300);
    }, duration);
    
    // Add animation keyframes
    if (!document.getElementById('toast-animations')) {
      const style = document.createElement('style');
      style.id = 'toast-animations';
      style.textContent = `
        @keyframes toast-enter {
          from { transform: translateX(100%); opacity: 0; }
          to { transform: translateX(0); opacity: 1; }
        }
        @keyframes toast-exit {
          from { transform: translateX(0); opacity: 1; }
          to { transform: translateX(100%); opacity: 0; }
        }
      `;
      document.head.appendChild(style);
    }
    
    // Return toast ID or some identifier that could be used to close it
    return Date.now();
  }, []);
  
  return { toast };
}

/**
 * Get color for variant
 */
function getVariantColor(variant: ToastVariant): string {
  switch (variant) {
    case 'destructive':
      return '#ef4444'; // Red
    case 'success':
      return '#10b981'; // Green
    case 'warning':
      return '#f59e0b'; // Amber
    case 'info':
      return '#3b82f6'; // Blue
    default:
      return '#ffffff'; // White
  }
}