import React, { useState, useEffect, useRef, createContext, useContext, useId } from "react";
import { X } from "lucide-react";
import { cn } from "../../lib/utils";
import { useCallback } from "react";

// Create Dialog context for managing state and accessibility connections
const DialogContext = createContext({
  open: false,
  onOpenChange: () => {},
  id: "",
  titleId: "",
  descriptionId: ""
});

/**
 * Enhanced Dialog component with improved accessibility
 */
const Dialog = ({ children, open, onOpenChange, defaultOpen = false }) => {
  const [isOpen, setIsOpen] = useState(open ?? defaultOpen ?? false);
  const id = useId();
  const titleId = `${id}-title`;
  const descriptionId = `${id}-description`;
  
  useEffect(() => {
    if (open !== undefined) {
      setIsOpen(open);
    }
  }, [open]);
  
  const handleOpenChange = useCallback((value) => {
    // If controlled externally, call onOpenChange
    if (onOpenChange) {
      onOpenChange(value);
    }
    
    // If not controlled, update internal state
    if (open === undefined) {
      setIsOpen(value);
    }
  }, [onOpenChange, open]);
  
  return (
    <DialogContext.Provider value={{ 
      open: isOpen, 
      onOpenChange: handleOpenChange,
      id,
      titleId,
      descriptionId
    }}>
      {children}
    </DialogContext.Provider>
  );
};

/**
 * Dialog trigger button that follows accessibility best practices
 */
const DialogTrigger = React.forwardRef(({ className, children, asChild, ...props }, ref) => {
  const { onOpenChange, open, id, titleId, descriptionId } = useContext(DialogContext);
  
  return (
    <button
      ref={ref}
      type="button"
      aria-haspopup="dialog"
      aria-expanded={open ? true : false}
      aria-controls={open ? id : undefined}
      aria-labelledby={titleId}
      aria-describedby={descriptionId}
      className={cn(
        "inline-flex items-center justify-center rounded-md text-sm font-medium transition-colors focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring focus-visible:ring-offset-2 disabled:opacity-50 disabled:pointer-events-none",
        className
      )}
      onClick={() => onOpenChange(true)}
      {...props}
    >
      {children}
    </button>
  );
});
DialogTrigger.displayName = "DialogTrigger";

/**
 * Dialog content component with improved accessibility and keyboard navigation
 */
const DialogContent = React.forwardRef(({ className, children, ...props }, ref) => {
  const { open, onOpenChange, id, titleId, descriptionId } = useContext(DialogContext);
  const contentRef = useRef(null);
  const initialFocusRef = useRef(null);
  
  // Trap focus inside dialog when open
  useEffect(() => {
    if (!open) return;
    
    const dialog = contentRef.current;
    if (!dialog) return;
    
    // Get all focusable elements
    const focusableElements = dialog.querySelectorAll(
      'button, [href], input, select, textarea, [tabindex]:not([tabindex="-1"])'
    );
    
    if (focusableElements.length === 0) return;
    
    const firstElement = focusableElements[0];
    const lastElement = focusableElements[focusableElements.length - 1];
    
    // Set initial focus
    const elementToFocus = initialFocusRef.current || firstElement;
    elementToFocus.focus();
    
    // Handle tab key to trap focus
    const handleTabKey = (e) => {
      if (e.key !== 'Tab') return;
      
      // If shift + tab on first element, move to last element
      if (e.shiftKey && document.activeElement === firstElement) {
        e.preventDefault();
        lastElement.focus();
        return;
      }
      
      // If tab on last element, move to first element
      if (!e.shiftKey && document.activeElement === lastElement) {
        e.preventDefault();
        firstElement.focus();
        return;
      }
    };
    
    document.addEventListener('keydown', handleTabKey);
    
    // Store previous active element to restore focus
    const previousActiveElement = document.activeElement;
    
    return () => {
      document.removeEventListener('keydown', handleTabKey);
      
      // Restore focus when dialog closes
      if (previousActiveElement && previousActiveElement.focus) {
        previousActiveElement.focus();
      }
    };
  }, [open]);
  
  // Close dialog on escape key or outside click
  useEffect(() => {
    if (!open) return;
    
    const handleEscape = (e) => {
      if (e.key === "Escape") {
        onOpenChange(false);
      }
    };
    
    const handleClickOutside = (e) => {
      if (contentRef.current && !contentRef.current.contains(e.target)) {
        onOpenChange(false);
      }
    };
    
    document.addEventListener("keydown", handleEscape);
    document.addEventListener("mousedown", handleClickOutside);
    
    // Add body scroll lock
    const originalStyle = window.getComputedStyle(document.body).overflow;
    document.body.style.overflow = "hidden";
    
    return () => {
      document.removeEventListener("keydown", handleEscape);
      document.removeEventListener("mousedown", handleClickOutside);
      
      // Restore body scroll
      document.body.style.overflow = originalStyle;
    };
  }, [open, onOpenChange]);
  
  if (!open) return null;
  
  return (
    <>
      {/* Overlay with appropriate role */}
      <div 
        className="fixed inset-0 z-50 bg-black/80 data-[state=open]:animate-in data-[state=closed]:animate-out data-[state=closed]:fade-out-0 data-[state=open]:fade-in-0" 
        aria-hidden="true"
      />
      
      {/* Content */}
      <div 
        className="fixed inset-0 z-50 flex items-center justify-center"
        onClick={(e) => e.stopPropagation()}
      >
        <div
          ref={mergeRefs(contentRef, ref)}
          role="dialog"
          aria-modal="true"
          id={id}
          aria-labelledby={titleId}
          aria-describedby={descriptionId}
          tabIndex={-1}
          className={cn(
            "w-full max-w-lg max-h-[85vh] overflow-y-auto gap-4 border bg-background p-6 shadow-lg sm:rounded-lg data-[state=open]:animate-in data-[state=closed]:animate-out data-[state=closed]:fade-out-0 data-[state=open]:fade-in-0 data-[state=closed]:zoom-out-95 data-[state=open]:zoom-in-95 data-[state=closed]:slide-out-to-center data-[state=open]:slide-in-from-center",
            className
          )}
          {...props}
        >
          {children}
          <button
            className="absolute right-4 top-4 rounded-sm opacity-70 ring-offset-background transition-opacity hover:opacity-100 focus:outline-none focus:ring-2 focus:ring-ring focus:ring-offset-2 disabled:pointer-events-none"
            onClick={() => onOpenChange(false)}
            aria-label="Close dialog"
            type="button"
          >
            <X className="h-4 w-4" />
            <span className="sr-only">Close</span>
          </button>
        </div>
      </div>
    </>
  );
});
DialogContent.displayName = "DialogContent";

/**
 * Dialog header component for consistent styling
 */
const DialogHeader = ({ className, ...props }) => (
  <div
    className={cn("flex flex-col space-y-1.5 text-center sm:text-left", className)}
    {...props}
  />
);
DialogHeader.displayName = "DialogHeader";

/**
 * Dialog footer component for consistent styling
 */
const DialogFooter = ({ className, ...props }) => (
  <div
    className={cn("flex flex-col-reverse sm:flex-row sm:justify-end sm:space-x-2", className)}
    {...props}
  />
);
DialogFooter.displayName = "DialogFooter";

/**
 * Dialog title component with proper accessibility attributes
 */
const DialogTitle = React.forwardRef(({ className, ...props }, ref) => {
  const { titleId } = useContext(DialogContext);
  
  return (
    <h2
      ref={ref}
      id={titleId}
      className={cn("text-lg font-semibold leading-none tracking-tight", className)}
      {...props}
    />
  );
});
DialogTitle.displayName = "DialogTitle";

/**
 * Dialog description component with proper accessibility attributes
 */
const DialogDescription = React.forwardRef(({ className, ...props }, ref) => {
  const { descriptionId } = useContext(DialogContext);
  
  return (
    <p
      ref={ref}
      id={descriptionId}
      className={cn("text-sm text-muted-foreground", className)}
      {...props}
    />
  );
});
DialogDescription.displayName = "DialogDescription";

/**
 * Utility function to merge multiple refs
 */
function mergeRefs(...refs) {
  return (value) => {
    refs.forEach((ref) => {
      if (typeof ref === "function") {
        ref(value);
      } else if (ref != null) {
        ref.current = value;
      }
    });
  };
}

export {
  Dialog,
  DialogTrigger,
  DialogContent,
  DialogHeader,
  DialogFooter,
  DialogTitle,
  DialogDescription,
};