import React, { memo, useId, useTransition, useCallback } from "react";
import { cn } from "../../lib/utils";

// Performance optimized button component
const Button = memo(React.forwardRef(({ 
  className, 
  variant = "default", 
  size = "default", 
  asChild = false,
  isLoading = false,
  loadingText = "Loading...",
  loadingIcon = true,
  onClick,
  children,
  ...props 
}, ref) => {
  const [isPending, startTransition] = useTransition();
  const id = useId();

  // Common base styles
  const baseStyles = "inline-flex items-center justify-center whitespace-nowrap rounded-md font-medium ring-offset-background transition-all focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring focus-visible:ring-offset-2 disabled:pointer-events-none disabled:opacity-50 active:scale-95";
  
  // Variant styles with hover and active states
  const variants = {
    default: "bg-primary text-primary-foreground hover:bg-primary/90 active:bg-primary/95",
    destructive: "bg-destructive text-destructive-foreground hover:bg-destructive/90 active:bg-destructive/95",
    outline: "border border-input bg-background hover:bg-accent hover:text-accent-foreground active:bg-accent/80",
    secondary: "bg-secondary text-secondary-foreground hover:bg-secondary/90 active:bg-secondary/95",
    ghost: "hover:bg-accent hover:text-accent-foreground active:bg-accent/80",
    link: "text-primary underline-offset-4 hover:underline active:text-primary/80",
  };
  
  // Size styles with better responsiveness
  const sizes = {
    default: "h-10 px-4 py-2 text-sm",
    sm: "h-8 rounded-md px-3 text-xs",
    md: "h-9 rounded-md px-3 text-sm",
    lg: "h-11 rounded-md px-5 text-base",
    xl: "h-12 rounded-md px-6 text-lg",
    icon: "h-10 w-10 p-0",
    "icon-sm": "h-8 w-8 p-0",
  };
  
  // Memoized click handler with transition for smoother UI updates
  const handleClick = useCallback((e) => {
    if (onClick) {
      if (props.disabled || isLoading || isPending) return;
      
      startTransition(() => {
        onClick(e);
      });
    }
  }, [onClick, props.disabled, isLoading, isPending]);
  
  // Show loading state
  const showLoading = isLoading || isPending;
  
  return (
    <button
      id={id}
      className={cn(
        baseStyles, 
        variants[variant], 
        sizes[size], 
        showLoading && "cursor-wait opacity-80",
        className
      )}
      ref={ref}
      onClick={handleClick}
      disabled={props.disabled || showLoading}
      {...props}
    >
      {showLoading ? (
        <>
          {loadingIcon && <span className="mr-2 h-4 w-4 rounded-full border-2 border-current border-t-transparent animate-spin"></span>}
          {loadingText || children}
        </>
      ) : children}
    </button>
  );
}));

Button.displayName = "Button";

export { Button };