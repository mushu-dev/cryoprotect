"use client"

import * as React from "react"
import * as ToastPrimitives from "@radix-ui/react-toast"
import { cva } from "class-variance-authority"
import { X, AlertCircle, CheckCircle, Info, AlertTriangle } from "lucide-react"

import { cn } from "../../lib/utils"
import { useId } from "react"

// Create context for toast accessibility
const ToastContext = React.createContext({ id: "", type: "default" })

const ToastProvider = ToastPrimitives.Provider

const ToastViewport = React.forwardRef(({ className, ...props }, ref) => (
  <ToastPrimitives.Viewport
    ref={ref}
    className={cn(
      "fixed top-0 z-[100] flex max-h-screen w-full flex-col-reverse p-4 sm:bottom-0 sm:right-0 sm:top-auto sm:flex-col md:max-w-[420px]",
      className
    )}
    // Improve a11y with proper ARIA labels
    aria-label="Notifications"
    role="region"
    aria-live="polite"
    {...props}
  />
))
ToastViewport.displayName = ToastPrimitives.Viewport.displayName

const toastVariants = cva(
  "group pointer-events-auto relative flex w-full items-center justify-between space-x-4 overflow-hidden rounded-md border p-6 pr-8 shadow-lg transition-all data-[swipe=cancel]:translate-x-0 data-[swipe=end]:translate-x-[var(--radix-toast-swipe-end-x)] data-[swipe=move]:translate-x-[var(--radix-toast-swipe-move-x)] data-[swipe=move]:transition-none data-[state=open]:animate-in data-[state=closed]:animate-out data-[swipe=end]:animate-out data-[state=closed]:fade-out-80 data-[state=closed]:slide-out-to-right-full data-[state=open]:slide-in-from-top-full data-[state=open]:sm:slide-in-from-bottom-full",
  {
    variants: {
      variant: {
        default: "border bg-background text-foreground",
        destructive: "destructive group border-destructive bg-destructive text-destructive-foreground",
        success: "success group border-green-500 bg-green-50 dark:bg-green-900/20 text-green-900 dark:text-green-300",
        warning: "warning group border-yellow-500 bg-yellow-50 dark:bg-yellow-900/20 text-yellow-900 dark:text-yellow-300",
        info: "info group border-blue-500 bg-blue-50 dark:bg-blue-900/20 text-blue-900 dark:text-blue-300",
      },
    },
    defaultVariants: {
      variant: "default",
    },
  }
)

// Enhanced Toast with accessibility features
const Toast = React.forwardRef(({ className, variant = "default", type, children, ...props }, ref) => {
  // Generate unique ID for each toast
  const id = useId();
  
  // Determine toast icon based on variant
  const getIcon = () => {
    switch (variant) {
      case "destructive":
        return <AlertCircle className="h-5 w-5 text-destructive-foreground" />;
      case "success":
        return <CheckCircle className="h-5 w-5 text-green-600 dark:text-green-400" />;
      case "warning":
        return <AlertTriangle className="h-5 w-5 text-yellow-600 dark:text-yellow-400" />;
      case "info":
        return <Info className="h-5 w-5 text-blue-600 dark:text-blue-400" />;
      default:
        return null;
    }
  };
  
  // Set appropriate ARIA roles based on variant
  const getRole = () => {
    switch (variant) {
      case "destructive":
        return "alert";
      case "warning":
        return "alert";
      default:
        return "status";
    }
  };
  
  return (
    <ToastContext.Provider value={{ id, type: variant }}>
      <ToastPrimitives.Root
        ref={ref}
        role={getRole()}
        aria-atomic="true"
        aria-live={variant === "destructive" ? "assertive" : "polite"}
        className={cn(toastVariants({ variant }), className)}
        {...props}
      >
        <div className="flex items-start gap-3">
          {getIcon()}
          <div className="flex-1">{children}</div>
        </div>
      </ToastPrimitives.Root>
    </ToastContext.Provider>
  )
})
Toast.displayName = ToastPrimitives.Root.displayName

const ToastAction = React.forwardRef(({ className, altText, ...props }, ref) => {
  const { id } = React.useContext(ToastContext);
  
  return (
    <ToastPrimitives.Action
      ref={ref}
      // Adds descriptive text for screen readers
      aria-label={altText}
      className={cn(
        "inline-flex h-8 shrink-0 items-center justify-center rounded-md border bg-transparent px-3 text-sm font-medium ring-offset-background transition-colors hover:bg-secondary focus:outline-none focus:ring-2 focus:ring-ring focus:ring-offset-2 disabled:pointer-events-none disabled:opacity-50 group-[.destructive]:border-muted/40 group-[.destructive]:hover:border-destructive/30 group-[.destructive]:hover:bg-destructive group-[.destructive]:hover:text-destructive-foreground group-[.destructive]:focus:ring-destructive",
        "group-[.success]:border-green-200 group-[.success]:hover:border-green-300 group-[.success]:hover:bg-green-100 dark:group-[.success]:hover:bg-green-900/30",
        "group-[.warning]:border-yellow-200 group-[.warning]:hover:border-yellow-300 group-[.warning]:hover:bg-yellow-100 dark:group-[.warning]:hover:bg-yellow-900/30",
        "group-[.info]:border-blue-200 group-[.info]:hover:border-blue-300 group-[.info]:hover:bg-blue-100 dark:group-[.info]:hover:bg-blue-900/30",
        className
      )}
      {...props}
    />
  )
})
ToastAction.displayName = ToastPrimitives.Action.displayName

const ToastClose = React.forwardRef(({ className, ...props }, ref) => {
  const { type } = React.useContext(ToastContext);
  
  return (
    <ToastPrimitives.Close
      ref={ref}
      className={cn(
        "absolute right-2 top-2 rounded-md p-1 text-foreground/50 opacity-0 transition-opacity hover:text-foreground focus:opacity-100 focus:outline-none focus:ring-2 group-hover:opacity-100",
        "group-[.destructive]:text-red-300 group-[.destructive]:hover:text-red-50 group-[.destructive]:focus:ring-red-400 group-[.destructive]:focus:ring-offset-red-600",
        "group-[.success]:text-green-600 group-[.success]:hover:text-green-700 group-[.success]:focus:ring-green-400",
        "group-[.warning]:text-yellow-600 group-[.warning]:hover:text-yellow-700 group-[.warning]:focus:ring-yellow-400",
        "group-[.info]:text-blue-600 group-[.info]:hover:text-blue-700 group-[.info]:focus:ring-blue-400",
        className
      )}
      // Add descriptive text for screen readers
      aria-label="Close notification"
      {...props}
    >
      <X className="h-4 w-4" />
    </ToastPrimitives.Close>
  )
})
ToastClose.displayName = ToastPrimitives.Close.displayName

const ToastTitle = React.forwardRef(({ className, ...props }, ref) => {
  const { id } = React.useContext(ToastContext);
  
  return (
    <ToastPrimitives.Title
      ref={ref}
      id={`${id}-title`}
      className={cn("text-sm font-semibold", className)}
      {...props}
    />
  )
})
ToastTitle.displayName = ToastPrimitives.Title.displayName

const ToastDescription = React.forwardRef(({ className, ...props }, ref) => {
  const { id } = React.useContext(ToastContext);
  
  return (
    <ToastPrimitives.Description
      ref={ref}
      id={`${id}-description`}
      className={cn("text-sm opacity-90", className)}
      {...props}
    />
  )
})
ToastDescription.displayName = ToastPrimitives.Description.displayName

export {
  ToastProvider,
  ToastViewport,
  Toast,
  ToastTitle,
  ToastDescription,
  ToastClose,
  ToastAction,
}