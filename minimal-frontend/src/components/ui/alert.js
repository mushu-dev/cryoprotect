import React from 'react';

/**
 * Alert component variants
 */
const alertVariants = {
  default: {
    container: "bg-blue-50 border-blue-200 text-blue-800",
    icon: "text-blue-600"
  },
  destructive: {
    container: "bg-red-50 border-red-200 text-red-800",
    icon: "text-red-600"
  },
  success: {
    container: "bg-green-50 border-green-200 text-green-800",
    icon: "text-green-600"
  },
  warning: {
    container: "bg-yellow-50 border-yellow-200 text-yellow-800",
    icon: "text-yellow-600"
  }
};

/**
 * Alert component for displaying important messages
 */
const Alert = React.forwardRef(({ 
  className = '', 
  variant = 'default', 
  children, 
  ...props 
}, ref) => {
  const variantClasses = alertVariants[variant] || alertVariants.default;
  
  return (
    <div
      ref={ref}
      role="alert"
      className={`relative w-full rounded-lg border p-4 ${variantClasses.container} ${className}`}
      {...props}
    >
      {children}
    </div>
  );
});

Alert.displayName = "Alert";

/**
 * Alert title component
 */
const AlertTitle = React.forwardRef(({ 
  className = '', 
  ...props 
}, ref) => (
  <h5
    ref={ref}
    className={`mb-1 font-medium leading-none tracking-tight ${className}`}
    {...props}
  />
));

AlertTitle.displayName = "AlertTitle";

/**
 * Alert description component
 */
const AlertDescription = React.forwardRef(({ 
  className = '', 
  ...props 
}, ref) => (
  <div
    ref={ref}
    className={`text-sm leading-relaxed ${className}`}
    {...props}
  />
));

AlertDescription.displayName = "AlertDescription";

export { Alert, AlertTitle, AlertDescription };