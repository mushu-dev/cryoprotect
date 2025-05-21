/**
 * Logger Hook
 * 
 * A React hook for accessing the structured logger in components.
 */

import { useMemo } from 'react';
import { useRouter } from 'next/router';
import defaultLogger, { StructuredLogger } from '@/services/logging/structured-logger';

interface UseLoggerOptions {
  /** Component name to include in logs */
  component?: string;
  
  /** Additional context to include in all logs */
  context?: Record<string, any>;
  
  /** Whether to include route information */
  includeRouteInfo?: boolean;
  
  /** Optional parent logger to derive from */
  parentLogger?: StructuredLogger;
}

/**
 * Hook to use the structured logger in React components
 */
export function useLogger(options: UseLoggerOptions = {}): StructuredLogger {
  const {
    component,
    context = {},
    includeRouteInfo = true,
    parentLogger = defaultLogger
  } = options;
  
  const router = useRouter();
  
  // Create the logger with component and route context
  const logger = useMemo(() => {
    const routeContext = includeRouteInfo ? {
      route: router.pathname,
      query: router.query,
      asPath: router.asPath,
    } : {};
    
    const loggerContext = {
      ...context,
      ...routeContext,
    };
    
    return component 
      ? parentLogger.component(component, loggerContext)
      : parentLogger.child(loggerContext);
  }, [component, context, includeRouteInfo, parentLogger, router.pathname, router.asPath]);
  
  return logger;
}

export default useLogger;