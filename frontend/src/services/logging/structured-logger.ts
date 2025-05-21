/**
 * Structured Logger
 * 
 * A comprehensive logging system that captures structured log events
 * with contextual information and supports multiple output targets.
 */

// Log levels with numeric values for comparison
export enum LogLevel {
  TRACE = 0,
  DEBUG = 1,
  INFO = 2,
  WARN = 3,
  ERROR = 4,
  FATAL = 5,
  SILENT = 6
}

// Log entry interface
export interface LogEntry {
  timestamp: string;
  level: LogLevel;
  levelName: string;
  message: string;
  context?: Record<string, any>;
  tags?: string[];
  sessionId?: string;
  userId?: string;
  requestId?: string;
  component?: string;
  duration?: number;
  error?: {
    message: string;
    name: string;
    stack?: string;
    cause?: any;
  };
  [key: string]: any;
}

// Logger configuration
export interface LoggerConfig {
  minLevel: LogLevel;
  enabled: boolean;
  globalContext?: Record<string, any>;
  outputs: LogOutput[];
  captureErrors?: boolean;
  captureRejections?: boolean;
  maskSensitiveData?: boolean;
  maskPatterns?: RegExp[];
  generateRequestId?: boolean;
  includeTimestamp?: boolean;
}

// Log output interface
export interface LogOutput {
  write(entry: LogEntry): void | Promise<void>;
}

// Console log output
export class ConsoleOutput implements LogOutput {
  private levelToConsoleMethod: Record<LogLevel, keyof Console> = {
    [LogLevel.TRACE]: 'debug',
    [LogLevel.DEBUG]: 'debug',
    [LogLevel.INFO]: 'info',
    [LogLevel.WARN]: 'warn',
    [LogLevel.ERROR]: 'error',
    [LogLevel.FATAL]: 'error',
    [LogLevel.SILENT]: 'log'
  };

  write(entry: LogEntry): void {
    const method = this.levelToConsoleMethod[entry.level] || 'log';
    
    // Format log for console
    const timestamp = entry.timestamp ? `[${entry.timestamp}]` : '';
    const level = `[${entry.levelName}]`;
    const component = entry.component ? `[${entry.component}]` : '';
    const message = entry.message;
    
    // Basic log
    console[method](`${timestamp} ${level} ${component} ${message}`);
    
    // If we have context or error details, log them as well
    if (entry.context && Object.keys(entry.context).length > 0) {
      console.groupCollapsed('Context');
      console.dir(entry.context);
      console.groupEnd();
    }
    
    if (entry.error) {
      console.groupCollapsed('Error Details');
      console.dir(entry.error);
      console.groupEnd();
    }
  }
}

// Local storage log output (for persistent logs in browser)
export class LocalStorageOutput implements LogOutput {
  private readonly storageKey: string;
  private readonly maxEntries: number;
  
  constructor(storageKey = 'application_logs', maxEntries = 1000) {
    this.storageKey = storageKey;
    this.maxEntries = maxEntries;
  }
  
  write(entry: LogEntry): void {
    if (typeof window === 'undefined' || !window.localStorage) {
      return;
    }
    
    try {
      // Get existing logs
      const existingLogsJson = localStorage.getItem(this.storageKey);
      const logs: LogEntry[] = existingLogsJson ? JSON.parse(existingLogsJson) : [];
      
      // Add new log
      logs.push(entry);
      
      // Trim if we exceed max entries
      if (logs.length > this.maxEntries) {
        logs.splice(0, logs.length - this.maxEntries);
      }
      
      // Save back to local storage
      localStorage.setItem(this.storageKey, JSON.stringify(logs));
    } catch (error) {
      console.error('Failed to write log to local storage:', error);
    }
  }
  
  // Get logs from local storage
  getLogs(): LogEntry[] {
    if (typeof window === 'undefined' || !window.localStorage) {
      return [];
    }
    
    try {
      const logsJson = localStorage.getItem(this.storageKey);
      return logsJson ? JSON.parse(logsJson) : [];
    } catch (error) {
      console.error('Failed to read logs from local storage:', error);
      return [];
    }
  }
  
  // Clear logs from local storage
  clearLogs(): void {
    if (typeof window === 'undefined' || !window.localStorage) {
      return;
    }
    
    localStorage.removeItem(this.storageKey);
  }
}

// Remote log output (sends logs to server endpoint)
export class RemoteOutput implements LogOutput {
  private readonly endpoint: string;
  private readonly batchSize: number;
  private readonly batchTimeMs: number;
  private logQueue: LogEntry[] = [];
  private flushTimeoutId: number | null = null;
  
  constructor(endpoint: string, batchSize = 10, batchTimeMs = 5000) {
    this.endpoint = endpoint;
    this.batchSize = batchSize;
    this.batchTimeMs = batchTimeMs;
  }
  
  async write(entry: LogEntry): Promise<void> {
    // Add to queue
    this.logQueue.push(entry);
    
    // If we've reached batch size, flush immediately
    if (this.logQueue.length >= this.batchSize) {
      this.flush();
      return;
    }
    
    // Otherwise set a timeout to flush soon
    if (this.flushTimeoutId === null) {
      this.flushTimeoutId = window.setTimeout(() => {
        this.flush();
        this.flushTimeoutId = null;
      }, this.batchTimeMs);
    }
  }
  
  private async flush(): Promise<void> {
    if (this.logQueue.length === 0) {
      return;
    }
    
    // Take current queue and reset
    const logsToSend = [...this.logQueue];
    this.logQueue = [];
    
    // Clear timeout if it exists
    if (this.flushTimeoutId !== null) {
      clearTimeout(this.flushTimeoutId);
      this.flushTimeoutId = null;
    }
    
    try {
      // Send logs to remote endpoint
      await fetch(this.endpoint, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json'
        },
        body: JSON.stringify({ logs: logsToSend }),
        // Use keepalive to ensure logs are sent even if page is unloading
        keepalive: true
      });
    } catch (error) {
      // If remote logging fails, don't lose the logs - push back to queue
      this.logQueue = [...logsToSend, ...this.logQueue];
      
      // Log error to console (but don't create an infinite loop)
      console.error('Failed to send logs to remote endpoint:', error);
    }
  }
}

// Main logger class
export class StructuredLogger {
  private config: LoggerConfig;
  private globalContext: Record<string, any>;
  
  constructor(config: Partial<LoggerConfig> = {}) {
    // Default configuration
    this.config = {
      minLevel: LogLevel.INFO,
      enabled: true,
      outputs: [new ConsoleOutput()],
      globalContext: {},
      captureErrors: true,
      captureRejections: true,
      maskSensitiveData: true,
      maskPatterns: [
        /password/i,
        /token/i,
        /credential/i,
        /secret/i,
        /key/i,
        /authorization/i,
        /auth/i,
        /credit_?card/i,
        /card_?number/i,
        /cvv/i,
        /ssn/i
      ],
      generateRequestId: true,
      includeTimestamp: true,
      ...config
    };
    
    this.globalContext = this.config.globalContext || {};
    
    // Set up global error handlers if enabled
    if (this.config.captureErrors && typeof window !== 'undefined') {
      this.setupGlobalErrorHandlers();
    }
  }
  
  // Create log entry
  private createLogEntry(
    level: LogLevel, 
    message: string, 
    context?: Record<string, any>,
    error?: Error
  ): LogEntry {
    const timestamp = this.config.includeTimestamp ? new Date().toISOString() : '';
    
    // Create basic log entry
    const entry: LogEntry = {
      timestamp,
      level,
      levelName: LogLevel[level],
      message,
      context: context ? { ...context } : undefined
    };
    
    // Add global context
    if (Object.keys(this.globalContext).length > 0) {
      entry.context = {
        ...this.globalContext,
        ...(entry.context || {})
      };
    }
    
    // Add request ID if enabled
    if (this.config.generateRequestId && !entry.requestId) {
      entry.requestId = this.generateId();
    }
    
    // Add error details if provided
    if (error) {
      entry.error = {
        message: error.message,
        name: error.name,
        stack: error.stack,
        cause: (error as any).cause
      };
    }
    
    // Mask sensitive data if enabled
    if (this.config.maskSensitiveData && this.config.maskPatterns) {
      this.maskSensitiveData(entry);
    }
    
    return entry;
  }
  
  // Log at specific level
  private log(
    level: LogLevel, 
    message: string, 
    context?: Record<string, any>,
    error?: Error
  ): void {
    // Check if logging is enabled and level is sufficient
    if (!this.config.enabled || level < this.config.minLevel) {
      return;
    }
    
    // Create log entry
    const entry = this.createLogEntry(level, message, context, error);
    
    // Send to all outputs
    for (const output of this.config.outputs) {
      try {
        output.write(entry);
      } catch (error) {
        console.error('Error writing to log output:', error);
      }
    }
  }
  
  // Public logging methods
  trace(message: string, context?: Record<string, any>): void {
    this.log(LogLevel.TRACE, message, context);
  }
  
  debug(message: string, context?: Record<string, any>): void {
    this.log(LogLevel.DEBUG, message, context);
  }
  
  info(message: string, context?: Record<string, any>): void {
    this.log(LogLevel.INFO, message, context);
  }
  
  warn(message: string, context?: Record<string, any>): void {
    this.log(LogLevel.WARN, message, context);
  }
  
  error(message: string, errorOrContext?: Error | Record<string, any>, context?: Record<string, any>): void {
    // Handle different argument patterns
    if (errorOrContext instanceof Error) {
      this.log(LogLevel.ERROR, message, context, errorOrContext);
    } else {
      this.log(LogLevel.ERROR, message, errorOrContext);
    }
  }
  
  fatal(message: string, errorOrContext?: Error | Record<string, any>, context?: Record<string, any>): void {
    // Handle different argument patterns
    if (errorOrContext instanceof Error) {
      this.log(LogLevel.FATAL, message, context, errorOrContext);
    } else {
      this.log(LogLevel.FATAL, message, errorOrContext);
    }
  }
  
  // Create child logger with additional context
  child(childContext: Record<string, any>, component?: string): StructuredLogger {
    const childLogger = new StructuredLogger(this.config);
    
    // Merge global contexts
    childLogger.globalContext = {
      ...this.globalContext,
      ...childContext
    };
    
    // Set component if provided
    if (component) {
      childLogger.globalContext.component = component;
    }
    
    return childLogger;
  }
  
  // Create component logger (convenience method)
  component(componentName: string, context?: Record<string, any>): StructuredLogger {
    return this.child(context || {}, componentName);
  }
  
  // Update logger configuration
  updateConfig(configUpdate: Partial<LoggerConfig>): void {
    this.config = {
      ...this.config,
      ...configUpdate
    };
    
    // Update global context if provided
    if (configUpdate.globalContext) {
      this.globalContext = {
        ...this.globalContext,
        ...configUpdate.globalContext
      };
    }
  }
  
  // Add context to the global context
  addGlobalContext(context: Record<string, any>): void {
    this.globalContext = {
      ...this.globalContext,
      ...context
    };
  }
  
  // Generate unique ID (for request ID)
  private generateId(): string {
    return Math.random().toString(36).substring(2, 15) + 
           Math.random().toString(36).substring(2, 15);
  }
  
  // Mask sensitive data in log entry
  private maskSensitiveData(entry: LogEntry): void {
    if (!this.config.maskPatterns || !this.config.maskSensitiveData) {
      return;
    }
    
    const maskValue = '***MASKED***';
    
    // Function to recursively mask sensitive data
    const maskObject = (obj: any): any => {
      if (!obj || typeof obj !== 'object') {
        return obj;
      }
      
      // Handle arrays
      if (Array.isArray(obj)) {
        return obj.map(item => maskObject(item));
      }
      
      // Handle objects
      const result: any = {};
      for (const [key, value] of Object.entries(obj)) {
        // Check if key matches any mask patterns
        const shouldMask = this.config.maskPatterns!.some(pattern => pattern.test(key));
        
        if (shouldMask) {
          result[key] = maskValue;
        } else if (typeof value === 'object' && value !== null) {
          result[key] = maskObject(value);
        } else {
          result[key] = value;
        }
      }
      
      return result;
    };
    
    // Mask context data
    if (entry.context) {
      entry.context = maskObject(entry.context);
    }
    
    // Mask error data
    if (entry.error && entry.error.cause && typeof entry.error.cause === 'object') {
      entry.error.cause = maskObject(entry.error.cause);
    }
  }
  
  // Setup global error handlers
  private setupGlobalErrorHandlers(): void {
    // Handle uncaught exceptions
    window.addEventListener('error', (event) => {
      this.error(
        `Uncaught error: ${event.message}`,
        event.error || new Error(event.message),
        {
          filename: event.filename,
          lineno: event.lineno,
          colno: event.colno
        }
      );
    });
    
    // Handle unhandled promise rejections
    window.addEventListener('unhandledrejection', (event) => {
      let error: Error;
      let message: string;
      
      if (event.reason instanceof Error) {
        error = event.reason;
        message = `Unhandled promise rejection: ${error.message}`;
      } else {
        error = new Error(String(event.reason));
        message = `Unhandled promise rejection: ${String(event.reason)}`;
      }
      
      this.error(message, error);
    });
  }
  
  // Create a performance timer for measuring durations
  startTimer(operation: string): () => void {
    const startTime = performance.now();
    
    return () => {
      const duration = performance.now() - startTime;
      this.info(`${operation} completed`, { 
        operation, 
        duration,
        durationMs: Math.round(duration)
      });
    };
  }
  
  // Log API request
  logApiRequest(
    method: string, 
    url: string, 
    requestData?: any,
    startTime: number = performance.now()
  ): (response?: any, error?: Error) => void {
    this.info(`API Request: ${method} ${url}`, {
      method,
      url,
      requestData: this.config.maskSensitiveData ? this.maskObject(requestData) : requestData
    });
    
    return (response?: any, error?: Error) => {
      const duration = performance.now() - startTime;
      
      if (error) {
        this.error(`API Error: ${method} ${url}`, error, {
          method,
          url,
          duration,
          durationMs: Math.round(duration),
          requestData: this.config.maskSensitiveData ? this.maskObject(requestData) : requestData
        });
      } else {
        this.info(`API Response: ${method} ${url}`, {
          method,
          url,
          duration,
          durationMs: Math.round(duration),
          status: response?.status,
          requestData: this.config.maskSensitiveData ? this.maskObject(requestData) : requestData,
          responseData: this.config.maskSensitiveData ? this.maskObject(response?.data) : response?.data
        });
      }
    };
  }
  
  // Helper to mask a single object
  private maskObject(obj: any): any {
    if (!obj || typeof obj !== 'object' || !this.config.maskPatterns) {
      return obj;
    }
    
    const maskValue = '***MASKED***';
    
    // Handle arrays
    if (Array.isArray(obj)) {
      return obj.map(item => this.maskObject(item));
    }
    
    // Handle objects
    const result: any = {};
    for (const [key, value] of Object.entries(obj)) {
      // Check if key matches any mask patterns
      const shouldMask = this.config.maskPatterns!.some(pattern => pattern.test(key));
      
      if (shouldMask) {
        result[key] = maskValue;
      } else if (typeof value === 'object' && value !== null) {
        result[key] = this.maskObject(value);
      } else {
        result[key] = value;
      }
    }
    
    return result;
  }
}

// Create and export default logger instance
const defaultLogger = new StructuredLogger({
  minLevel: process.env.NODE_ENV === 'production' ? LogLevel.INFO : LogLevel.DEBUG,
  outputs: [new ConsoleOutput()]
});

export default defaultLogger;