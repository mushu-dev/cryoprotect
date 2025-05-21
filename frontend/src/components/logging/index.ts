/**
 * Logging Components
 * 
 * This module exports components for the structured logging system
 */

export { default as LoggingProvider, useLoggingContext } from './logging-provider';
export { default as LogViewer } from './log-viewer';

// Re-export types from the logger
export {
  LogLevel,
  type LogEntry,
  type LoggerConfig,
  type LogOutput,
  ConsoleOutput,
  LocalStorageOutput,
  RemoteOutput
} from '@/services/logging/structured-logger';