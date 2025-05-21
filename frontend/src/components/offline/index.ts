/**
 * Offline Components
 * 
 * This module exports all components and hooks related to offline functionality.
 */

import OfflineWrapper from '../offline-wrapper';
import OfflineStatusIndicator from '../offline-status-indicator';
import useOfflineMode from '@/hooks/use-offline-mode';

export {
  OfflineWrapper,
  OfflineStatusIndicator,
  useOfflineMode
};

// Export the offline components as a group
const OfflineComponents = {
  Wrapper: OfflineWrapper,
  StatusIndicator: OfflineStatusIndicator,
};

export default OfflineComponents;