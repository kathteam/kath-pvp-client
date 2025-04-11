import { useContext } from 'react';
import { SetupContext } from './SetupProvider';

export default function useSetupState() {
  const context = useContext(SetupContext);
  if (context === undefined) {
    throw new Error('useSetupState must be used within a SetupProvider');
  }
  return context;
}
