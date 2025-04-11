import { createContext, useEffect, useState, ReactNode, JSX } from 'react';

type SetupStatus = 'idle' | 'started' | 'progress' | 'completed' | 'failed';

export interface SetupContextType {
  status: SetupStatus;
  progress: number;
  message: string;
}

// eslint-disable-next-line react-refresh/only-export-components
export const SetupContext = createContext<SetupContextType | undefined>(undefined);

interface SetupProviderProps {
  children: ReactNode;
  initialState?: SetupContextType;
}

export function SetupProvider({
  children,
  initialState = { status: 'idle', progress: 0, message: '' }
}: SetupProviderProps): JSX.Element {
  const [status, setStatus] = useState<SetupStatus>(initialState.status);
  const [progress, setProgress] = useState<number>(initialState.progress);
  const [message, setMessage] = useState<string>(initialState.message);

  useEffect(() => {
    const handleSetupStatus = (event: CustomEvent) => {
      const { status, progress, message } = event.detail;

      switch (status) {
        case 'started':
          setStatus(status);
          setProgress(0);
          setMessage('Setup started...');
          break;

        case 'progress':
          setStatus(status);
          setProgress(progress ?? 0);
          setMessage(message ?? 'Processing...');
          break;

        case 'completed':
          setStatus(status);
          setProgress(1);
          setMessage('Setup completed successfully');
          break;

        case 'failed':
          setStatus(status);
          setProgress(0);
          setMessage('Setup failed');
          break;

        default:
          setStatus('idle');
          setProgress(0);
          setMessage('');
      }
    };

    window.addEventListener('setup-status', handleSetupStatus as EventListener);

    return () => {
      window.removeEventListener('setup-status', handleSetupStatus as EventListener);
    };
  }, []);

  return (
    <SetupContext.Provider value={{ status, progress, message }}>
      {children}
    </SetupContext.Provider>
  );
};
