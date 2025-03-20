import { StrictMode } from 'react';
import { createRoot } from 'react-dom/client';
import { createHashRouter, RouterProvider } from 'react-router-dom';
import routes from '@/routes';
import '@/index.css';

declare global {
  interface Window {
    pywebview: {
      api: {
        // Known methods
        fullscreen: () => Promise<void>;

        // Generic type definition
        [key: string]: (...args: unknown[]) => Promise<unknown>;
      };
      [key: string]: unknown;
    };
  }
}

const router = createHashRouter(routes);

createRoot(document.getElementById('root')!).render(
  <StrictMode>
    <RouterProvider router={router} />
  </StrictMode>,
);
