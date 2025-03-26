import { StrictMode } from 'react';
import { createRoot } from 'react-dom/client';
import { createHashRouter, RouterProvider } from 'react-router-dom';
import routes from '@/routes';
import '@/index.css';

declare global {
  interface Window {
    pywebview: {
      api: {
        // Known services
        ui_controller: {
          fullscreen: () => Promise<void>;
        },

        file_manager: {
          list_files: (...args: unknown[]) => Promise<{ //path string as arg
            filename: string; 
            type: string; 
            size_kb: number | null; 
            item_count: number | null; 
          }[]>;
          upload_file: (...args: unknown[]) => Promise<unknown>; //path: str, file_name: str, file_content: bytes          [method: string]: (...args: unknown[]) => Promise<unknown>;
        },
        
        fasta_service: {
          create_disease_download: () => Promise<JSON>
          download_reference_genome_grch38: () => Promise<JSON>;
          [method: string]: (...args: unknown[]) => Promise<unknown>;
        },
        // eslint-disable-next-line @typescript-eslint/consistent-indexed-object-style
        [service: string]: {
          [method: string]: (...args: unknown[]) => Promise<unknown>;
        },
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
