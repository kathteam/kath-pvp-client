import { StrictMode } from 'react';
import { createRoot } from 'react-dom/client';
import { createHashRouter, RouterProvider } from 'react-router-dom';
import routes from '@/routes';
import '@/index.css';

declare global {
  interface Window {
    pywebview: {
      api: {

        // --------------
        // Known services
        // --------------

        ui_controller: {
          fullscreen: () => Promise<void>;
        },

        file_controller: {
          list_files: (path: string) => Promise<{
            filename: string; 
            type: string; 
            size_kb: number | null; 
            item_count: number | null; 
          }[]>;
          upload_file: (path: string, file_name: string, file_content: number[]) => Promise<void>;
          rename_file: (path: string, old_name: string, new_name: string) => Promise<void>;
          delete_file: (path: string, file_name: string) => Promise<void>;
          [method: string]: (...args: any[]) => Promise<any>;
          get_kath_directory: () => Promise<string>;
          create_vcf_database: (db_path: string) => Promise<void>;
        },
        
        fasta_service: {
          create_disease_download: (disease: string, ref_max: number) => Promise<{
            status: string;
            disease_term: string;
            max_results: number;
            downloaded_files: string[];
            count: number;
          }>;
          download_reference_genome_grch38: () => Promise<{
            status: string;
            reference_genome: string;
          }>;
          [method: string]: (...args: any[]) => Promise<any>;
        },

        blast_service: {
          align_mutations: () => Promise<{
            status: string;
            result_file: string;
          }>;
          perform_blast_analysis: () => Promise<{
            status: string;
            result_file: string;
          }>;
          disease_extraction: (fasta_file: string) => Promise<{
            status: string;
            result_file: string;
          }>;
          [method: string]: (...args: any[]) => Promise<any>;
        },

        disease_service: {
          get_disease_data: (file_path: string) => Promise<{
            clinicalSignificance: string;
            disease: string;
          }[]>;
          [method: string]: (...args: any[]) => Promise<any>;
        },

        // --------------
        // Unknown services
        // --------------

        // eslint-disable-next-line @typescript-eslint/consistent-indexed-object-style
        [service: string]: {
          [method: string]: (...args: any[]) => Promise<any>;
        },
      };
      [key: string]: any;
    };
  }
}

const router = createHashRouter(routes);

createRoot(document.getElementById('root')!).render(
  <StrictMode>
    <RouterProvider router={router} />
  </StrictMode>,
);
