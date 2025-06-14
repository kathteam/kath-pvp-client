import { Navigate, RouteObject } from 'react-router-dom';
import App from '@/App';
import { MainLayout } from '@/layouts';
import { Dashboard, GVATool, FileManager, Macros, Manual, AnalysisHistory } from '@/pages';

const publicRoutes: RouteObject[] = [
  {
    Component: App,
    children: [
      {
        path: '/',
        Component: MainLayout,
        children: [
          {
            index: true,
            element: <Navigate to="dashboard" replace />
          },
          {
            path: 'dashboard',
            Component: Dashboard
          },
          {
            path: 'features/gvatool',
            Component: GVATool
          },
          {
            path: 'features/analysis_history',
            Component: AnalysisHistory
          },
          {
            path: 'system/file_manager',
            Component: FileManager
          },
          {
            path: 'system/macros',
            Component: Macros
          },
          {
            path: 'resources/manual',
            Component: Manual
          },
        ],
      },
    ],
  },
];

export default publicRoutes;
