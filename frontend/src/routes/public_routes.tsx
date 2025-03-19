import { RouteObject } from "react-router-dom";
import App from "@/App";
import { MainLayout } from "@/layouts";
import { Dashboard, GVATool, FileManager, Macros, Manual } from "@/pages";

const publicRoutes: RouteObject[] = [
  {
    Component: App,
    children: [
      {
        path: "/",
        Component: MainLayout,
        children: [
          {
            path: "dashboard",
            Component: Dashboard
          },
          {
            path: 'features/gvatool',
            Component: GVATool
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