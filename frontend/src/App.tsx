import { JSX } from 'react';
import { Outlet } from 'react-router-dom';
import { ReactRouterAppProvider } from '@toolpad/core/react-router';
import { Navigation } from '@toolpad/core/AppProvider';
import { 
  Dashboard as DashboardIcon,
  Plagiarism as PlagiarismIcon,
  Inventory as InventoryIcon,
  AutoAwesome as AutoAwesomeIcon,
  MenuBook as MenuBookIcon
} from '@mui/icons-material';
import { Box } from '@mui/material';

const NAVIGATION: Navigation = [
  {
    segment: 'index.html',
    title: 'Dashboard',
    icon: <DashboardIcon />,
  },
  {
    kind: 'divider',
  },
  {
    kind: 'header',
    title: 'Features',
  },
  {
    segment: 'features/gvatool',
    title: 'GVATool',
    icon: <PlagiarismIcon />,
  },
  {
    kind: 'divider',
  },
  {
    kind: 'header',
    title: 'System',
  },
  {
    segment: 'system/file_manager',
    title: 'File Manager',
    icon: <InventoryIcon />,
  },
  {
    segment: 'system/macros',
    title: 'Macros',
    icon: <AutoAwesomeIcon />,
  },
  {
    kind: 'divider',
  },
  {
    kind: 'header',
    title: 'Resources',
  },
  {
    segment: 'resources/manual',
    title: 'Manual',
    icon: <MenuBookIcon />,
  },
];

// TODO: Fix the logo
const BRANDING = {
  logo: <Box sx={{p: '8px'}}><img src='logo.svg' width={'24px'} height={'24px'}></img></Box>,
  title: "Kath",
  homeUrl: '/index.html'
};

export default function App(): JSX.Element {
  return (
    <ReactRouterAppProvider navigation={NAVIGATION} branding={BRANDING}>
      <Outlet />
    </ReactRouterAppProvider>
  );
}