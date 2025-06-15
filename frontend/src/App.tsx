import { JSX } from 'react';
import { Outlet } from 'react-router-dom';
import { ReactRouterAppProvider } from '@toolpad/core/react-router';
import { Navigation } from '@toolpad/core/AppProvider';
import { 
  Dashboard as DashboardIcon,
  Plagiarism as PlagiarismIcon,
  Inventory as InventoryIcon,
  AutoAwesome as AutoAwesomeIcon,
  MenuBook as MenuBookIcon,
  Storage as StorageIcon
} from '@mui/icons-material';
import { createTheme } from '@mui/material';
import { SetupProvider } from '@/states/setup';

const navigation: Navigation = [
  {
    segment: 'dashboard',
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
    segment: 'features/analysis_history',
    title: 'Analysis History',
    icon: <StorageIcon />,
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

const theme = createTheme({
  cssVariables: {
    colorSchemeSelector: 'data-toolpad-color-scheme',
  },
  defaultColorScheme: 'light',
  colorSchemes: {
    light: {
      palette: {
        primary: {
          main: '#1E1E24',
        },
        secondary: {
          main: '#5D8D9D',
        },
        background: {
          default: '#FFFFFB',
          paper: '#FFF9FB',
        },
        divider: '#ebecf0',
        text: {
          primary: '#404040',
          secondary: '#bdbdbd',
        },
      }
    },
    dark: {
      palette: {
        mode: 'dark',
        primary: {
          main: '#FCE4D8',
        },
        secondary: {
          main: '#5D8D9D',
        },
        background: {
          default: '#252627',
          paper: '#191716',
        },
        divider: '#2B2B2B',
        text: {
          primary: '#F6F6F6',
          secondary: '#757575',
        },
      },
    },
  },
  typography: {
    fontFamily: 'Inter',
  },
});

export default function App(): JSX.Element {
  return (
    <SetupProvider>
      <ReactRouterAppProvider navigation={navigation} theme={theme}>
        <Outlet />
      </ReactRouterAppProvider>
    </SetupProvider>
  );
}
