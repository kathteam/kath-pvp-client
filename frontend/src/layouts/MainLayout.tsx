import { JSX } from 'react';
import { Outlet } from 'react-router-dom';
import { DashboardLayout, ThemeSwitcher } from '@toolpad/core/DashboardLayout';
import { Box, IconButton, Stack } from '@mui/material';
import { LogoIcon, TitleIcon } from '@/components/icons';
import FullscreenIcon from '@mui/icons-material/Fullscreen';

function customAppTitle(): JSX.Element {
  return (
    <Stack direction='row' alignItems='center' spacing={2} px={0.5}>
      <LogoIcon />
      <TitleIcon />
    </Stack>
  );
}

function customToolbarActions(): JSX.Element {
  return (
    <Stack direction='row' alignItems='center' spacing={2} px={0.5}>
      <ThemeSwitcher />
      <IconButton
        color="primary"
        onClick={() => window.pywebview.api.ui_controller.fullscreen()}
        sx={{ mb: 4 }}
      >
        <FullscreenIcon />
      </IconButton>
    </Stack>
  );
}

export default function MainLayout(): JSX.Element {
  return (
    <DashboardLayout slots={{
      appTitle: customAppTitle,
      toolbarActions: customToolbarActions,
    }}>
      <Box sx={{ flex: 1, border: 0 }}>
        <Outlet />
      </Box>
    </DashboardLayout>
  );
}
