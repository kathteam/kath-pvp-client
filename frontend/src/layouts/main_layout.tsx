import { JSX } from "react";
import { Outlet } from "react-router-dom";
import { DashboardLayout } from '@toolpad/core/DashboardLayout';
import { Box, Stack } from "@mui/material";
import { LogoIcon, TitleIcon} from "@/components/icons";

function customAppTitle(): JSX.Element {
  return (
    <Stack direction='row' alignItems='center' spacing={2} px={0.5}>
      <LogoIcon />
      <TitleIcon />
    </Stack>
  )
}

export default function MainLayout(): JSX.Element {
  return (
    <DashboardLayout slots={{
      appTitle: customAppTitle,
    }}>
      <Box sx={{ flex: 1, border: 0}}>
        <Outlet />
      </Box>
    </DashboardLayout>
  );
}
