import { JSX } from "react";
import { Outlet } from "react-router-dom";
import { DashboardLayout } from '@toolpad/core/DashboardLayout';
import { Box } from "@mui/material";

export default function MainLayout(): JSX.Element {
  return (
    <DashboardLayout>
      <Box sx={{ flex: 1, border: 0}}>
        <Outlet />
      </Box>
    </DashboardLayout>
  );
}
