import { JSX } from 'react';
import { useNavigate } from 'react-router-dom';
import {
  Typography,
  Box,
  Button,
  Container,
  Stack,
  Paper
} from '@mui/material';
import DashboardIcon from '@mui/icons-material/Dashboard';
import PlagiarismIcon from '@mui/icons-material/Plagiarism';
import InventoryIcon from '@mui/icons-material/Inventory';
import AutoAwesomeIcon from '@mui/icons-material/AutoAwesome';
import MenuBookIcon from '@mui/icons-material/MenuBook';
import StorageIcon from '@mui/icons-material/Storage';

export default function Dashboard(): JSX.Element {
  const navigate = useNavigate();

  return (
    <Box
      sx={{
        minHeight: '100vh',
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'center',
        py: 4,
      }}
    >
      <Container maxWidth="sm">
        <Paper
          elevation={6}
          sx={{
            p: 5,
            borderRadius: 4,
            textAlign: 'center',
            boxShadow: 6,
          }}
        >
          <DashboardIcon color="primary" sx={{ fontSize: 48, mb: 1 }} />
          <Typography variant="h4" component="h1" gutterBottom fontWeight={700}>
            Welcome to the Dashboard
          </Typography>
          <Typography variant="body1" color="text.secondary" sx={{ mb: 3 }}>
            Access all features and resources from here.
          </Typography>
          <Stack
            direction="column"
            spacing={2}
            alignItems="center"
            sx={{ mt: 2 }}
          >
            <Button
              variant="outlined"
              startIcon={<PlagiarismIcon />}
              fullWidth
              size="large"
              onClick={() => navigate('/features/gvatool')}
            >
              GVATool
            </Button>
            <Button
              variant="outlined"
              startIcon={<StorageIcon />}
              fullWidth
              size="large"
              onClick={() => navigate('/features/analysis_history')}
            >
              Analysis History
            </Button>
            <Button
              variant="outlined"
              startIcon={<InventoryIcon />}
              fullWidth
              size="large"
              onClick={() => navigate('/system/file_manager')}
            >
              File Manager
            </Button>
            <Button
              variant="outlined"
              startIcon={<AutoAwesomeIcon />}
              fullWidth
              size="large"
              onClick={() => navigate('/system/macros')}
            >
              Macros
            </Button>
            <Button
              variant="outlined"
              startIcon={<MenuBookIcon />}
              fullWidth
              size="large"
              onClick={() => navigate('/resources/manual')}
            >
              Manual
            </Button>
          </Stack>
        </Paper>
      </Container>
    </Box>
  );
}
