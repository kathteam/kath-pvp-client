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

export default function Dashboard(): JSX.Element {
  const navigate = useNavigate();

  return (
    <Container
      maxWidth="md"
      sx={{
        display: 'flex',
        justifyContent: 'center',
        py: 4
      }}
    >
      <Paper elevation={0} sx={{ p: 3, width: '100%', textAlign: 'center' }}>
        <Typography variant="h4" component="h1" gutterBottom>
          Dashboard Page
        </Typography>
        <Typography variant="body1">
          This is the dashboard page of our application.
        </Typography>
        <Button
          variant="contained"
          color="primary"
          onClick={() => window.pywebview.api.fullscreen()}
          sx={{ mt: 2 }}
        >
          Test fullscreen
        </Button>
        <Box sx={{ p: 2, mt: 3 }}>
          <Stack
            direction="row"
            spacing={2}
            justifyContent="center"
            flexWrap="wrap"
          >
            <Button
              variant="outlined"
              onClick={() => navigate('/features/gvatool')}
            >
              GVATool
            </Button>
            <Button
              variant="outlined"
              onClick={() => navigate('/system/file_manager')}
            >
              File Manager
            </Button>
            <Button
              variant="outlined"
              onClick={() => navigate('/system/macros')}
            >
              Macros
            </Button>
            <Button
              variant="outlined"
              onClick={() => navigate('/resources/manual')}
            >
              Manual
            </Button>
          </Stack>
        </Box>
      </Paper>
    </Container>
  );
}
