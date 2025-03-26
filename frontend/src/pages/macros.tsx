import { JSX } from 'react';
import { useNavigate } from 'react-router-dom';
import {
  Typography,
  Box,
  Button,
  Container,
  Paper
} from '@mui/material';

export default function Macros(): JSX.Element {
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
          Macros Page
        </Typography>
        <Typography variant="body1">
          This is the macros page of our application.
        </Typography>
        <Box sx={{ mt: 3 }}>
          <Button
            variant="contained"
            color="primary"
            onClick={() => navigate('/dashboard')}
          >
            Back to Dashboard
          </Button>
        </Box>
      </Paper>
    </Container>
  );
}
