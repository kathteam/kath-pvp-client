import { SetupContextType } from '@/states/setup';
import { Button, CircularProgress, LinearProgress, Paper, Typography } from '@mui/material';
import { JSX } from 'react';

export default function SetupCard({
  status,
  progress,
  message
}: SetupContextType): JSX.Element {
  return (
    <Paper elevation={0} sx={{ p: 3, width: '100%', height: '100%', textAlign: 'center', alignContent: 'center' }}>
      <Typography variant="h5" mb={3}>
        Setting Up Application
      </Typography>
      
      {(status === 'started' || status === 'progress') && (
        <>
          <LinearProgress 
            variant="determinate" 
            value={progress * 100} 
            sx={{ width: '80%', maxWidth: 400, mb: 2 }} 
          />
          <Typography variant="body1" mb={4}>
            {message}
          </Typography>
          <CircularProgress size={40} />
        </>
      )}
      
      {status === 'idle' && (
        <>
          <CircularProgress size={40} />
          <Typography variant="body1" mt={2}>
            Initializing...
          </Typography>
        </>
      )}
      
      {status === 'failed' && (
        <>
          <Typography color="error" variant="h6" mb={2}>
            Setup Failed
          </Typography>
          <Typography color="error" variant="body2" sx={{ maxWidth: 500 }}>
            {message}
          </Typography>
          <Button 
            variant="contained" 
            color="primary" 
            sx={{ mt: 3 }}
            onClick={() => window.location.reload()}
          >
            Retry
          </Button>
        </>
      )}
    </Paper>
  );
}
