import { JSX, useState } from 'react';
import { useTheme } from '@mui/material/styles';
import {
  Container,
  Typography,
  Card,
  CardContent,
  Button,
  Grid,
  Box,
  Divider,
  Chip,
  CircularProgress,
  Paper,
//   Alert
} from '@mui/material';

export default function BlastCard(): JSX.Element {
  const theme = useTheme();

  // State for align mutations operation
  const [alignResult, setAlignResult] = useState<{
    status: string;
    result_file: string;
  }>();

  // State for blast analysis operation
  const [blastResult, setBlastResult] = useState<{
    status: string;
    result_file: string;
  }>();

  // Loading states
  const [alignLoading, setAlignLoading] = useState<boolean>(false);
  const [blastLoading, setBlastLoading] = useState<boolean>(false);

  // Error states
  const [alignError, setAlignError] = useState<string>('');
  const [blastError, setBlastError] = useState<string>('');

  // Function to handle align mutations
  const handleAlignMutations = async () => {
    setAlignLoading(true);
    setAlignError('');
    
    try {
      const response = await window.pywebview.api.blast_service.align_mutations();
      setAlignResult(response);
    } catch (err) {
      setAlignError('Failed to align mutations. Please try again.');
      console.error('Error aligning mutations:', err);
    } finally {
      setAlignLoading(false);
    }
  };

  // Function to handle blast analysis
  const handleBlastAnalysis = async () => {
    setBlastLoading(true);
    setBlastError('');
    
    try {
      const response = await window.pywebview.api.blast_service.perform_blast_analysis();
      setBlastResult(response);
    } catch (err) {
      setBlastError('Failed to perform BLAST analysis. Please try again.');
      console.error('Error performing BLAST analysis:', err);
    } finally {
      setBlastLoading(false);
    }
  };

  // Helper function to determine button color
  const getButtonColor = (result?: { status: string }) => {
    if (!result) {
      return 'primary';
    }
    if (result.status === 'processing') {
      return 'warning';
    }
    if (result.status === 'completed' || result.status === 'success') {
      return 'success';
    }
    return 'error';
  };

  // Helper function to get button text
  const getButtonText = (result: { status: string } | undefined, loading: boolean, operation: string) => {
    if (loading) {
      return 'Processing...';
    }
    if (!result) {
      return operation;
    }
    if (result.status === 'completed' || result.status === 'success') {
      return 'Completed';
    }
    if (result.status === 'processing') {
      return 'Processing...';
    }
    return 'Failed';
  };

  return (
    <Container maxWidth="md" sx={{ py: 4 }}>
      <Typography variant="h4" component="h1" gutterBottom sx={{ mb: 3 }}>
        BLAST Analysis Tools
      </Typography>
      <Typography variant="body1" component="h2" gutterBottom sx={{ mb: 3 }}>
        Caution: files are taken only from 'uploads' folder.
      </Typography>
      
      <Grid container spacing={3}>
        {/* Align Mutations Card */}
        <Grid item xs={12} md={6}>
          <Card elevation={3} sx={{ height: '100%' }}>
            <CardContent>
              <Typography variant="h6" gutterBottom>
                Align Mutations
              </Typography>
              <Typography variant="body2" color="text.secondary" sx={{ mb: 3 }}>
                Align gene mutations with the reference genome
              </Typography>
              
              <Button
                variant="contained"
                color={getButtonColor(alignResult)}
                onClick={handleAlignMutations}
                disabled={alignLoading}
                fullWidth
                sx={{ mb: 2 }}
              >
                {alignLoading ? (
                  <CircularProgress size={24} color="inherit" />
                ) : (
                  getButtonText(alignResult, alignLoading, 'Align Mutations')
                )}
              </Button>
              
              {alignError && (
                <Paper 
                  sx={{ 
                    p: 2, 
                    mt: 2, 
                    bgcolor: theme.palette.mode === 'dark' ? 'rgba(255, 0, 0, 0.1)' : '#fdeded' 
                  }}
                >
                  <Typography color="error">{alignError}</Typography>
                </Paper>
              )}
              
              {alignResult && (
                <Box sx={{ mt: 2 }}>
                  <Divider sx={{ my: 1 }} />
                  <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mt: 1 }}>
                    <Typography variant="body2" color="text.secondary">Status:</Typography>
                    <Chip 
                      label={alignResult.status}
                      color={alignResult.status === 'completed' || alignResult.status === 'success' ? 'success' : 'primary'}
                      size="small"
                    />
                  </Box>
                  <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mt: 1 }}>
                    <Typography variant="body2" color="text.secondary">Result File:</Typography>
                    <Typography variant="body2" noWrap sx={{ maxWidth: '70%', textAlign: 'right' }}>
                      {alignResult.result_file}
                    </Typography>
                  </Box>
                </Box>
              )}
            </CardContent>
          </Card>
        </Grid>

        {/* BLAST Analysis Card */}
        <Grid item xs={12} md={6}>
          <Card elevation={3} sx={{ height: '100%' }}>
            <CardContent>
              <Typography variant="h6" gutterBottom>
                BLAST Analysis
              </Typography>
              <Typography variant="body2" color="text.secondary" sx={{ mb: 3 }}>
                Perform BLAST analysis on sequence data
              </Typography>
              
              <Button
                variant="contained"
                color={getButtonColor(blastResult)}
                onClick={handleBlastAnalysis}
                disabled={blastLoading}
                fullWidth
                sx={{ mb: 2 }}
              >
                {blastLoading ? (
                  <CircularProgress size={24} color="inherit" />
                ) : (
                  getButtonText(blastResult, blastLoading, 'Perform BLAST Analysis')
                )}
              </Button>
              
              {blastError && (
                <Paper 
                  sx={{ 
                    p: 2, 
                    mt: 2, 
                    bgcolor: theme.palette.mode === 'dark' ? 'rgba(255, 0, 0, 0.1)' : '#fdeded' 
                  }}
                >
                  <Typography color="error">{blastError}</Typography>
                </Paper>
              )}
              
              {blastResult && (
                <Box sx={{ mt: 2 }}>
                  <Divider sx={{ my: 1 }} />
                  <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mt: 1 }}>
                    <Typography variant="body2" color="text.secondary">Status:</Typography>
                    <Chip 
                      label={blastResult.status}
                      color={blastResult.status === 'completed' || blastResult.status === 'success' ? 'success' : 'primary'}
                      size="small"
                    />
                  </Box>
                  <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mt: 1 }}>
                    <Typography variant="body2" color="text.secondary">Result File:</Typography>
                    <Typography variant="body2" noWrap sx={{ maxWidth: '70%', textAlign: 'right' }}>
                      {blastResult.result_file}
                    </Typography>
                  </Box>
                </Box>
              )}
            </CardContent>
          </Card>
        </Grid>
      </Grid>
    </Container>
  );
}
