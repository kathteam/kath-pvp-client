import { JSX, useState } from 'react';
import { useNavigate } from 'react-router-dom';
import {
  Typography,
  Box,
  Button,
  Container,
  Paper,
  Snackbar,
  Stack,
  Alert,
  CircularProgress,
  Fade
} from '@mui/material';
import CloudDownloadIcon from '@mui/icons-material/CloudDownload';
import CheckCircleIcon from '@mui/icons-material/CheckCircle';
import ReplayIcon from '@mui/icons-material/Replay';
import ArrowBackIcon from '@mui/icons-material/ArrowBack';

import DiseaseDownloadCard from '@/components/cards/DiseaseDownloadCard';
import DiseaseOptionCard from '@/components/cards/DiseaseOptionCard';
// import { SetupCard } from '@/components/cards';
// import { useSetupState } from '@/states/setup';

export default function GVATool(): JSX.Element {
  const navigate = useNavigate();

  const [referenceGenomePath, setReferenceGenomePath] = useState<{
    status: string;
    reference_genome_path: string;
  }>();

  const downloadReferenceGenome = async () => {
    setReferenceGenomePath({ status: 'processing', reference_genome_path: '' });
    try {
      const response = await window.pywebview.api.fasta_service.download_reference_genome_grch38();
      setReferenceGenomePath(response);
    }
    catch (error) {
      console.error(error);
      setReferenceGenomePath({ status: 'error', reference_genome_path: '' });
    }
  };

  // const setupState = useSetupState();
  // const [showMainContent, setShowMainContent] = useState(() => {
  //   if (setupState.status === 'completed') {
  //     return true;
  //   }
  //   return false;
  // });

  // useEffect(() => {
  //   if (setupState.status === 'completed' && !showMainContent) {
  //     setTimeout(() => {
  //       setShowMainContent(true);
  //     }, 1000);
  //   }
  // }, [setupState.status, showMainContent]);

  const getButtonColor = () => {
    if (!referenceGenomePath) {
      return 'info';
    }
    if (referenceGenomePath?.status === 'processing') {
      return 'warning';
    }
    if (referenceGenomePath?.status === 'success') {
      return 'success';
    }
    return 'error';
  };

  // Show setup screen if not completed
  // if (setupState.status !== 'completed' || !showMainContent) {
  //   return <SetupCard {...setupState} />;
  // }

  // Show main content if setup is completed
  return (
    <Box
      sx={{
        minHeight: '100vh',
        py: 6,
      }}
    >
      <Fade in>
        <Container maxWidth="sm">
          <Paper
            elevation={4}
            sx={{
              p: 5,
              borderRadius: 4,
              boxShadow: 6,
              textAlign: 'center',
            }}
          >
            <Typography variant="h3" component="h1" gutterBottom fontWeight={700}>
              GVATool
            </Typography>
            <Typography variant="subtitle1" color="text.secondary" gutterBottom>
              Gene Variation Analysis Tool
            </Typography>
            <Typography variant="body1" sx={{ mb: 3 }}>
              Analyze gene variations with reference genome data. Download the latest reference genome and start your analysis.
            </Typography>
            <Stack direction="row" spacing={2} justifyContent="center" sx={{ mb: 3 }}>
              <Button
                variant="outlined"
                color="primary"
                startIcon={<ArrowBackIcon />}
                onClick={() => navigate('/dashboard')}
              >
                Dashboard
              </Button>
              <Button
                variant="contained"
                color={getButtonColor()}
                startIcon={
                  referenceGenomePath?.status === 'success' ? <CheckCircleIcon /> :
                    referenceGenomePath?.status === 'processing' ? <CircularProgress size={20} color="inherit" /> :
                      referenceGenomePath?.status === 'error' ? <ReplayIcon /> :
                        <CloudDownloadIcon />
                }
                onClick={downloadReferenceGenome}
                sx={{ minWidth: 220 }}
                disabled={referenceGenomePath?.status === 'processing'}
              >
                {referenceGenomePath?.status === 'processing'
                  ? 'Downloading...'
                  : referenceGenomePath?.status === 'success'
                    ? 'Downloaded'
                    : referenceGenomePath?.status === 'error'
                      ? 'Retry Download'
                      : 'Download Reference Genome'}
              </Button>
            </Stack>
            <Snackbar
              open={Boolean(referenceGenomePath)}
              autoHideDuration={6000}
              onClose={() => setReferenceGenomePath(undefined)}
              anchorOrigin={{ vertical: 'top', horizontal: 'center' }}
            >
              <Alert
                severity={
                  referenceGenomePath?.status === 'success' ? 'success' :
                    referenceGenomePath?.status === 'processing' ? 'warning' : 'error'
                }
                onClose={() => setReferenceGenomePath(undefined)}
                sx={{ width: '100%' }}
              >
                {!referenceGenomePath ? 'Processing...' :
                  referenceGenomePath.status === 'success'
                    ? `Reference genome downloaded to: ${referenceGenomePath.reference_genome_path}`
                    : referenceGenomePath.status === 'processing'
                      ? 'Downloading reference genome...'
                      : 'Failed to download reference genome'}
              </Alert>
            </Snackbar>
          </Paper>
          <Box sx={{ mt: 4 }}>
            <DiseaseDownloadCard />
            <DiseaseOptionCard />
          </Box>
        </Container>
      </Fade>
    </Box>
  );
}
