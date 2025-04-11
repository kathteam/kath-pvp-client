import { JSX, useEffect, useState } from 'react';
import { useNavigate } from 'react-router-dom';
import {
  Typography,
  Box,
  Button,
  Container,
  Paper,
  Snackbar,
  Stack,
  Alert
} from '@mui/material';

import DiseaseDownloadCard from '@/components/cards/DiseaseDownloadCard';
import DiseaseOptionCard from '@/components/cards/DiseaseOptionCard';
import { SetupCard } from '@/components/cards';
import { useSetupState } from '@/states/setup';

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

  const setupState = useSetupState();
  const [showMainContent, setShowMainContent] = useState(() => {
    if (setupState.status === 'completed') {
      return true;
    }
    return false;
  });
  // const [playAnimation, setPlayAnimation] = useState(false);

  useEffect(() => {
    if (setupState.status === 'completed' && !showMainContent) {
      setTimeout(() => {
        setShowMainContent(true);
      }, 1500);
    }
  }, [setupState.status, showMainContent]);

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
  if (setupState.status !== 'completed' || !showMainContent) {
    return <SetupCard {...setupState} />;
  }

  // Show main content if setup is completed
  return (
    <>
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
          GVATool
          </Typography>
          <Typography variant="body1">
          This is the gene variation analysis tool page of our application.
          </Typography>
          <Box sx={{ mt: 3 }}>
            <Button
              variant="contained"
              color="primary"
              onClick={() => navigate('/dashboard')}
            >
            Back to Dashboard
            </Button>
            <Box sx={{ mt: 3, display: 'flex', gap: 2, justifyContent: 'center' }}>
              <Button
                variant="contained"
                color={getButtonColor()}
                onClick={() => downloadReferenceGenome()}
                sx={{ mt: 2 }}
                disabled={referenceGenomePath?.status === 'processing'}
              >
                {referenceGenomePath?.status === 'processing'
                  ? 'Downloading...'
                  : referenceGenomePath?.status === 'success'
                    ? 'Downloaded'
                    : 'Download reference Genome'}
              </Button>
              <Button
                variant="contained"
                color="info"
                onClick={() => NaN}
                sx={{ mt: 2 }}
              >
              Test perform aligning
              </Button>
            </Box>
          </Box>
        </Paper>
        <Stack>
          <Snackbar
            open={Boolean(referenceGenomePath)}
            autoHideDuration={6000}
            onClose={() => setReferenceGenomePath(undefined)}
            message={referenceGenomePath?.status}
          >
            <Alert 
              severity={
                referenceGenomePath?.status === 'success' ? 'success' : 
                  referenceGenomePath?.status === 'processing' ? 'warning' : 'error'
              }
              onClose={() => setReferenceGenomePath(undefined)}
            >
              {!referenceGenomePath ? 'Processing...' :
                referenceGenomePath.status === 'success' 
                  ? `Reference genome downloaded to: ${referenceGenomePath.reference_genome_path}` 
                  : referenceGenomePath.status === 'processing'
                    ? 'Downloading reference genome...'
                    : 'Failed to download reference genome'}
            </Alert>
          </Snackbar>
        </Stack>
      </Container>
      <DiseaseDownloadCard/>
      <DiseaseOptionCard/>
    </>
  );
}
