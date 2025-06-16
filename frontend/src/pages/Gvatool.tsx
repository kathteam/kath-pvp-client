import { Fragment, JSX, useEffect, useState } from 'react';
import {
  Typography,
  Box,
  Snackbar,
  Alert,
  CircularProgress,
} from '@mui/material';
import {
  Plagiarism as PlagiarismIcon,
  CloudDownload as CloudDownloadIcon,
  CheckCircle as CheckCircleIcon,
  Replay as ReplayIcon,
} from '@mui/icons-material';

import DiseaseDownloadCard from '@/components/cards/DiseaseDownloadCard';
import DiseaseOptionCard from '@/components/cards/DiseaseOptionCard';
import { SetupCard } from '@/components/cards';
import { useSetupState } from '@/states/setup';
import { handleScroll } from '@/utils';
import { RouteHeader } from '@/components';
import { Button, Row } from '@/components/core';

export default function GVATool(): JSX.Element {
  const [referenceGenomePath, setReferenceGenomePath] = useState<{
		status: string;
		reference_genome: string;
	}>();

  const downloadReferenceGenome = async () => {
    setReferenceGenomePath({ status: 'processing', reference_genome: '' });
    try {
      const response =
				await window.pywebview.api.fasta_service.download_reference_genome_grch38();
      setReferenceGenomePath(response);
    } catch (error) {
      console.error(error);
      setReferenceGenomePath({ status: 'error', reference_genome: '' });
    }
  };

  const howTo = {
    media: 'DownloadFasta.gif',
    title: 'Human reference genome download',
    description: 'You can download the reference genome from the button below. The .fasta file will be present in the /.kath/shared/data/fasta_files/reference folder.'
  };

  const setupState = useSetupState();
  const [showMainContent, setShowMainContent] = useState(() => {
    if (setupState.status === 'completed') {
      return true;
    }
    return false;
  });

  useEffect(() => {
    if (setupState.status === 'completed' && !showMainContent) {
      setTimeout(() => {
        setShowMainContent(true);
      }, 1000);
    }
  }, [setupState.status, showMainContent]);

  // Reset scroll position on initial load
  useEffect(() => {
    handleScroll('GVATool');
  }, []);

  Show setup screen if not completed
  if (setupState.status !== 'completed' || !showMainContent) {
    return <SetupCard {...setupState} />;
  }

  // Show main content if setup is completed
  return (
    <Fragment>
      <RouteHeader
        icon={PlagiarismIcon}
        title="GVATool"
        description="Analyze genetic data efficiently, detect mutations, and access clear visualizations for informed decisions."
        howTo={howTo}
      />
      <Row sx={{ borderBottom: 1, borderColor: 'divider' }}>
        <Typography>
          Download the latest reference genome to ensure your data is aligned to the most current sequence. Once downloaded, you can begin your analysis.
        </Typography>
        <Box sx={{ display: 'flex', flex: 1, justifyContent: 'right' }}>
          <Button
            variant="contained"
            startIcon={
              referenceGenomePath?.status === 'success' ? (
                <CheckCircleIcon />
              ) : referenceGenomePath?.status === 'processing' ? (
                <CircularProgress size={20} color="inherit" />
              ) : referenceGenomePath?.status === 'error' ? (
                <ReplayIcon />
              ) : (
                <CloudDownloadIcon />
              )
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
        </Box>
      </Row>
      <Snackbar
        open={Boolean(referenceGenomePath)}
        autoHideDuration={6000}
        onClose={() => setReferenceGenomePath(undefined)}
        anchorOrigin={{ vertical: 'top', horizontal: 'center' }}
      >
        <Alert
          severity={
            referenceGenomePath?.status === 'success'
              ? 'success'
              : referenceGenomePath?.status === 'processing'
                ? 'warning'
                : 'error'
          }
          onClose={() => setReferenceGenomePath(undefined)}
          sx={{ width: '100%' }}
        >
          {!referenceGenomePath
            ? 'Processing...'
            : referenceGenomePath.status === 'success'
              ? `Reference genome downloaded to: ${referenceGenomePath.reference_genome}`
              : referenceGenomePath.status === 'processing'
                ? 'Downloading reference genome...'
                : 'Failed to download reference genome'}
        </Alert>
      </Snackbar>
      <DiseaseDownloadCard />
      <DiseaseOptionCard />
    </Fragment>
  );
}
