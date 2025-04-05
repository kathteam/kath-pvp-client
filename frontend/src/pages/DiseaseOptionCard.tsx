import { useState, JSX } from 'react';
import {
  Typography,
  Box,
  Button,
  Paper,
  Alert,
  TextField,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  useTheme
} from '@mui/material';

export default function DiseaseOptionCard(): JSX.Element {
  const [fastaFilePath, setFastaFilePath] = useState<string>('');
  const [extractionResult, setExtractionResult] = useState<{
    status: string;
    message: string;
  }>();
  const [isExtracting, setIsExtracting] = useState<boolean>(false);
  const [isDialogOpen, setIsDialogOpen] = useState<boolean>(false);

  const theme = useTheme(); // Access Material-UI theme for dark mode compatibility

  const handleDiseaseExtraction = async () => {
    if (!fastaFilePath) {
      setExtractionResult({
        status: 'error',
        message: 'Please specify a FASTA file path'
      });
      return;
    }

    setIsExtracting(true);
    setExtractionResult({
      status: 'processing',
      message: 'Extracting disease information...'
    });

    try {
      // Ensure the path is properly formatted
      const formattedPath = fastaFilePath.trim();

      const response = await window.pywebview.api.blast_service.disease_extraction(formattedPath);
      setExtractionResult({
        status: response.status || 'success',
        message: response.result_file || `Extraction completed for ${formattedPath}`
      });
    } catch (error) {
      console.error('Disease extraction failed:', error);
      setExtractionResult({
        status: 'error',
        message: `Failed to extract disease information: ${error || 'Unknown error'}`
      });
    } finally {
      setIsExtracting(false);
    }
  };

  const openDialog = () => {
    setIsDialogOpen(true);
  };

  const closeDialog = () => {
    setIsDialogOpen(false);
  };

  const handleDialogSubmit = () => {
    if (fastaFilePath) {
      closeDialog();
    } else {
      setExtractionResult({
        status: 'error',
        message: 'Please specify a valid folder path.'
      });
    }
  };

  return (
    <Paper
      elevation={3}
      sx={{
        p: 3,
        mt: 4,
        width: '90%',
        maxWidth: 800,
        mx: 'auto',
        textAlign: 'center',
        backgroundColor: theme.palette.background.paper, // Theme-aware background
        color: theme.palette.text.primary // Theme-aware text color
      }}
    >
      <Typography variant="h5" component="h2" gutterBottom>
        Disease Extraction
      </Typography>
      <Typography variant="body2" sx={{ mb: 2 }}>
        Specify a FASTA file path to extract disease information
      </Typography>

      <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'center', mb: 2 }}>
        <TextField
          fullWidth
          label="FASTA File Path"
          variant="outlined"
          value={fastaFilePath}
          onChange={(e) => setFastaFilePath(e.target.value)}
          sx={{
            mr: 1,
            flexGrow: 1,
            input: { color: theme.palette.text.primary }, // Theme-aware input text
            label: { color: theme.palette.text.secondary } // Theme-aware label
          }}
        />
        <Button
          variant="contained"
          onClick={openDialog}
          sx={{ flexShrink: 0 }}
        >
          Browse
        </Button>
      </Box>

      <Button
        variant="contained"
        color="primary"
        onClick={handleDiseaseExtraction}
        disabled={isExtracting || !fastaFilePath}
        sx={{ mt: 1 }}
      >
        {isExtracting ? 'Processing...' : 'Extract Disease Information'}
      </Button>

      {extractionResult && (
        <Alert
          severity={
            extractionResult.status === 'success'
              ? 'success'
              : extractionResult.status === 'processing'
              ? 'info'
              : 'error'
          }
          sx={{ mt: 2 }}
          onClose={() => setExtractionResult(undefined)}
        >
          {extractionResult.message}
        </Alert>
      )}

      {/* Dialog for folder selection */}
      <Dialog open={isDialogOpen} onClose={closeDialog}>
        <DialogTitle>Select Folder Path</DialogTitle>
        <DialogContent>
          <TextField
            fullWidth
            label="Folder Path"
            variant="outlined"
            value={fastaFilePath}
            onChange={(e) => setFastaFilePath(e.target.value)}
            sx={{
              mt: 2,
              input: { color: theme.palette.text.primary }, // Theme-aware input text
              label: { color: theme.palette.text.secondary } // Theme-aware label
            }}
          />
        </DialogContent>
        <DialogActions>
          <Button onClick={closeDialog} color="secondary">
            Cancel
          </Button>
          <Button onClick={handleDialogSubmit} color="primary">
            Submit
          </Button>
        </DialogActions>
      </Dialog>
    </Paper>
  );
}