import { useState, JSX } from 'react';
import {
  Typography,
  Box,
  Button,
  Paper,
  Alert,
  TextField,
  useTheme,
  Stack,
  CircularProgress,
} from '@mui/material';
import ScienceIcon from '@mui/icons-material/Science';
import VisibilityIcon from '@mui/icons-material/Visibility';
import TableModal from '@/components/modals/tableModal';

export default function DiseaseOptionCard(): JSX.Element {
  const [fastaFilePath, setFastaFilePath] = useState<string>('');
  const [extractionResult, setExtractionResult] = useState<{
    status: string;
    result_file: string;
  }>();
  const [showModal, setShowModal] = useState<boolean>(false);
  const [isExtracting, setIsExtracting] = useState<boolean>(false);
  const [tableData, setTableData] = useState<any[]>([]);

  const theme = useTheme();

  const handleDiseaseExtraction = async () => {
    if (!fastaFilePath) {
      setExtractionResult({
        status: 'error',
        result_file: 'Please specify a FASTA file path',
      });
      return;
    }

    setIsExtracting(true);
    setExtractionResult({
      status: 'processing',
      result_file: 'Extracting disease information...',
    });

    try {
      const formattedPath = fastaFilePath
        .trim()
        .replace(/\\/g, '/')
        .replace(/"/g, '');

      const response =
        await window.pywebview.api.blast_service.disease_extraction(
          formattedPath
        );
      setExtractionResult({
        status: response.status || 'success',
        result_file:
          // TODO REMOVE
          'Extraction completed',
        // TODO UNCOMMENT
        // response.result_file || `Extraction completed for ${formattedPath}`,
      });
    } catch (error) {
      console.error('Disease extraction failed:', error);
      setExtractionResult({
        status: 'error',
        result_file: `Failed to extract disease information: ${error || 'Unknown error'}`,
      });
    } finally {
      setIsExtracting(false);
    }
  };

  const handleDisplayGeneticDisease = async () => {
    if (!extractionResult?.result_file) {
      alert('Please extract disease information first.');
      return;
    }

    try {
      // Fetch your table data from your backend or API
      // Example: fetch from Flask API or your backend endpoint
      const res = await fetch('http://localhost:5000/api/variants'); // Adjust URL as needed
      const data = await res.json();
      setTableData(data);
      setShowModal(true);
    } catch (error) {
      console.error('Failed to fetch table data:', error);
      alert('Failed to fetch table data.');
    }
  };

  return (
    <>
      <Paper
        elevation={6}
        sx={{
          p: 4,
          mt: 4,
          width: '95%',
          maxWidth: 800,
          mx: 'auto',
          textAlign: 'center',
          borderRadius: 4,
          boxShadow: 8,
        }}
      >
        <Typography variant="h5" component="h2" gutterBottom fontWeight={700}>
          Disease Extraction
        </Typography>
        <Typography variant="body2" sx={{ mb: 3 }}>
          Specify a FASTA file path to extract disease information.
        </Typography>

        <Box
          sx={{
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'center',
            mb: 2,
            gap: 2,
          }}
        >
          <TextField
            fullWidth
            label="FASTA File Path"
            variant="outlined"
            value={fastaFilePath}
            onChange={(e) => setFastaFilePath(e.target.value)}
            sx={{
              input: { color: theme.palette.text.primary },
              label: { color: theme.palette.text.secondary },
              borderRadius: 2,
            }}
          />
        </Box>

        <Stack direction={{ xs: 'column', sm: 'row' }} spacing={2} justifyContent="center" sx={{ mt: 2 }}>
          <Button
            variant="contained"
            color="primary"
            onClick={handleDiseaseExtraction}
            disabled={isExtracting || !fastaFilePath}
            startIcon={isExtracting ? <CircularProgress size={20} color="inherit" /> : <ScienceIcon />}
            sx={{
              fontWeight: 600,
              minWidth: 220,
              background: 'linear-gradient(90deg, #4C7380 0%, #5D8D9D 100%)',
              color: '#fff',
              '&:hover': {
                background: 'linear-gradient(90deg, #3a5a68 0%, #4c7380 100%)',
              },
            }}
          >
            {isExtracting ? 'Processing...' : 'Extract Disease Information'}
          </Button>

          <Button
            variant="contained"
            color="secondary"
            disabled={!extractionResult?.result_file}
            startIcon={<VisibilityIcon />}
            onClick={handleDisplayGeneticDisease}
            sx={{
              fontWeight: 600,
              minWidth: 220,
              background: 'linear-gradient(90deg, #5D8D9D 0%, #4C7380 100%)',
              color: '#fff',
              '&:hover': {
                background: 'linear-gradient(90deg, #4C7380 0%, #5D8D9D 100%)',
              },
            }}
          >
            Display Genetic Disease Info
          </Button>
        </Stack>

        {extractionResult && (
          <Alert
            severity={
              extractionResult.status === 'success'
                ? 'success'
                : extractionResult.status === 'processing'
                  ? 'info'
                  : 'error'
            }
            sx={{ mt: 3, fontWeight: 500 }}
            onClose={() => setExtractionResult(undefined)}
          >
            {extractionResult.result_file}
          </Alert>
        )}
      </Paper>
      {showModal && (
        <TableModal
          open={showModal}
          onClose={() => setShowModal(false)}
          data={tableData}
          title="All Variants"
        />
      )}
    </>
  );
}