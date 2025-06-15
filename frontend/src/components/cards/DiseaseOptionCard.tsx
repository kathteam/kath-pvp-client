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
import DiseaseModal from '@/components/modals/DiseaseModal';
import ScienceIcon from '@mui/icons-material/Science';
import VisibilityIcon from '@mui/icons-material/Visibility';
import HowTo from '../buttons/HowTo';

export default function DiseaseOptionCard(): JSX.Element {
  const [fastaFilePath, setFastaFilePath] = useState<string>('');
  const [extractionResult, setExtractionResult] = useState<{
    status: string;
    result_file: string;
  }>();
  const [geneticDiseaseData, setGeneticDiseaseData] = useState<{
    clinicalSignificance: string;
    disease: string;
  }[]>([]);

  const [showModal, setShowModal] = useState<boolean>(false);
  const [isExtracting, setIsExtracting] = useState<boolean>(false);

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
        result_file: response.result_file || `Extraction completed for ${formattedPath}`,
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
      const diseaseData =
        await window.pywebview.api.disease_service.get_disease_data(
          extractionResult.result_file
        );

      if (!diseaseData || diseaseData.length === 0) {
        // TODO UNCOMMENT
        alert('No disease data found.');
        setGeneticDiseaseData([
          {
            clinicalSignificance: 'Pathogenic',
            disease: 'Type 1 Diabetes Mellitus',
          },
          {
            clinicalSignificance: 'Likely Pathogenic',
            disease: 'Maturity-Onset Diabetes of the Young (MODY)',
          },
          {
            clinicalSignificance: 'Benign',
            disease: 'Gestational Diabetes',
          },
          {
            clinicalSignificance: 'Uncertain Significance',
            disease: 'Latent Autoimmune Diabetes in Adults (LADA)',
          },
          {
            clinicalSignificance: 'Likely Benign',
            disease: 'Prediabetes',
          },
        ]);
      } else {
        setGeneticDiseaseData(diseaseData.map((item: any) => ({
          clinicalSignificance: item.clinical_significance,
          disease: item.disease_name,
        })));
      }
      setShowModal(true);
    } catch (error) {
      console.error('Failed to fetch disease data:', error);
      alert('Failed to fetch disease data.');
    }
  };

  const howTo = 
  {
    media: 'ExtractString.gif',
    title: 'Extract genetic diseases from found variations',
    description: 'Select the .fasta file and copy it as a string (You are free to use our File Manager as a substitute for this step), with the analyze button you will receive a list of genetic diseases that are associated with the variations found. WARNING: File may not have any variations or diseases.'
  };

  return (
    <>
      <Paper
        elevation={4}
        sx={{
          p: 5,
          textAlign: 'center',
          borderRadius: 4,
          boxShadow: 6,
          position: 'relative'
        }}
      >
        <Box sx={{ position: 'absolute', top: 5, right: 5 }}>
          <HowTo media={howTo.media} title={howTo.title} description={howTo.description} />
        </Box>
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
        <DiseaseModal
          diseases={geneticDiseaseData}
          onClose={() => setShowModal(false)}
        />
      )}
    </>
  );
}
