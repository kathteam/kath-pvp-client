import { useState, JSX } from 'react';
import {
  Typography,
  Box,
  Button,
  Paper,
  Alert,
  TextField,
  useTheme,
} from '@mui/material';
import DiseaseModal from '@/components/modals/DiseaseModal';

interface DiseaseOptionCardProps {
  initialFastaPath?: string;
}

export default function DiseaseOptionCard({ initialFastaPath = '' }: DiseaseOptionCardProps): JSX.Element {
  const [fastaFilePath, setFastaFilePath] = useState<string>(initialFastaPath);
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
      const diseaseData =
				await window.pywebview.api.disease_service.get_disease_data(
				  extractionResult.result_file
				);

      if (!diseaseData || !diseaseData.disease_data) {
        // TODO UNCOMMENT
        // alert('No disease data found.');
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
        setGeneticDiseaseData(diseaseData.disease_data.map((item: any) => ({
          clinicalSignificance: item.clinical_significance,
          disease: item.disease_name,
        })));
      }
      console.log(diseaseData.disease_data);
      setShowModal(true);
    } catch (error) {
      console.error('Failed to fetch disease data:', error);
      alert('Failed to fetch disease data.');
    }
  };

  return (
    <>
      <Paper
        elevation={3}
        sx={{
          p: 3,
          mt: 4,
          width: '90%',
          maxWidth: 800,
          mx: 'auto',
          textAlign: 'center',
          backgroundColor: theme.palette.background.paper,
          color: theme.palette.text.primary,
        }}
      >
        <Typography variant="h5" component="h2" gutterBottom>
					Disease Extraction
        </Typography>
        <Typography variant="body2" sx={{ mb: 2 }}>
					Specify a FASTA file path to extract disease information
        </Typography>

        <Box
          sx={{
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'center',
            mb: 2,
          }}
        >
          <TextField
            fullWidth
            label="FASTA File Path"
            variant="outlined"
            value={fastaFilePath}
            onChange={(e) => setFastaFilePath(e.target.value)}
            sx={{
              mr: 1,
              flexGrow: 1,
              input: { color: theme.palette.text.primary },
              label: { color: theme.palette.text.secondary },
            }}
          />
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

        <Button
          variant="contained"
          color="secondary"
          disabled={!extractionResult?.result_file}
          sx={{ mt: 1, ml: 2 }}
          onClick={handleDisplayGeneticDisease}
        >
					Display genetic disease information
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
