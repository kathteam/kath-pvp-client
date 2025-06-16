import { Fragment, JSX, useState } from 'react';
import {
  TextField,
  Typography,
  CircularProgress,
  Chip,
  useTheme,
  Box,
} from '@mui/material';
import SearchIcon from '@mui/icons-material/Search';
import ErrorOutlineIcon from '@mui/icons-material/ErrorOutline';
import { Button, Column, Row } from '@/components/core';
import { RouteHeader } from '@/components';

export default function DiseaseDownloadCard(): JSX.Element {
  const [queryParams, setQueryParams] = useState<{
    disease: string;
    ref_max: number;
  }>({ disease: '', ref_max: 10 });
  const theme = useTheme();

  const [diseaseData, setDiseaseData] = useState<{
    status: string;
    disease_term: string;
    max_results: number;
    downloaded_files: string[];
    count: number;
  }>();

  const [loading, setLoading] = useState<boolean>(false);
  const [error, setError] = useState<string>('');

  const handleInputChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const { name, value } = e.target;
    setQueryParams({
      ...queryParams,
      [name]: name === 'ref_max' ? Number(value) : value,
    });
  };

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setLoading(true);
    setError('');

    try {
      const response = await window.pywebview.api.fasta_service.create_disease_download(
        queryParams.disease,
        queryParams.ref_max
      );
      setDiseaseData(response);
    } catch (err) {
      setError('Failed to fetch disease data. Please try again.');
      console.error('Error fetching disease data:', err);
    } finally {
      setLoading(false);
    }
  };

  const getStatusColor = (status: string) => {
    if (status === 'completed') return 'success';
    if (status === 'processing') return 'warning';
    if (status === 'error') return 'error';
    return 'default';
  };

  const howTo = {
    media: 'GeneDownload.gif',
    title: 'Download disease fasta files',
    description: 'Download specific disease fasta files for testing or analysis purposes. You can select how many references from found variations you want to download.'
  };

  return (
    <Fragment>
      <RouteHeader
        title="Disease Query Tool"
        description="Search for references related to a specific disease. Results are in FASTA format."
        howTo={howTo}
        sx={true}
      />
      
      <form onSubmit={handleSubmit}>
        <Row >
          <TextField
            fullWidth
            label="Disease Name"
            name="disease"
            value={queryParams.disease}
            onChange={handleInputChange}
            required
            helperText="e.g. 'diabetes', 'alzheimer'"
            variant="outlined"
            sx={{ flex: 4 }}
          />
          <TextField
            fullWidth
            label="Max References"
            name="ref_max"
            type="number"
            value={queryParams.ref_max}
            onChange={handleInputChange}
            required
            helperText="Maximum number of references to download"
            variant="outlined"
            sx={{ flex: 1 }}
          />
        </Row>
        <Column sx={{ pt: 0, borderBottom: 1, borderColor: 'divider' }}>
          <Button
            variant="contained"
            type="submit"
            fullWidth
            size="large"
            startIcon={loading ? <CircularProgress size={20} color="inherit" /> : <SearchIcon />}
            disabled={loading || !queryParams.disease}
            sx={{
              fontWeight: 600,
              maxWidth: '50%'
            }}
          >
            {loading ? 'Searching...' : 'Search'}
          </Button>
        </Column>
      </form>
      
      {error && (
        <Row sx={{ pt: 0, borderBottom: 1, borderColor: 'divider' }}>
          <Row sx={{ flex: 1, p: theme.spacing(2), bgcolor: 'background.default', borderRadius: theme.shape.borderRadius }}>
            <ErrorOutlineIcon color="error" />
            <Typography color="error">{error}</Typography>
          </Row>
        </Row>
      )}

      {diseaseData && (
        <Fragment>
          <Row sx={{ pb: 0 }}>
            <Typography variant="h6" gutterBottom sx={{ fontWeight: 600 }}>
              Results for "{diseaseData.disease_term}"
            </Typography>
          </Row>
          <Column sx={{ borderBottom: 1, borderColor: 'divider' }}>
            <Row sx={{ p: 0, width: 1 }}>
              <Typography variant="body2" color="text.secondary" sx={{ flex: 1 }}>
                Status
              </Typography>
              <Typography variant="body2" color="text.secondary" sx={{ flex: 1 }}>
                Max Results
              </Typography>
              <Typography variant="body2" color="text.secondary" sx={{ flex: 1 }}>
                References Found
              </Typography>
            </Row>
            <Row sx={{ p: 0, width: 1 }}>
              <Box sx={{ flex: 1 }}>
                <Chip
                  label={diseaseData.status}
                  color={getStatusColor(diseaseData.status)}
                  size="small"
                  sx={{ fontWeight: 600, textTransform: 'capitalize' }}
                />
              </Box>
              <Typography variant="body1" sx={{ flex: 1 }}>{diseaseData.max_results}</Typography>
              <Typography variant="body1" sx={{ flex: 1 }}>{diseaseData.count}</Typography>
            </Row>
          </Column>
        </Fragment>
      )}
    </Fragment>
  );
}
