import { JSX, useState } from 'react';
import {
  Container,
  TextField,
  Button,
  Typography,
  Card,
  CardContent,
  Grid,
  CircularProgress,
  Divider,
  Paper,
  Chip,
  Box,
  Fade,
} from '@mui/material';
import SearchIcon from '@mui/icons-material/Search';
import ErrorOutlineIcon from '@mui/icons-material/ErrorOutline';
import HowTo from '../buttons/HowTo';

export default function DiseaseDownloadCard(): JSX.Element {
  const [queryParams, setQueryParams] = useState<{
    disease: string;
    ref_max: number;
  }>({ disease: '', ref_max: 10 });

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

  const howTo = 
  {
    media: 'DownloadFasta.gif',
    title: 'Download disease fasta files',
    description: 'Download specific disease fasta files for testing or analysis purposes. You can select how many references from found variations you want to download.'
  };

  return (
    <Container maxWidth="md" sx={{ py: 4 }}>
      <Fade in>
        <Card
          elevation={4}
          sx={{
            p: 5,
            borderRadius: 4,
            boxShadow: 6,
            position: 'relative'
          }}
        >
          <CardContent>
            <Box sx={{ position: 'absolute', top: 5, right: 5 }}>
              <HowTo media={howTo.media} title={howTo.title} description={howTo.description} />
            </Box>
            <Typography
              variant="h4"
              component="h1"
              gutterBottom
              sx={{ mb: 2, fontWeight: 700, textAlign: 'center' }}
            >
              Disease Query Tool
            </Typography>
            <Typography
              variant="body1"
              color="text.secondary"
              gutterBottom
              sx={{ textAlign: 'center', mb: 3 }}
            >
              Search for references related to a specific disease. Results are in FASTA format.
            </Typography>
            <Divider sx={{ mb: 3 }} />
            <form onSubmit={handleSubmit}>
              <Grid container spacing={3} alignItems="center">
                <Grid item xs={12} sm={8}>
                  <TextField
                    fullWidth
                    label="Disease Name"
                    name="disease"
                    value={queryParams.disease}
                    onChange={handleInputChange}
                    required
                    helperText="e.g. 'diabetes', 'alzheimer'"
                    variant="outlined"
                  />
                </Grid>
                <Grid item xs={12} sm={4}>
                  <TextField
                    fullWidth
                    label="Max References"
                    name="ref_max"
                    type="number"
                    value={queryParams.ref_max}
                    onChange={handleInputChange}
                    required
                    inputProps={{ min: 1 }}
                    variant="outlined"
                  />
                </Grid>
                <Grid item xs={12}>
                  <Button
                    variant="contained"
                    type="submit"
                    fullWidth
                    size="large"
                    startIcon={loading ? <CircularProgress size={20} color="inherit" /> : <SearchIcon />}
                    disabled={loading || !queryParams.disease}
                    sx={{
                      fontWeight: 600,
                      background: 'linear-gradient(90deg, #4C7380 0%, #5D8D9D 100%)',
                      color: '#fff',
                      '&:hover': {
                        background: 'linear-gradient(90deg, #3a5a68 0%, #4c7380 100%)',
                      },
                    }}
                  >
                    {loading ? 'Searching...' : 'Search'}
                  </Button>
                </Grid>
              </Grid>
            </form>
            {error && (
              <Paper
                sx={{
                  p: 2,
                  mt: 3,
                  mb: 1,
                  bgcolor: '#fdeded',
                  display: 'flex',
                  alignItems: 'center',
                  gap: 1,
                }}
                elevation={0}
              >
                <ErrorOutlineIcon color="error" />
                <Typography color="error">{error}</Typography>
              </Paper>
            )}
          </CardContent>
        </Card>
      </Fade>

      {diseaseData && (
        <Fade in>
          <Card elevation={4} sx={{ borderRadius: 4, mb: 2 }}>
            <CardContent>
              <Typography variant="h6" gutterBottom sx={{ fontWeight: 600 }}>
                Results for "{diseaseData.disease_term}"
              </Typography>
              <Divider sx={{ my: 2 }} />
              <Box sx={{ mb: 2 }}>
                <Grid container spacing={2}>
                  <Grid item xs={12} sm={4}>
                    <Typography variant="body2" color="text.secondary">
                      Status
                    </Typography>
                    <Chip
                      label={diseaseData.status}
                      color={getStatusColor(diseaseData.status)}
                      size="small"
                      sx={{ fontWeight: 600, textTransform: 'capitalize' }}
                    />
                  </Grid>
                  <Grid item xs={6} sm={4}>
                    <Typography variant="body2" color="text.secondary">
                      Max Results
                    </Typography>
                    <Typography variant="body1">{diseaseData.max_results}</Typography>
                  </Grid>
                  <Grid item xs={6} sm={4}>
                    <Typography variant="body2" color="text.secondary">
                      References Found
                    </Typography>
                    <Typography variant="body1">{diseaseData.count}</Typography>
                  </Grid>
                </Grid>
              </Box>
              {/* You can add a list of downloaded files or further actions here */}
            </CardContent>
          </Card>
        </Fade>
      )}
    </Container>
  );
}
