import { JSX, useState } from "react";
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
  Box
} from "@mui/material";

export default function DiseaseCard(): JSX.Element {
  const [queryParams, setQueryParams] = useState<{
    disease: string;
    ref_max: number;
  }>({ disease: "", ref_max: 10 });
  
  const [diseaseData, setDiseaseData] = useState<{
    status: string;
    disease_term: string;
    max_results: number;
    downloaded_files: string[];
    count: number;
  }>();
  
//   const [disease, setDisease] = useState<string[]>([]);
  const [loading, setLoading] = useState<boolean>(false);
  const [error, setError] = useState<string>("");

  const handleInputChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const { name, value } = e.target;
    setQueryParams({
      ...queryParams,
      [name]: name === "ref_max" ? Number(value) : value,
    });
  };

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setLoading(true);
    setError("");
    
    try {
      // Replace with your actual API endpoint
    const response = await window.pywebview.api.fasta_service.create_disease_download(queryParams.disease, queryParams.ref_max);
      setDiseaseData(response);
    //   setDisease(response.downloaded_files || []);
    } catch (err) {
      setError("Failed to fetch disease data. Please try again.");
      console.error("Error fetching disease data:", err);
    } finally {
      setLoading(false);
    }
  };

  return (
    <Container maxWidth="md" sx={{ py: 4 }}>
      <Typography variant="h4" component="h1" gutterBottom sx={{ mb: 3 }}>
        Disease Query Tool
      </Typography>
      <Typography variant="body1" color="text.secondary" gutterBottom>
        Use this tool to search for references related to a specific disease. Results in FASTA format.
      </Typography>
      
      <Card elevation={3} sx={{ mb: 4 }}>
        <CardContent>
          <Typography variant="h6" gutterBottom>
            Search Parameters
          </Typography>
          
          <form onSubmit={handleSubmit}>
            <Grid container spacing={3}>
              <Grid item xs={12} sm={8}>
                <TextField
                  fullWidth
                  label="Disease Name"
                  name="disease"
                  value={queryParams.disease}
                  onChange={handleInputChange}
                  required
                  helperText="Enter a disease name (e.g. 'diabetes', 'alzheimer')"
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
                  color="primary" 
                  type="submit"
                  disabled={loading || !queryParams.disease}
                  fullWidth
                >
                  {loading ? <CircularProgress size={24} color="inherit" /> : "Search"}
                </Button>
              </Grid>
            </Grid>
          </form>
        </CardContent>
      </Card>

      {error && (
        <Paper sx={{ p: 2, mb: 3, bgcolor: "#fdeded" }}>
          <Typography color="error">{error}</Typography>
        </Paper>
      )}

      {diseaseData && (
        <Card elevation={3}>
          <CardContent>
            <Typography variant="h6" gutterBottom>
              Results for "{diseaseData.disease_term}"
            </Typography>
            
            <Box sx={{ mb: 2 }}>
              <Grid container spacing={2}>
                <Grid item xs={6} sm={4}>
                  <Typography variant="body2" color="text.secondary">Status</Typography>
                  <Chip 
                    label={diseaseData.status} 
                    color={diseaseData.status === "completed" ? "success" : "primary"} 
                    size="small"
                  />
                </Grid>
                <Grid item xs={6} sm={4}>
                  <Typography variant="body2" color="text.secondary">Max Results</Typography>
                  <Typography variant="body1">{diseaseData.max_results}</Typography>
                </Grid>
                <Grid item xs={6} sm={4}>
                  <Typography variant="body2" color="text.secondary">References Found</Typography>
                  <Typography variant="body1">{diseaseData.count}</Typography>
                </Grid>
              </Grid>
            </Box>
            
            <Divider sx={{ my: 2 }} />
            
            {/* <Typography variant="subtitle1" gutterBottom>
              Downloaded Files
            </Typography>
            
            {disease.length > 0 ? (
              <List sx={{ bgcolor: "#f5f5f5", borderRadius: 1 }}>
                {disease.map((file, index) => (
                  <ListItem key={index} divider={index < disease.length - 1}>
                    <ListItemText primary={file} />
                  </ListItem>
                ))}
              </List>
            ) : (
              <Typography variant="body2" color="text.secondary">
                No files downloaded yet.
              </Typography>
            )} */} 
            {/* Will be implemented in the next step */}
          </CardContent>
        </Card>
      )}
    </Container>
  );
}