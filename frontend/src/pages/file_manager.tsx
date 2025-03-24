import { useEffect, useState } from 'react';
import { useNavigate } from 'react-router-dom';
import {
  Typography,
  Box,
  Button,
  Container,
  Paper,
  Grid
} from '@mui/material';

export default function FileManager(){
  const navigate = useNavigate();

  const [fileList, setFileList] = useState<{ name: string; size: number }[]>([]);

  useEffect(() => {
    const fetchFiles = async () => {
      try {
        // Fetch the files from pywebview API
        const files = await window.pywebview.api.list_files();
        const formattedFiles = files.map((file: string) => ({
          name: file,
          size: 0, // Replace with actual size if available
        }));
        setFileList(formattedFiles); // Store files in state to display them
      } catch (error) {
        console.error('Error fetching files:', error);
      }
    };

    fetchFiles(); // Call the function when the component mounts
  }, []);

  return (
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
          File Manager
        </Typography>
        <Typography variant="body1">
          Here is a list of the files retrieved from the server:
        </Typography>

        <Grid container spacing={2} sx={{ mt: 3 }}>
          {fileList.length > 0 ? (
            fileList.map((file, index) => (
              <Grid item xs={12} sm={6} md={4} key={index}>
                <Box
                  sx={{
                    p: 2,
                    border: '1px solid #ddd',
                    borderRadius: 2,
                    boxShadow: 2,
                    textAlign: 'center',
                  }}
                >
                  <Typography variant="h6" component="h2" sx={{ mb: 2 }}>
                    {file.name}
                  </Typography>
                  <Typography variant="body2" color="textSecondary">
                    {file.size} KB
                  </Typography>
                </Box>
              </Grid>
            ))
          ) : (
            <Typography variant="body1">No files available.</Typography>
          )}
        </Grid>
        <Box sx={{ mt: 3 }}>
          <Button
            variant="contained"
            color="primary"
            onClick={() => navigate('/dashboard')}
          >
            Back to Dashboard
          </Button>
        </Box>
      </Paper>
    </Container>
  );
}
