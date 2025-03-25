import { useEffect, useState } from 'react';
import { useNavigate } from 'react-router-dom';
import {
  Typography,
  Box,
  Button,
  Container,
  Paper
} from '@mui/material';

export default function FileManager() {
  const navigate = useNavigate();
  const [fileList, setFileList] = useState<{
    filename: string;
    type: string;
    size_kb: number | null;
    item_count: number | null;
  }[]>([]);

  useEffect(() => {
    const fetchFiles = async () => {
      try {
        const files = await window.pywebview.api.list_files();
        setFileList(files);
      } catch (error) {
        console.error('Error fetching files:', error);
      }
    };

    fetchFiles();
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
      <Paper elevation={0} sx={{ p: 3, width: '100%' }}>
        <Typography variant="h4" component="h1" gutterBottom>
          File Manager
        </Typography>

        <Box sx={{ mt: 3 }}>
          {fileList.length > 0 ? (
            fileList.map((file, index) => (
              <Box
                key={index}
                sx={{
                  display: 'flex',
                  justifyContent: 'space-between',
                  p: 2,
                  borderBottom: '1px solid #ddd',
                }}
              >
                <Typography variant="body1" sx={{ fontWeight: 'bold' }}>
                  {file.filename}
                </Typography>
                <Typography variant="body2" color="textSecondary">
                  {file.type === 'folder'
                    ? `${file.item_count} items`
                    : `${file.size_kb?.toFixed(2)} KB`}
                </Typography>
                <Typography variant="body2" color="textSecondary">
                  {file.type}
                </Typography>
              </Box>
            ))
          ) : (
            <Typography variant="body1">No files available.</Typography>
          )}
        </Box>

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
