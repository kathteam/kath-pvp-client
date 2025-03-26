import { useEffect, useState } from 'react';
import { useNavigate } from 'react-router-dom';
import {
  Typography,
  Box,
  Button,
  Container,
  Paper,
  TextField
} from '@mui/material';

export default function FileManager() {
  const navigate = useNavigate();
  const [fileList, setFileList] = useState<{
    filename: string;
    type: string;
    size_kb: number | null;
    item_count: number | null;
  }[]>([]);
  const [currentPath, setCurrentPath] = useState<string>('/'); // Track the current path
  const [inputPath, setInputPath] = useState<string>(currentPath); // Track the user input for the path

  useEffect(() => {
    const fetchFiles = async (path: string) => {
      try {
        const files = await window.pywebview.api.list_files(path);
        setFileList(files);
      } catch (error) {
        console.error('Error fetching files:', error);
      }
    };

    fetchFiles(currentPath);
  }, [currentPath]);

  const handleNavigateUp = () => {
    const parentPath = currentPath.substring(0, currentPath.lastIndexOf('/')) || '/';
    setCurrentPath(parentPath);
    setInputPath(parentPath); // Update the input field as well
  };

  const handleFolderClick = (folderName: string) => {
    const newPath = `${currentPath}/${folderName}`.replace('//', '/');
    setCurrentPath(newPath);
    setInputPath(newPath); // Update the input field as well
  };

  const handlePathChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    setInputPath(event.target.value); // Update the input field value
  };

  // Update the current path when <Enter> is pressed
  const handlePathSubmit = (event: React.KeyboardEvent<HTMLInputElement>) => {
    if (event.key === 'Enter') {
      setCurrentPath(inputPath); 
    }
  };

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

        {/* Path Bar */}
        <Box sx={{ mt: 2, mb: 3 }}>
          <TextField
            fullWidth
            value={inputPath}
            onChange={handlePathChange}
            onKeyDown={handlePathSubmit}
            variant="outlined"
            size="small"
          />
        </Box>

        <Box sx={{ mt: 3 }}>
          {fileList.length > 0 ? (
            <>
              {/* .. folder for navigation */}
              {currentPath !== '/' && (
                <Box
                  sx={{
                    display: 'flex',
                    justifyContent: 'space-between',
                    p: 2,
                    borderBottom: '1px solid #ddd',
                    cursor: 'pointer',
                  }}
                  onClick={handleNavigateUp}
                >
                  <Typography variant="body1" sx={{ fontWeight: 'bold' }}>
                    ..
                  </Typography>
                  <Typography variant="body2" color="textSecondary">
                    Parent Directory
                  </Typography>
                </Box>
              )}
              {fileList.map((file, index) => (
                <Box
                  key={index}
                  sx={{
                    display: 'flex',
                    justifyContent: 'space-between',
                    p: 2,
                    borderBottom: '1px solid #ddd',
                    cursor: file.type === 'folder' ? 'pointer' : 'default',
                  }}
                  onClick={() => file.type === 'folder' && handleFolderClick(file.filename)}
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
              ))}
            </>
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
