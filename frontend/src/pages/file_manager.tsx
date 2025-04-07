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
import { KeyboardReturn,
  Folder,InsertDriveFile, Description, TableChart, Storage,
  BlurOn, PictureAsPdf, Terminal, Coronavirus, Image,
  Movie, Audiotrack, Archive } from '@mui/icons-material';
import { useDropzone } from 'react-dropzone';

export default function FileManager() {
  const navigate = useNavigate();
  const [fileList, setFileList] = useState<{
    filename: string;
    type: string;
    size_kb: number | null;
    item_count: number | null;
  }[]>([]);
  const [currentPath, setCurrentPath] = useState<string>('/');
  const [inputPath, setInputPath] = useState<string>(currentPath);
  const [searchQuery, setSearchQuery] = useState<string>(''); // State for search query

  useEffect(() => {
    const fetchFiles = async (path: string) => {
      try {
        const files = await window.pywebview.api.file_manager.list_files(path);
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

  const handleSearchChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    setSearchQuery(event.target.value); // Update search query
  };

  const filteredFiles = fileList.filter((file) =>
    file.filename.toLowerCase().includes(searchQuery.toLowerCase())
  ); // Filter files by search query

  const getFileIcon = (fileType: string) => {
    switch (fileType) {
      case 'text file':
        return <Description sx={{ mr: 1, color: 'text.secondary' }} />; // Icon for .txt
      case 'CSV':
        return <TableChart sx={{ mr: 1, color: 'text.secondary' }} />; // Icon for .csv
      case 'fasta':
        return <Coronavirus sx={{ mr: 1, color: 'text.secondary' }} />; // Icon for .fasta or .fa
      case 'VCF':
        return <BlurOn sx={{ mr: 1, color: 'text.secondary' }} />; // Icon for .vcf
      case 'database':
        return <Storage sx={{ mr: 1, color: 'text.secondary' }} />; // Icon for .db or .sqlite
      case 'PDF':
        return <PictureAsPdf sx={{ mr: 1, color: 'error.main' }} />; // Icon for .pdf
      case 'executable':
        return <Terminal sx={{ mr: 1, color: 'success.main' }} />; // Icon for executables (.exe, .sh, etc.)
      case 'folder':
        return <Folder sx={{ mr: 1, color: 'primary.main' }} />; // Icon for folders
      case 'image':
        return <Image sx={{ mr: 1, color: 'info.main' }} />; // Icon for image files
      case 'video':
        return <Movie sx={{ mr: 1, color: 'info.main' }} />; // Icon for video files
      case 'audio':
        return <Audiotrack sx={{ mr: 1, color: 'info.main' }} />; // Icon for audio files
      case 'archive':
        return <Archive sx={{ mr: 1, color: 'warning.main' }} />; // Icon for archive files
      default:
        return <InsertDriveFile sx={{ mr: 1, color: 'text.secondary' }} />; // Default file icon
    }
  };

  const onDrop = async (acceptedFiles: File[]) => {
    try {
      for (const file of acceptedFiles) {
        const fileContent = await file.arrayBuffer();
        const fileName = file.name;

        await window.pywebview.api.file_manager.upload_file(
          currentPath,
          fileName,
          Array.from(new Uint8Array(fileContent))
        );
      }

      // Refresh the file list after upload
      const updatedFiles = await window.pywebview.api.file_manager.list_files(currentPath);
      setFileList(updatedFiles);
    } catch (error) {
      console.error('Error uploading files:', error);
    }
  };

  const { getRootProps, getInputProps, isDragActive } = useDropzone({
    onDrop,
    onDragEnter: () => console.warn('Drag entered'),
    onDragOver: () => console.warn('Dragging over'),
    onDragLeave: () => console.warn('Drag left'),
    multiple: false
  });

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
        <Box sx={{ mt: 2, mb: 3, display: 'flex', alignItems: 'center' }}>
          <TextField
            fullWidth
            value={inputPath}
            onChange={handlePathChange}
            onKeyDown={handlePathSubmit}
            variant="outlined"
            size="small"
            sx={{ flexGrow: 1 }}
          />
          <KeyboardReturn
            sx={{
              ml: 1,
              cursor: 'pointer',
              color: 'primary.main',
            }}
            onClick={() => setCurrentPath(inputPath)}
          />
        </Box>

{/* Drag-and-Drop Area */}
<Box
          {...getRootProps()}
          sx={{
            border: '2px dashed #ddd',
            borderRadius: 2,
            p: 3,
            textAlign: 'center',
            backgroundColor: isDragActive ? '#f0f0f0' : 'transparent',
            cursor: 'pointer',
          }}
        >
          <input {...getInputProps()} type="file" />
          {isDragActive ? (
            <Typography variant="body1" color="primary">
              Drop the files here...
            </Typography>
          ) : (
            <Typography variant="body1" color="textSecondary">
              Drag and drop files here, or click to select files
            </Typography>
          )}
        </Box>
        {/* Search Bar */}
        <Box sx={{ mb: 3 }}>
          <TextField
            fullWidth
            value={searchQuery}
            onChange={handleSearchChange}
            placeholder="Search by filename..."
            variant="outlined"
            size="small"
          />
        </Box>

                <Box sx={{ mt: 3 }}>
          {fileList.length > 0 ? (
            <>
              {/* Column Headers */}
              <Box
                sx={{
                  display: 'flex',
                  alignItems: 'center',
                  justifyContent: 'space-between',
                  p: 2,
                  fontWeight: 'bold',
                }}
              >
                <Typography variant="body1" sx={{ flex: 2 }}>
                  Name
                </Typography>
                <Typography variant="body1" sx={{ flex: 1 }}>
                  Type
                </Typography>
                <Typography variant="body1" sx={{ flex: 1 }}>
                  Size
                </Typography>
              </Box>

              {/* Parent Directory */}
              {currentPath !== '/' && (
                <Box
                  sx={{
                    display: 'flex',
                    alignItems: 'center',
                    justifyContent: 'space-between',
                    p: 2,
                    borderBottom: '1px solid #ddd',
                    cursor: 'pointer',
                  }}
                  onClick={handleNavigateUp}
                >
                  <Box sx={{ display: 'flex', alignItems: 'center', flex: 2 }}>
                    <Folder sx={{ mr: 1, color: 'primary.main' }} />
                    <Typography variant="body1" sx={{ fontWeight: 'bold' }}>
                      ..
                    </Typography>
                  </Box>
                  <Typography variant="body2" sx={{ flex: 1 }} color="textSecondary">
                    Parent Directory
                  </Typography>
                  <Typography variant="body2" sx={{ flex: 1 }} color="textSecondary">
                    Folder
                  </Typography>
                </Box>
              )}

              {/* Filtered Files */}
              {filteredFiles.map((file, index) => (
                <Box
                  key={index}
                  sx={{
                    display: 'flex',
                    alignItems: 'center',
                    justifyContent: 'space-between',
                    p: 2,
                    borderBottom: '1px solid #ddd',
                    cursor: file.type === 'folder' ? 'pointer' : 'default',
                  }}
                  onClick={() => file.type === 'folder' && handleFolderClick(file.filename)}
                >
                  <Box sx={{ display: 'flex', alignItems: 'center', flex: 2 }}>
                    {getFileIcon(file.type)}
                    <Typography variant="body1" sx={{ fontWeight: 'bold' }}>
                      {file.filename}
                    </Typography>
                  </Box>
                  <Typography variant="body2" sx={{ flex: 1 }} color="textSecondary">
                    {file.type}
                  </Typography>
                  <Typography variant="body2" sx={{ flex: 1 }} color="textSecondary">
                    {file.type === 'folder'
                      ? `${file.item_count} items`
                      : `${file.size_kb?.toFixed(2)} KB`}
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
