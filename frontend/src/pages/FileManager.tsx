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
  BlurOn, PictureAsPdf, Terminal, Coronavirus, Image, Movie, Audiotrack, Archive } from '@mui/icons-material';
import { Menu, MenuItem, Dialog, DialogActions, DialogTitle, DialogContent, DialogContentText } from '@mui/material';
import DragAndDrop from '../components/cards/DragAndDrop';
import FileList from '../components/cards/FileListCard';

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
  const [searchQuery] = useState<string>(''); // State for search query
  const [typeFilter] = useState<string>(''); // Filter for file type
  const [sizeFilter] = useState<string>(''); // Filter for file size
  const [sortConfig, setSortConfig] = useState<{ column: string; direction: 'asc' | 'desc' | null }>({
    column: '',
    direction: null,
  });

  const [openDbDialog, setOpenDbDialog] = useState(false);
  const [dbFilename, setDbFilename] = useState<string>('personalized_gene_database.db');

  const handleCreateDatabase = async () => {
    try {
      const dbPath = `${currentPath}/${dbFilename}`;
      await window.pywebview.api.file_controller.create_vcf_database(dbPath);
      console.log(`Database created at: ${dbPath}`);
      alert('Personalized gene database created successfully!');

      // Refresh the file list after creating the database
      const updatedFiles = await window.pywebview.api.file_controller.list_files(currentPath);
      setFileList(updatedFiles);
    } catch (error) {
      console.error('Error creating personalized gene database:', error);
      alert('Failed to create personalized gene database.');
    }
    setOpenDbDialog(false);
  };

  useEffect(() => {
    const fetchFiles = async (path: string) => {
      try {
        const files = await window.pywebview.api.file_controller.list_files(path);
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
    setInputPath(parentPath);
  };

  const handleFolderClick = (folderName: string) => {
    const newPath = `${currentPath}/${folderName}`.replace('//', '/');
    setCurrentPath(newPath);
    setInputPath(newPath);
  };

  const handlePathChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    setInputPath(event.target.value);
  };

  const handlePathSubmit = (event: React.KeyboardEvent<HTMLInputElement>) => {
    if (event.key === 'Enter') {
      setCurrentPath(inputPath);
    }
  };

  const handleSort = (column: string) => {
    setSortConfig((prev) => {
      if (prev.column === column) {
        // Toggle sort direction
        const newDirection = prev.direction === 'asc' ? 'desc' : prev.direction === 'desc' ? null : 'asc';
        return { column, direction: newDirection };
      }
      return { column, direction: 'asc' }; // Default to ascending
    });
  };

  const filteredFiles = fileList.filter((file) => {
    const matchesName = file.filename.toLowerCase().includes(searchQuery.toLowerCase());
    const matchesType = typeFilter ? file.type.toLowerCase().includes(typeFilter.toLowerCase()) : true;
    
    const sizeColumn = file.type === 'folder'
      ? `${file.item_count} items`
      : `${file.size_kb?.toFixed(2)} KB`;
    const matchesSize = sizeFilter
      ? sizeColumn.toLowerCase().includes(sizeFilter.toLowerCase())
      : true;

    return matchesName && matchesType && matchesSize;
  });

  const sortedFiles = [...filteredFiles].sort((a, b) => {
    if (!sortConfig.direction) return 0; // No sorting
    const isAsc = sortConfig.direction === 'asc';
    if (sortConfig.column === 'Name') {
      return isAsc ? a.filename.localeCompare(b.filename) : b.filename.localeCompare(a.filename);
    }
    if (sortConfig.column === 'Type') {
      return isAsc ? a.type.localeCompare(b.type) : b.type.localeCompare(a.type);
    }
    if (sortConfig.column === 'Size') {
      const sizeA = a.type === 'folder'
        ? `${a.item_count} items`
        : `${a.size_kb?.toFixed(2)} KB`;
      const sizeB = b.type === 'folder'
        ? `${b.item_count} items`
        : `${b.size_kb?.toFixed(2)} KB`;
      return isAsc ? sizeA.localeCompare(sizeB) : sizeB.localeCompare(sizeA);
    }
    return 0;
  });

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

        await window.pywebview.api.file_controller.upload_file(
          currentPath,
          fileName,
          Array.from(new Uint8Array(fileContent))
        );
      }

      // Refresh the file list after upload
      const updatedFiles = await window.pywebview.api.file_controller.list_files(currentPath);
      setFileList(updatedFiles);
    } catch (error) {
      console.error('Error uploading files:', error);
    }
  };

  const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);
  const open = Boolean(anchorEl);

  const handleClose = () => {
    setAnchorEl(null);
    setSelectedFileForMenu(null); // Reset the selected file
  };

  const [selectedFileForMenu, setSelectedFileForMenu] = useState<string | null>(null);

  const handleOptionsClick = (event: React.MouseEvent<SVGSVGElement>, filename: string) => {
    setAnchorEl(event.currentTarget as unknown as HTMLElement);
    setSelectedFileForMenu(filename); // Set the file associated with the menu
  };

  const [openRenameDialog, setOpenRenameDialog] = useState(false);
  const [openDeleteDialog, setOpenDeleteDialog] = useState(false);
  const [selectedFile, setSelectedFile] = useState<string | null>(null);  

  const handleRenameFile = (oldName: string) => async (event: React.MouseEvent<HTMLElement>) => {
    event.stopPropagation(); // Prevent the menu from closing before the action is complete
    setSelectedFile(oldName);
    setOpenRenameDialog(true);
  };

  const handleRenameConfirm = async () => {
    if (selectedFile) {
      const newName = (document.getElementById('new-name') as HTMLInputElement).value;
      await window.pywebview.api.file_controller.rename_file(currentPath, selectedFile, newName);
      const updatedFiles = await window.pywebview.api.file_controller.list_files(currentPath);
      setFileList(updatedFiles);
    }
    setOpenRenameDialog(false);
  };

  const handleDeleteFile = (filename: string) => async (event: React.MouseEvent<HTMLElement>) => {
    event.stopPropagation();
    setSelectedFile(filename);
    setOpenDeleteDialog(true);
  };

  const handleDeleteConfirm = async () => {
    if (selectedFile) {
      try {
        await window.pywebview.api.file_controller.delete_file(currentPath, selectedFile);
        const updatedFiles = await window.pywebview.api.file_controller.list_files(currentPath);
        setFileList(updatedFiles);
      } catch (error) {
        console.error('Error deleting file:', error);
      }
    }
    setOpenDeleteDialog(false);
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
        <DragAndDrop onDrop={onDrop} />

        {/* File List Component */}
        <Box sx={{ mt: 3 }}>
          <FileList
            files={sortedFiles}
            currentPath={currentPath}
            sortConfig={sortConfig}
            onSort={handleSort}
            onNavigateUp={handleNavigateUp}
            onFolderClick={handleFolderClick}
            onOptionsClick={handleOptionsClick}
            getFileIcon={getFileIcon}
          />
        </Box>

        <Menu
          id="options-menu"
          anchorEl={anchorEl}
          open={open}
          onClose={handleClose}
          MenuListProps={{
            'aria-labelledby': 'options-button',
          }}
        >
          {selectedFileForMenu && (
            <>
              <MenuItem onClick={(e) => handleRenameFile(selectedFileForMenu)(e)}>Rename</MenuItem>
              <MenuItem onClick={(e) => handleDeleteFile(selectedFileForMenu)(e)}>Delete</MenuItem>
            </>
          )}
        </Menu>

        <Box sx={{ mt: 3 }}>
          <Button
            variant="contained"
            color="primary"
            onClick={() => navigate('/dashboard')}
          >
            Back to Dashboard
          </Button>
        </Box>
        <Box sx={{ mt: 3 }}>
          <Button
            variant="contained"
            color="primary"
            onClick={async () => {
              try {
                const kathDirectory = await window.pywebview.api.file_controller.get_kath_directory();
                setCurrentPath(kathDirectory);
                setInputPath(kathDirectory);
                const files = await window.pywebview.api.file_controller.list_files(kathDirectory);
                setFileList(files);
              } catch (error) {
                console.error('Error changing directory to kath:', error);
              }
            }}
          >
            Change directory to kath
          </Button>
        </Box>
        <Box sx={{ mt: 3 }}>
          <Button
            variant="contained"
            color="primary"
            onClick={() => setOpenDbDialog(true)}
          >
            Create personalized gene database
          </Button>
        </Box>

        <Dialog open={openRenameDialog} onClose={() => setOpenRenameDialog(false)}>
          <DialogTitle>Rename File</DialogTitle>
          <DialogContent>
            <DialogContentText>
              Enter a new name for the file:
            </DialogContentText>
            <TextField
              autoFocus
              margin="dense"
              id="new-name"
              label="New Name"
              type="text"
              fullWidth
              variant="standard"
            />
          </DialogContent>
          <DialogActions>
            <Button onClick={() => setOpenRenameDialog(false)}>Cancel</Button>
            <Button onClick={handleRenameConfirm}>Rename</Button>
          </DialogActions>
        </Dialog>
        <Dialog open={openDeleteDialog} onClose={() => setOpenDeleteDialog(false)}>
          <DialogTitle>Delete File</DialogTitle>
          <DialogContent>
            <DialogContentText>
              Are you sure you want to delete the file?
            </DialogContentText>
          </DialogContent>
          <DialogActions>
            <Button onClick={() => setOpenDeleteDialog(false)}>Cancel</Button>
            <Button onClick={handleDeleteConfirm}>Delete</Button>
          </DialogActions>
        </Dialog>
        <Dialog open={openDbDialog} onClose={() => setOpenDbDialog(false)}>
          <DialogTitle>Create Database</DialogTitle>
          <DialogContent>
            <DialogContentText>
              Enter a name for the database file:
            </DialogContentText>
            <TextField
              autoFocus
              margin="dense"
              id="db-filename"
              label="Database Filename"
              type="text"
              fullWidth
              variant="standard"
              value={dbFilename}
              onChange={(e) => setDbFilename(e.target.value)}
            />
          </DialogContent>
          <DialogActions>
            <Button onClick={() => setOpenDbDialog(false)}>Cancel</Button>
            <Button onClick={handleCreateDatabase}>Create</Button>
          </DialogActions>
        </Dialog>
      </Paper>
    </Container>
  );
}
