import { useEffect, useState } from 'react';
import { useNavigate } from 'react-router-dom';
import {
  Typography,
  Box,
  Button,
  Container,
  Paper,
  TextField,
  Menu,
  MenuItem,
} from '@mui/material';

import { KeyboardReturn } from '@mui/icons-material';
import DragAndDrop from '../components/cards/DragAndDrop';
import FileList from '../components/cards/FileListCard';
import { getFileIcon } from '@/utils/FileManagerUtils';
import { RenameDialog, DeleteDialog, CreateDatabaseDialog } from '../components/Dialogs';
import { Dialog, DialogTitle, DialogContent, CircularProgress } from '@mui/material';
import { DiseaseModal } from '@/components/modals';
import HowTo from '@/components/buttons/HowTo';

const howTo = 
{
  media: 'ExtractButton.gif',
  title: 'Extract genetic diseases from found variations',
  description: 'Select the .fasta file with three dots (as shown in the video), press \'Analyze\', you will receive a list of genetic diseases that are associated with the variations found. WARNING: File may not have any variations or diseases.'
};

export default function FileManager() {
  const navigate = useNavigate();
  const [fileList, setFileList] = useState<{
    filename: string;
    type: string;
    size_kb: number | null;
    item_count: number | null;
  }[]>([]);
  const [currentPath, setCurrentPath] = useState<string>('');
  const [inputPath, setInputPath] = useState<string>('');
  const [kathRootDirectory, setKathRootDirectory] = useState<string>('');
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
    const initializeKathDirectory = async () => {
      try {
        const kathDirectory = await window.pywebview.api.file_controller.get_kath_directory();
        setCurrentPath(kathDirectory);
        setInputPath(kathDirectory);
        setKathRootDirectory(kathDirectory);
      } catch (error) {
        console.error('Error getting kath directory:', error);
        setCurrentPath('/');
        setInputPath('/');
        setKathRootDirectory('/');
      }
    };
    
    initializeKathDirectory();
  }, []);

  useEffect(() => {
    const fetchFiles = async (path: string) => {
      try {
        if (path) {
          const files = await window.pywebview.api.file_controller.list_files(path);
          setFileList(files);
        }
      } catch (error) {
        console.error('Error fetching files:', error);
      }
    };

    fetchFiles(currentPath);
  }, [currentPath]);

  const handleNavigateUp = () => {
    const parentPath = currentPath.substring(0, currentPath.lastIndexOf('/')) || '/';
    
    if (parentPath.length >= kathRootDirectory.length && 
        parentPath.startsWith(kathRootDirectory)) {
      setCurrentPath(parentPath);
      setInputPath(parentPath);
    } else {
      // If trying to navigate above Kath directory reset to Kath directory
      setCurrentPath(kathRootDirectory);
      setInputPath(kathRootDirectory);
    }
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
      if (inputPath.startsWith(kathRootDirectory)) {
        setCurrentPath(inputPath);
      } else {
        // Reset to Kath directory if trying to navigate elsewhere
        setCurrentPath(kathRootDirectory);
        setInputPath(kathRootDirectory);
      }
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

  const [analysisDialogOpen, setAnalysisDialogOpen] = useState(false);
  const [validity, setValidity] = useState(false);
  const [analysisResult, setAnalysisResult] = useState<{
      status: string;
      result_file: string;
    }>();
  const [isAnalyzing, setIsAnalyzing] = useState(false);

  const handleAnalyzeFasta = async (filename: string) => {
    const filePath = `${currentPath}/${filename}`;
    setIsAnalyzing(true);
    setAnalysisDialogOpen(true);
  
    try {
      const result = await window.pywebview.api.blast_service.disease_extraction(filePath);
      setAnalysisResult(result);
      setValidity(true);
    } catch (error) {
      console.error('FASTA analysis failed:', error);
      setAnalysisResult({
        status: 'error',
        result_file: `Failed to extract disease information: ${error || 'Unknown error'}`,
      });
      setValidity(false);
    } finally {
      setIsAnalyzing(false);
    }
  };

  // Genetical disease data management
  const [showModal, setShowModal] = useState<boolean>(false);
  const [geneticDiseaseData, setGeneticDiseaseData] = useState<{
      clinicalSignificance: string;
      disease: string;
    }[]>([]);
  
  const handleDisplayGeneticDisease = async () => {
    if (!analysisResult?.result_file) {
      alert('Please extract disease information first.');
      return;
    }
  
    try {
      const diseaseData =
          await window.pywebview.api.disease_service.get_disease_data(
            analysisResult.result_file
          );
  
      if (!diseaseData) {
        // TODO UNCOMMENT
        // alert('No disease data found.');
        setGeneticDiseaseData([
          {
            clinicalSignificance: 'Pathogenic',
            disease: 'VERY EXAMPLE DISEASE',
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

  return (
    <Container
      maxWidth="md"
      sx={{
        display: 'flex',
        justifyContent: 'center',
        py: 4
      }}
    >
      <Paper 
        elevation={4}
        sx={{
          p: 5,
          width: '100%',
          borderRadius: 4,
          boxShadow: 6,
          textAlign: 'center',
          position: 'relative',
        }}>
        <Box sx={{ position: 'absolute', top: 5, right: 5 }}>
          <HowTo media={howTo.media} title={howTo.title} description={howTo.description} />
        </Box>
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
              {selectedFileForMenu?.endsWith('.fasta') && (
                <MenuItem onClick={() => {
                  handleAnalyzeFasta(selectedFileForMenu);
                  handleClose();
                }}>
                  Analyze
                </MenuItem>
              )}
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
            onClick={() => setOpenDbDialog(true)}
          >
            Create personalized gene database
          </Button>
        </Box>

        <RenameDialog
          open={openRenameDialog}
          onClose={() => setOpenRenameDialog(false)}
          onConfirm={handleRenameConfirm}
        />
        <DeleteDialog
          open={openDeleteDialog}
          onClose={() => setOpenDeleteDialog(false)}
          onConfirm={handleDeleteConfirm}
        />
        <CreateDatabaseDialog
          open={openDbDialog}
          onClose={() => setOpenDbDialog(false)}
          dbFilename={dbFilename}
          setDbFilename={setDbFilename}
          onConfirm={handleCreateDatabase}
        />
        <Dialog open={analysisDialogOpen} onClose={() => setAnalysisDialogOpen(false)} maxWidth="sm" fullWidth>
          <DialogTitle>FASTA Analysis Result</DialogTitle>
          <DialogContent dividers>
            {isAnalyzing ? (
              <Box sx={{ display: 'flex', justifyContent: 'center', my: 4 }}>
                <CircularProgress />
              </Box>
            ) : analysisResult ? (
              <>
                <Typography variant="h6">Status: {analysisResult.status}</Typography>
                <Typography variant="body1">Result File: {analysisResult.result_file}</Typography>
              </>
            ) : (
              <Typography variant="body1">No result.</Typography>
            )}
            {validity && (
              <Button
                variant="contained"
                color="primary"
                onClick={handleDisplayGeneticDisease}>
                            Display Genetic Disease Information
              </Button>)}
          </DialogContent>
        </Dialog>
      </Paper>
      {showModal && (
        <DiseaseModal
          diseases={geneticDiseaseData}
          onClose={() => setShowModal(false)}
        />)}
    </Container>
  );
}
