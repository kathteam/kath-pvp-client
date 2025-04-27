import React from 'react';
import { Box, Typography } from '@mui/material';
import { useDropzone } from 'react-dropzone';

interface DragAndDropProps {
  onDrop: (acceptedFiles: File[]) => Promise<void>;
}

const DragAndDrop: React.FC<DragAndDropProps> = ({ onDrop }) => {
  const { getRootProps, getInputProps, isDragActive } = useDropzone({
    onDrop,
    onDragEnter: () => console.warn('Drag entered'),
    onDragOver: () => console.warn('Dragging over'),
    onDragLeave: () => console.warn('Drag left'),
    multiple: false,
  });

  return (
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
  );
};

export default DragAndDrop;
