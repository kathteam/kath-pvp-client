import React from 'react';
import { Box, Typography } from '@mui/material';
import { MoreVert, ArrowDropUp, ArrowDropDown, Folder } from '@mui/icons-material';

interface File {
  filename: string;
  type: string;
  size_kb: number | null;
  item_count: number | null;
}

interface FileListProps {
  files: File[];
  currentPath: string;
  sortConfig: { column: string; direction: 'asc' | 'desc' | null };
  onSort: (column: string) => void;
  onNavigateUp: () => void;
  onFolderClick: (folderName: string) => void;
  onOptionsClick: (event: React.MouseEvent<SVGSVGElement>, filename: string) => void;
  getFileIcon: (fileType: string) => React.ReactNode;
}

const FileList: React.FC<FileListProps> = ({
  files,
  currentPath,
  sortConfig,
  onSort,
  onNavigateUp,
  onFolderClick,
  onOptionsClick,
  getFileIcon,
}) => {
  return (
    <Box sx={{ mt: 3 }}>
      {files.length > 0 ? (
        <>
          {/* Column Headers */}
          <Box
            sx={{
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'space-between',
              p: 2,
              fontWeight: 'bold',
              borderTop: '1px solid #ddd',
            }}
          >
            <Typography
              variant="body1"
              sx={{ flex: 2, cursor: 'pointer', display: 'flex', alignItems: 'center' }}
              onClick={() => onSort('Name')}
            >
              Name
              {sortConfig.column === 'Name' && sortConfig.direction && (
                sortConfig.direction === 'asc' ? <ArrowDropUp fontSize="small" /> : <ArrowDropDown fontSize="small" />
              )}
            </Typography>
            <Typography
              variant="body1"
              sx={{ flex: 1, cursor: 'pointer', display: 'flex', alignItems: 'center' }}
              onClick={() => onSort('Type')}
            >
              Type
              {sortConfig.column === 'Type' && sortConfig.direction && (
                sortConfig.direction === 'asc' ? <ArrowDropUp fontSize="small" /> : <ArrowDropDown fontSize="small" />
              )}
            </Typography>
            <Typography
              variant="body1"
              sx={{ flex: 1, cursor: 'pointer', display: 'flex', alignItems: 'center' }}
              onClick={() => onSort('Size')}
            >
              Size
              {sortConfig.column === 'Size' && sortConfig.direction && (
                sortConfig.direction === 'asc' ? <ArrowDropUp fontSize="small" /> : <ArrowDropDown fontSize="small" />
              )}
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
              onClick={onNavigateUp}
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

          {/* Files */}
          {files.map((file, index) => (
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
              onClick={() => file.type === 'folder' && onFolderClick(file.filename)}
            >
              <Box sx={{ display: 'flex', alignItems: 'center', flex: 2, whiteSpace: 'nowrap', overflow: 'hidden', textOverflow: 'ellipsis' }}>
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
              <MoreVert
                onClick={(e) => onOptionsClick(e, file.filename)}
                sx={{ cursor: 'pointer' }}
              />
            </Box>
          ))}
        </>
      ) : (
        <Typography variant="body1">No files available.</Typography>
      )}
    </Box>
  );
};

export default FileList;
