import {
  Description,
  TableChart,
  Coronavirus,
  BlurOn,
  Storage,
  PictureAsPdf,
  Terminal,
  Folder,
  Image,
  Movie,
  Audiotrack,
  Archive,
  InsertDriveFile,
} from '@mui/icons-material';

export const getFileIcon = (fileType: string) => {
  switch (fileType) {
    case 'text file':
      return <Description sx={{ mr: 1, color: 'text.secondary' }} />;
    case 'CSV':
      return <TableChart sx={{ mr: 1, color: 'text.secondary' }} />;
    case 'fasta':
      return <Coronavirus sx={{ mr: 1, color: 'text.secondary' }} />;
    case 'VCF':
      return <BlurOn sx={{ mr: 1, color: 'text.secondary' }} />;
    case 'database':
      return <Storage sx={{ mr: 1, color: 'text.secondary' }} />;
    case 'PDF':
      return <PictureAsPdf sx={{ mr: 1, color: 'error.main' }} />;
    case 'executable':
      return <Terminal sx={{ mr: 1, color: 'success.main' }} />;
    case 'folder':
      return <Folder sx={{ mr: 1, color: 'primary.main' }} />;
    case 'image':
      return <Image sx={{ mr: 1, color: 'info.main' }} />;
    case 'video':
      return <Movie sx={{ mr: 1, color: 'info.main' }} />;
    case 'audio':
      return <Audiotrack sx={{ mr: 1, color: 'info.main' }} />;
    case 'archive':
      return <Archive sx={{ mr: 1, color: 'warning.main' }} />;
    default:
      return <InsertDriveFile sx={{ mr: 1, color: 'text.secondary' }} />;
  }
};