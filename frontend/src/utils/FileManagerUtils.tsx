import {
  Description as DescriptionIcon,
  TableChart as TableChartIcon,
  Coronavirus as CoronavirusIcon,
  BlurOn as BlurOnIcon,
  Storage as StorageIcon,
  PictureAsPdf as PictureAsPdfIcon,
  Terminal as TerminalIcon,
  Folder as FolderIcon,
  Image as ImageIcon,
  Movie as MovieIcon,
  Audiotrack as AudiotrackIcon,
  Archive as ArchiveIcon,
  InsertDriveFile as InsertDriveFileIcon,
} from '@mui/icons-material';

export const getFileIcon = (fileType: string) => {
  switch (fileType) {
    case 'text file':
      return <DescriptionIcon sx={{ mr: 1, color: 'text.secondary' }} />;
    case 'CSV':
      return <TableChartIcon sx={{ mr: 1, color: 'text.secondary' }} />;
    case 'fasta':
      return <CoronavirusIcon sx={{ mr: 1, color: 'text.secondary' }} />;
    case 'VCF':
      return <BlurOnIcon sx={{ mr: 1, color: 'text.secondary' }} />;
    case 'database':
      return <StorageIcon sx={{ mr: 1, color: 'text.secondary' }} />;
    case 'PDF':
      return <PictureAsPdfIcon sx={{ mr: 1, color: 'error.main' }} />;
    case 'executable':
      return <TerminalIcon sx={{ mr: 1, color: 'success.main' }} />;
    case 'folder':
      return <FolderIcon sx={{ mr: 1, color: 'primary.main' }} />;
    case 'image':
      return <ImageIcon sx={{ mr: 1, color: 'info.main' }} />;
    case 'video':
      return <MovieIcon sx={{ mr: 1, color: 'info.main' }} />;
    case 'audio':
      return <AudiotrackIcon sx={{ mr: 1, color: 'info.main' }} />;
    case 'archive':
      return <ArchiveIcon sx={{ mr: 1, color: 'warning.main' }} />;
    default:
      return <InsertDriveFileIcon sx={{ mr: 1, color: 'text.secondary' }} />;
  }
};
