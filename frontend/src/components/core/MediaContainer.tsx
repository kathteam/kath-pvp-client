import { Box as MuiBox, styled } from '@mui/material';

const MediaContainer = styled(MuiBox)(({ theme }) => ({
  width: '100%',
  maxWidth: '1000px',
  backgroundColor: theme.palette.background.paper,
}));

export default MediaContainer;
