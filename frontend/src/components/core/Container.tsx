import { Box as MuiBox, styled } from '@mui/material';

const Container = styled(MuiBox)(({ theme }) => ({
  width: '100%',
  maxWidth: '2000px',
  backgroundColor: theme.palette.background.paper,
}));

export default Container;
