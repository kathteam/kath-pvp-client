import { Box as MuiBox, styled } from '@mui/material';

const Column = styled(MuiBox)(({ theme }) => ({
  display: 'flex',
  flexDirection: 'column',
  alignItems: 'center',
  flex: 1,
  padding: theme.spacing(4),
  gap: theme.spacing(2),
  backgroundColor: theme.palette.background.paper,
}));

export default Column;
