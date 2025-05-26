import { Box as MuiBox, styled } from '@mui/material';

const Row = styled(MuiBox)(({ theme }) => ({
  display: 'flex',
  flexDirection: 'row',
  alignItems: 'center',
  padding: theme.spacing(4),
  gap: theme.spacing(2),
  backgroundColor: theme.palette.background.paper,
}));

export default Row;
