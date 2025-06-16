import { Button as MuiButton, styled } from '@mui/material';

const Button = styled(MuiButton)(({ theme }) => ({
  backgroundColor: theme.palette.mode === 'dark' ? theme.palette.primary.dark : theme.palette.primary.main,
}));

export default Button;
