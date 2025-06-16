import { Typography as MuiTypography, styled } from '@mui/material';

const Header = styled(MuiTypography)(({ theme }) => ({
  color: theme.palette.mode === 'dark' ? theme.palette.primary.dark : theme.palette.primary.main,
}));

export default Header;
