import { JSX } from 'react';
import { Typography, useTheme } from '@mui/material';
import { Column, Container, Row } from '@/components/core';

export interface RouteHeaderProps {
  icon: React.ElementType;
  title: string;
  description?: string;
}

export default function RouteHeader({ icon: Icon, title, description }: RouteHeaderProps): JSX.Element {
  const theme = useTheme();

  return (
    <Column id={title}>
      <Container>
        <Row sx={{ p: 0 }}>
          <Icon sx={{ fontSize: 48, color: theme.palette.mode === 'dark' ? theme.palette.primary.dark : theme.palette.primary.main }}/>
          <Typography variant="h3" component="h1" fontWeight={700}>
            {title}
          </Typography>
        </Row>
        <Row sx={{ px: 0, pb: 0 }}>
          <Typography variant="body1" color='text.secondary'>
            {description}
          </Typography>
        </Row>
      </Container>
    </Column>
  );
}
