import { JSX } from 'react';
import { Box, Typography, useTheme } from '@mui/material';
import { Column, Container, Row } from '@/components/core';
import HowTo, { HowToProps } from './buttons/HowTo';

export interface RouteHeaderProps {
  icon?: React.ElementType;
  title: string;
  description?: string;
  howTo?: HowToProps;
  sx?: boolean;
}

export default function RouteHeader({ icon: Icon, title, description, howTo, sx = false }: RouteHeaderProps): JSX.Element {
  const theme = useTheme();

  return (
    <Column id={title} sx={{ borderBottom: 1, borderColor: 'divider' }}>
      <Container>
        <Row sx={{ p: 0 }}>
          {Icon && (
            <Icon sx={{ fontSize: sx ? 36 : 48 , color: theme.palette.mode === 'dark' ? theme.palette.primary.dark : theme.palette.primary.main }}/>
          )}
          <Typography variant={sx ? 'h4' : 'h3' } component="h1" fontWeight={700}>
            {title}
          </Typography>
          {howTo && (
            <Box sx={{ display: 'flex', flexDirection: 'row', flex: 1, justifyContent: 'flex-end' }}>
              <HowTo
                media={howTo.media}
                title={howTo.title}
                description={howTo.description}
              />
            </Box>
          )}
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
