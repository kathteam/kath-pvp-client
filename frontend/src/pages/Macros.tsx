import { JSX } from 'react';
import { useNavigate } from 'react-router-dom';
import {
  Typography,
  Button
} from '@mui/material';
import { Column, Row } from '@/components/core';

export default function Macros(): JSX.Element {
  const navigate = useNavigate();

  return (
    <Column>
      <Typography variant="h4" component="h1">
        Macros Page
      </Typography>
      <Typography variant="body1">
        This is the macros page of our application.
      </Typography>
      <Row>
        <Button
          variant="contained"
          color="primary"
          onClick={() => navigate('/dashboard')}
        >
          Back to Dashboard
        </Button>
      </Row>
    </Column>
  );
}
