import { JSX } from "react";
import { Modal, Box, Typography, Chip, Button, Divider, Paper } from "@mui/material";
import HealthAndSafetyIcon from '@mui/icons-material/HealthAndSafety';
import CloseIcon from '@mui/icons-material/Close';

export default function DiseaseModal(props: {
  clinicalSignificance: string;
  disease: string;
  onClose?: () => void;
}): JSX.Element {
  
  const handleClose = () => {
    if (props.onClose) {
      props.onClose(); 
    }
  };

  const getSeverityColor = (significance: string) => {
    const lowerSignificance = significance.toLowerCase();
    if (lowerSignificance.includes('pathogenic') || lowerSignificance.includes('high')) return 'error';
    if (lowerSignificance.includes('likely') || lowerSignificance.includes('moderate')) return 'warning';
    if (lowerSignificance.includes('benign') || lowerSignificance.includes('low')) return 'success';
    return 'default';
  };

  return (
    <Modal
      open={true}
      onClose={handleClose}
      aria-labelledby="disease-information"
      aria-describedby="genetic-disease-details"
    >
      <Paper sx={{
        position: 'absolute',
        top: '50%',
        left: '50%',
        transform: 'translate(-50%, -50%)',
        width: 450,
        bgcolor: 'background.paper',
        borderRadius: 2,
        boxShadow: 24,
        p: 0,
        overflow: 'hidden'
      }}>
        <Box sx={{ 
          bgcolor: 'primary.main', 
          color: 'primary.contrastText', 
          py: 2, 
          px: 3, 
          display: 'flex', 
          alignItems: 'center',
          justifyContent: 'space-between'
        }}>
          <Box sx={{ display: 'flex', alignItems: 'center' }}>
            <HealthAndSafetyIcon sx={{ mr: 1 }} />
            <Typography variant="h6" component="h2">
              Genetic Disease Information
            </Typography>
          </Box>
          <Button 
            color="inherit" 
            onClick={handleClose} 
            sx={{ minWidth: 'auto', p: 0.5 }}
          >
            <CloseIcon />
          </Button>
        </Box>

        <Box sx={{ p: 3 }}>
          <Box sx={{ mb: 3 }}>
            <Typography variant="subtitle2" color="text.secondary" gutterBottom>
              CLINICAL SIGNIFICANCE
            </Typography>
            <Chip 
              label={props.clinicalSignificance} 
              color={getSeverityColor(props.clinicalSignificance) as any}
              sx={{ fontWeight: 'bold', fontSize: '0.9rem' }}
            />
          </Box>

          <Divider sx={{ my: 2 }} />

          <Box>
            <Typography variant="subtitle2" color="text.secondary" gutterBottom>
              ASSOCIATED DISEASE
            </Typography>
            <Typography variant="body1" fontWeight="500">
              {props.disease}
            </Typography>
          </Box>

          <Box sx={{ mt: 4, display: 'flex', justifyContent: 'flex-end' }}>
            <Button 
              variant="contained" 
              onClick={handleClose}
            >
              Close
            </Button>
          </Box>
        </Box>
      </Paper>
    </Modal>
  );
}