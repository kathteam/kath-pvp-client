import { JSX } from 'react';
import { useNavigate } from 'react-router-dom';
import {
  Typography,
  Box,
  Button
} from '@mui/material';
import DashboardIcon from '@mui/icons-material/Dashboard';
import PlagiarismIcon from '@mui/icons-material/Plagiarism';
import InventoryIcon from '@mui/icons-material/Inventory';
import AutoAwesomeIcon from '@mui/icons-material/AutoAwesome';
import MenuBookIcon from '@mui/icons-material/MenuBook';
import { Column, Row } from '@/components/core';

export default function Dashboard(): JSX.Element {
  const navigate = useNavigate();

  return (
    <Box>
      <Column>
        <DashboardIcon color="primary" sx={{ fontSize: 48 }}/>
        <Typography variant="h3" component="h1" fontWeight={700}>
          Welcome to the Dashboard
        </Typography>
        <Typography variant="body1" color="text.secondary">
          Access all features and resources from here.
        </Typography>
      </Column>
      <Row sx={{ pb: 0, borderTop: 1, borderColor: 'divider' }}>
        <Typography variant="h5" color="primary" fontWeight={600}>
          Features
        </Typography>
      </Row>
      <Row>
        <Column sx={{ padding: 0 }}>
          <Typography variant="h4" component="h1">
            Gene Variation Analysis Tool
          </Typography>
          <Typography variant="body1" color="text.secondary">
            GVATool is an advanced genetic data processing tool designed to simplify the work of genetics 
            professionals when analyzing patients’ genetic information. The platform enables efficient detection 
            of critical genetic mutations, supports multiple file formats (.fasta, .fastq, .vcf), and converts 
            them into a standardized .vcf format. GVATool also integrates gene variation databases for mutation 
            analysis, offering interactive visualizations that help users interpret results and make critical healthcare decisions.
          </Typography>
          <Row>
            <Button
              variant="contained"
              startIcon={<PlagiarismIcon />}
              onClick={() => navigate('/features/gvatool')}
            >
              Open GVATool
            </Button>
          </Row>
        </Column>
      </Row>
      <Row sx={{ pb: 0, borderTop: 1, borderColor: 'divider' }}>
        <Typography variant="h5" color="primary" fontWeight={600}>
          System
        </Typography>
      </Row>
      <Row>
        <Column sx={{ padding: 0 }}>
          <Typography variant="h4" component="h1">
            File Manager
          </Typography>
          <Typography variant="body1" color="text.secondary">
            The File Manager is a tool for conveniently managing genetic data within the platform. It allows 
            users to upload, view, rename, delete, and organize genetic files in one place. The File Manager 
            simplifies working with large datasets, ensuring that genetics professionals can easily access and 
            manage the files they are analyzing, providing a clear and streamlined workflow.
          </Typography>
          <Row>
            <Button
              variant="contained"
              startIcon={<InventoryIcon />}
              onClick={() => navigate('/system/file_manager')}
            >
              Open File Manager
            </Button>
          </Row>
        </Column>
        <Column sx={{ padding: 0 }}>
          <Typography variant="h4" component="h1">
            Macros
          </Typography>
          <Typography variant="body1" color="text.secondary">
            The Macros functionality enables users to automate frequently repeated analysis steps. Genetics professionals 
            can record a sequence of actions - such as applying filters or performing mutation analyses and reuse them for 
            similar queries in the future. This feature saves time, reduces the risk of errors, and ensures consistency when 
            working with genetic data.
          </Typography>
          <Row>
            <Button
              variant="contained"
              startIcon={<AutoAwesomeIcon />}
              onClick={() => navigate('/system/macros')}
            >
              Open Macros
            </Button>
          </Row>
        </Column>
      </Row>
      <Row sx={{ pb: 0, borderTop: 1, borderColor: 'divider' }}>
        <Typography variant="h5" color="primary" fontWeight={600}>
          Resources
        </Typography>
      </Row>
      <Row>
        <Column sx={{ padding: 0 }}>
          <Typography variant="h4" component="h1">
            Manual
          </Typography>
          <Typography variant="body1" color="text.secondary">
            The Manual is a comprehensive user guide designed to help users quickly understand how to use GVATool's 
            features and workflows. It provides clear explanations on how to upload genetic data, perform mutation 
            analysis, apply filters, generate reports, and manage data files. The Manual is useful for both beginners 
            and experienced professionals, ensuring that every user can confidently and effectively use the platform.
          </Typography>
          <Row>
            <Button
              variant="contained"
              startIcon={<MenuBookIcon />}
              onClick={() => navigate('/resources/manual')}
            >
              Open Manual
            </Button>
          </Row>
        </Column>
      </Row>
    </Box>
  );
}
