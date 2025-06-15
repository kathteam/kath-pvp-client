import { JSX, Fragment, useEffect } from 'react';
import { useNavigate } from 'react-router-dom';
import {
  Dashboard as DashboardIcon,
  Plagiarism as PlagiarismIcon,
  Inventory as InventoryIcon,
  AutoAwesome as AutoAwesomeIcon,
  MenuBook as MenuBookIcon,
  Storage as StorageIcon,
} from '@mui/icons-material';
import { Typography } from '@mui/material';
import { Button, Column, Row } from '@/components/core';
import { RouteHeader } from '@/components';
import { handleScroll } from '@/utils';

const sections = [
  {
    header: 'Features',
    items: [
      {
        title: 'Gene Variation Analysis Tool',
        icon: <PlagiarismIcon />,
        buttonText: 'Open GVATool',
        link: '/features/gvatool',
        description: `GVATool is an advanced genetic data processing tool designed to simplify the work of genetics 
          professionals when analyzing patients’ genetic information. The platform enables efficient detection 
          of critical genetic mutations, supports multiple file formats (.fasta, .fastq, .vcf), and converts 
          them into a standardized .vcf format. GVATool also integrates gene variation databases for mutation 
          analysis, offering interactive visualizations that help users interpret results and make critical healthcare decisions.`,
      },
      {
        title: 'Analysis History',
        icon: <StorageIcon />,
        buttonText: 'Open Analysis History',
        link: '/features/analysis_history',
        description: `The Analysis History feature provides a comprehensive log of all previously conducted analyses, 
          enabling users to easily track, review, and manage their past work. It offers quick access to detailed 
          records. This ensures reproducibility, simplifies  audit processes, and helps users maintain a clear workflow history 
          for informed decision-making and continuous improvement.`,
      }
    ],
  },
  {
    header: 'System',
    items: [
      {
        title: 'File Manager',
        icon: <InventoryIcon />,
        buttonText: 'Open File Manager',
        link: '/system/file_manager',
        description: `The File Manager is a tool for conveniently managing genetic data within the platform. It allows 
          users to upload, view, rename, delete, and organize genetic files in one place. The File Manager 
          simplifies working with large datasets, ensuring that genetics professionals can easily access and 
          manage the files they are analyzing, providing a clear and streamlined workflow.`,
      },
      {
        title: 'Macros',
        icon: <AutoAwesomeIcon />,
        buttonText: 'Open Macros',
        link: '/system/macros',
        description: `The Macros functionality enables users to automate frequently repeated analysis steps. Genetics professionals 
          can record a sequence of actions - such as applying filters or performing mutation analyses and reuse them for 
          similar queries in the future. This feature saves time, reduces the risk of errors, and ensures consistency when 
          working with genetic data.`,
      },
    ],
  },
  {
    header: 'Resources',
    items: [
      {
        title: 'Manual',
        icon: <MenuBookIcon />,
        buttonText: 'Open Manual',
        link: '/resources/manual',
        description: `The Manual is a comprehensive user guide designed to help users quickly understand how to use GVATool's 
          features and workflows. It provides clear explanations on how to upload genetic data, perform mutation 
          analysis, apply filters, generate reports, and manage data files. The Manual is useful for both beginners 
          and experienced professionals, ensuring that every user can confidently and effectively use the platform.`,
      },
    ],
  },
];

export default function Dashboard(): JSX.Element {
  const navigate = useNavigate();

  // Reset scroll position on initial load
  useEffect(() => {
    handleScroll('Welcome to the Dashboard!');
  }, []);

  return (
    <Fragment>
      <RouteHeader
        icon={DashboardIcon}
        title="Welcome to the Dashboard!"
        description="Access all features and resources from here."
      />

      {sections.map((section) => (
        <Fragment key={section.header}>
          <Row sx={{ pt: 2, pb: 0 }}>
            <Typography variant="h5" color="primary" fontWeight={600}>
              {section.header}
            </Typography>
          </Row>
          <Row sx={{ borderBottom: 1, borderColor: 'divider' }}>
            {section.items.map((item) => (
              <Column sx={{ padding: 0 }}>
                <Typography variant="h4" component="h2" gutterBottom>
                  {item.title}
                </Typography>
                <Typography variant="body1" color="text.primary" textAlign="justify">
                  {item.description}
                </Typography>
                <Button
                  variant="contained"
                  startIcon={item.icon}
                  onClick={() => navigate(item.link)}
                >
                  {item.buttonText}
                </Button>
              </Column>
            ))}
          </Row>
        </Fragment>
      ))}
    </Fragment>
  );
}
