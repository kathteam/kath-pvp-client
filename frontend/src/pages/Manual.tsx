import { JSX, Fragment, useEffect, useState } from 'react';
import { useLocation, useNavigate } from 'react-router-dom';
import {
  Plagiarism as PlagiarismIcon,
  Inventory as InventoryIcon,
  AutoAwesome as AutoAwesomeIcon,
  MenuBook as MenuBookIcon,
  ExpandMore as ExpandMoreIcon,
  ExpandLess as ExpandLessIcon,
} from '@mui/icons-material';
import { Collapse, Divider, IconButton, List, ListItemButton, ListItemIcon, Typography } from '@mui/material';
import { Column, Row, Container, Button, Header } from '@/components/core';
import { RouteHeader, Video } from '@/components';
import { VideoPlaceholder } from '@/assets';
import { handleScroll } from '@/utils';

const data = [
  {
    header: 'Features',
    items: [
      {
        title: 'Gene Variation Analysis Tool',
        icon: <PlagiarismIcon />,
        buttonText: 'Open GVATool',
        link: '/features/gvatool',
        description: 'GVATool is an advanced genetic data processing tool designed to simplify the work of genetics professionals when analyzing patients’ genetic information. The platform enables efficient detection of critical genetic mutations, supports multiple file formats (.fasta, .fastq, .vcf), and converts them into a standardized .vcf format. GVATool also integrates gene variation databases for mutation analysis, offering interactive visualizations that help users interpret results and make critical healthcare decisions.',
        subitems: [
          {
            text: 'Placeholder for video.',
            video: VideoPlaceholder,
          }
        ]
      },
    ]
  },
  {
    header: 'System',
    items: [
      {
        title: 'File Manager',
        icon: <InventoryIcon />,
        buttonText: 'Open File Manager',
        link: '/system/file_manager',
        description: 'The File Manager is a tool for conveniently managing genetic data within the platform. It allows users to upload, view, rename, delete, and organize genetic files in one place. The File Manager simplifies working with large datasets, ensuring that genetics professionals can easily access and manage the files they are analyzing, providing a clear and streamlined workflow.',
        subitems: [
          {
            text: 'Placeholder for video.',
            video: VideoPlaceholder,
          }
        ]
      },
      {
        title: 'Macros',
        icon: <AutoAwesomeIcon />,
        buttonText: 'Open Macros',
        link: '/system/macros',
        description: 'The Macros functionality enables users to automate frequently repeated analysis steps. Genetics professionals can record a sequence of actions - such as applying filters or performing mutation analyses and reuse them for similar queries in the future. This feature saves time, reduces the risk of errors, and ensures consistency when working with genetic data.',
        subitems: [
          {
            text: 'Placeholder for video.',
            video: VideoPlaceholder,
          }
        ]
      }
    ]
  }
];

export default function Manual(): JSX.Element {
  const navigate = useNavigate();
  const location = useLocation();
  
  const [isExpanded, setIsExpanded] = useState(false);

  useEffect(() => {
    const params = new URLSearchParams(location.search);
    const scrollTo = params.get('scrollTo');

    if (scrollTo) {
      handleScroll(scrollTo);
    }
  }, [location]);

  // Reset scroll position on initial load
  useEffect(() => {
    handleScroll('Manual');
  }, []);

  return (
    <Fragment>
      <RouteHeader
        icon={MenuBookIcon}
        title="Manual"
        description="A guide to help you explore the application's tools and features step by step."
      />

      <Column sx={{ borderTop: 1, borderColor: 'divider' }}>
        <Container>
          <Row sx={{ p: 0, justifyContent: 'space-between' }}>
            <Header variant="h5" fontWeight={600}>
              Contents
            </Header>
            <IconButton onClick={() => setIsExpanded((prev) => !prev)}>
              {isExpanded ? <ExpandLessIcon /> : <ExpandMoreIcon />}
            </IconButton>
          </Row>
          <Collapse in={isExpanded} timeout="auto" unmountOnExit>
            <List dense>
              {data.map((section) => (
                <Fragment key={section.header}>
                  <ListItemButton sx={{ py: 1 }} onClick={() => handleScroll(section.header)} >
                    <Typography variant='body1' fontWeight={700}>
                      {section.header}
                    </Typography>
                  </ListItemButton>
                  {section.items.map((item) => (
                    <ListItemButton
                      key={item.title}
                      sx={{ py: 1 }}
                      onClick={() => handleScroll(`${section.header} ${item.title}`)}
                    >
                      <ListItemIcon sx={{ minWidth: 40 }}>
                        {item.icon}
                      </ListItemIcon>
                      <Typography variant='body1'>
                        {item.title}
                      </Typography>
                    </ListItemButton>
                  ))}
                </Fragment>
              ))}
            </List>
          </Collapse>
        </Container>
      </Column>

      {data.map((section, idx) => (
        <Column key={`section${idx}`} id={section.header} sx={{ borderTop: 1, borderColor: 'divider' }}>
          <Container>
            <Row sx={{ p: 0 }}>
              <Header variant="h5" fontWeight={600}>
                {section.header}
              </Header>
            </Row>
            {section.items.map((item, itemidx) => (
              <Fragment key={`section${idx}-item${itemidx}`}>
                {itemidx > 0 && <Divider /> }
                <Row id={`${section.header} ${item.title}`} sx={{ px: 0, justifyContent: 'space-between' }}>
                  <Typography variant="h4" component="h1">
                    {item.title}
                  </Typography>
                  <Button
                    variant="contained"
                    startIcon={item.icon}
                    onClick={() => navigate(item.link)}
                  >
                    {item.buttonText}
                  </Button>
                </Row>
                <Row sx={{ px: 0, pt: 0 }}>
                  <Typography variant="body1">
                    {item.description}
                  </Typography>
                </Row>
                {item.subitems.map((subitem, subidx) => (
                  <Fragment key={`section${idx}-item${itemidx}-subitem${subidx}`}>
                    <Row sx={{ p: 0 }}>
                      <Typography variant="body1">
                        {subitem.text}
                      </Typography>
                    </Row>
                    <Video src={subitem.video} />
                  </Fragment>
                ))}
              </Fragment>
            ))}
          </Container>
        </Column>
      ))}
    </Fragment>
  );
}
