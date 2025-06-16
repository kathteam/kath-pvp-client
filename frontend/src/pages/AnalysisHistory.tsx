import { Fragment, JSX, useEffect, useState } from 'react';
import {
  Typography,
  Box,
  CircularProgress,
  Chip,
  Card,
  CardContent,
  CardActions,
  Button,
  Divider,
  useTheme,
} from '@mui/material';
import {
  Dns as DnsIcon,
  Description as DescriptionIcon,
  LocationOn as LocationIcon,
  Link as LinkIcon,
  Storage as StorageIcon,
} from '@mui/icons-material';
import Ephasize from '@/components/text/Ephasize';
import { RouteHeader } from '@/components';
import { handleScroll } from '@/utils';
import { Column, Row, Button as ButtonCore } from '@/components/core';

export interface MutationEntry {
  file_name: string;
  clinical_significance: string;
  disease_name: string;
  synonyms: string;
  chromosome: string;
  position: string;
  reference: string;
  alternate: string;
  hgvs_id: string;
}

// Function to determine chip color based on clinical significance using semantic colors
const getSignificanceColor = (significance: string): 'error' | 'success' | 'warning' | 'default' => {
  const significanceLower = significance.toLowerCase();
  if (significanceLower.includes('pathogenic')) return 'error';
  if (significanceLower.includes('benign')) return 'success';
  if (significanceLower.includes('uncertain')) return 'warning';
  return 'default';
};

export default function AnalysisHistory(): JSX.Element {
  const theme = useTheme();
  const [mutationEntries, setMutationEntries] = useState<MutationEntry[]>([]);
  const [loading, setLoading] = useState<boolean>(true);
  const [pdfGeneratingStatus, setPdfGeneratingStatus] = useState<Record<string, 'idle' | 'processing' | 'success' | 'error'>>({});

  useEffect(() => {
    const fetchMutationEntries = async () => {
      try {
        const entries = await window.pywebview.api.file_controller.get_mutation_entries();
        setMutationEntries(entries);
      } catch (error) {
        console.error('Error fetching mutation entries:', error);
        setMutationEntries([]);
      } finally {
        setLoading(false);
      }
    };

    fetchMutationEntries();
  }, []);

  // Group mutationEntries by file_name
  const groupedEntries = mutationEntries.reduce((acc, entry) => {
    if (!acc[entry.file_name]) {
      acc[entry.file_name] = [];
    }
    acc[entry.file_name].push(entry);
    return acc;
  }, {} as Record<string, MutationEntry[]>);

  const handlePdfGeneration = async (entries: MutationEntry[], fileName: string) => {
    try {
      setPdfGeneratingStatus(prev => ({ ...prev, [fileName]: 'processing' }));
      await window.pywebview.api.ui_controller.generate_pdf(entries, `${fileName}_report.pdf`);
      setPdfGeneratingStatus(prev => ({ ...prev, [fileName]: 'success' }));
      window.pywebview.api.ui_controller.open_pdf_in_browser(`${fileName}_report.pdf`);
    } catch (error) {
      alert('Failed to generate PDF report. Please try again later.');
      setPdfGeneratingStatus(prev => ({ ...prev, [fileName]: 'error' }));
    }
  };

  // Reset scroll position on initial load
  useEffect(() => {
    handleScroll('Analysis History');
  }, []);

  return (
    <Fragment>
      <RouteHeader
        icon={StorageIcon}
        title="Analysis History"
        description="Review previous genetic analyses with detailed records, including results, parameters, and timestamps, ensuring easy tracking and reproducibility."
      />
      {loading ? (
        <Column>
          <CircularProgress color="primary" />
        </Column>
      ) : mutationEntries.length > 0 ? (
        <>
          {Object.entries(groupedEntries).map(([fileName, entries]) => (
            <Column key={fileName} sx={{ p: 0, borderBottom: 1, borderColor: 'divider', display: '-webkit-box', WebkitLineClamp: 2, WebkitBoxOrient: 'vertical', overflow: 'hidden', textOverflow: 'ellipsis', maxWidth: '100%' }}>
              <Row>
                <Typography variant="h5" sx={{ mb: 2, fontWeight: 600 }}>
                  {fileName}
                </Typography>
              </Row>
              <Row sx={{ py: 0 }}>
                <ButtonCore
                  key={fileName}
                  variant="contained"
                  startIcon={pdfGeneratingStatus[fileName] == 'processing' ? <CircularProgress size={20} color="inherit" /> : <DescriptionIcon/>}
                  onClick={() => handlePdfGeneration(entries, fileName)}
                >
                  Generate Report
                </ButtonCore>
                <ButtonCore
                  key={fileName}
                  variant="contained"
                  sx={{ visibility: pdfGeneratingStatus[fileName] === 'success' ? 'visible' : 'hidden' }}
                  onClick={() => window.pywebview.api.ui_controller.open_pdf_in_browser(`${fileName}_report.pdf`)}
                >
                  Open Report
                </ButtonCore>
              </Row>
              <Row sx={{ flexWrap: 'wrap' }}>
                {entries.map((entry, index) => (
                  <Card
                    key={index}
                    sx={{ 
                      height: '100%',
                      display: 'flex',
                      width: { xs: '100%', sm: '48%' },
                      flexDirection: 'column',
                      bgcolor: 'background.paper',
                      borderRadius: theme.shape.borderRadius
                    }}
                  >
                    <CardContent sx={{ flexGrow: 1 }}>
                      <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start', mb: 2, overflow: 'hidden', textOverflow: 'ellipsis', maxWidth: '100%'  }}>
                        <Typography 
                          variant="h6" 
                          gutterBottom 
                          component="div" 
                          sx={{ 
                            fontWeight: 500,
                            color: theme.palette.text.primary
                          }}
                        >
                          {entry.disease_name}
                        </Typography>
                        <Chip 
                          label={entry.clinical_significance} 
                          color={getSignificanceColor(entry.clinical_significance)}
                          size="small"
                          sx={{
                            fontWeight: 500,
                          }}
                        />
                      </Box>
                            
                      <Box sx={{ display: 'flex', alignItems: 'center', mb: 1 }}>
                        <DescriptionIcon fontSize="small" sx={{ mr: 1, color: theme.palette.text.primary }} />
                        <Ephasize label="File:" text={entry.file_name} />
                      </Box>
                            
                      <Box sx={{ display: 'flex', alignItems: 'center', mb: 1 }}>
                        <LocationIcon fontSize="small" sx={{ mr: 1, color: theme.palette.primary.main }} />
                        <Ephasize label="Chr:" text={entry.chromosome} />
                        <Ephasize label="Pos:" text={entry.position} />
                        <Ephasize label="Ref:" text={entry.reference} />
                        <Ephasize label="Alt:" text={entry.alternate} />
                      </Box>
                            
                      <Box sx={{ display: 'flex', alignItems: 'center', mb: 1 }}>
                        <DnsIcon fontSize="small" sx={{ mr: 1, color: theme.palette.info.main }} />
                        <Ephasize label="HGVS:" text={entry.hgvs_id} />
                      </Box>
                    </CardContent>
                          
                    <Divider sx={{ bgcolor: 'divider' }} />
                          
                    <CardActions sx={{ justifyContent: 'space-between', px: 2, py: 1.5 }}>
                      <Button 
                        size="small" 
                        startIcon={<LinkIcon />}
                        variant="outlined"
                        color="primary"
                        href={`https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr${entry.chromosome}%3A${entry.position}`}
                        target="_blank"
                        rel="noopener noreferrer"
                      >
                              UCSC
                      </Button>
                      <Button
                        size="small"
                        startIcon={<LinkIcon />}
                        variant="outlined"
                        color="primary"
                        href={`https://www.ensembl.org/Homo_sapiens/Location/View?db=core;r=${entry.chromosome}:${entry.position}`}
                        target="_blank"
                        rel="noopener noreferrer"
                      >
                              Ensembl
                      </Button>
                    </CardActions>
                  </Card>
                ))}
              </Row>
            </Column>
          ))}
        </>
      ) : (
        <Column sx={{ textAlign: 'center', py: 4, borderBottom: 1, borderColor: 'divider' }}>
          <Typography variant="body1" sx={{ fontWeight: 500, color: 'text.secondary' }}>
            No mutation entries found.
          </Typography>
        </Column>
      )}
    </Fragment>
  );
};
