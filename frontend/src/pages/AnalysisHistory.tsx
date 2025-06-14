import React, { useEffect, useState } from 'react';
import {
  Container,
  Typography,
  Paper,
  Box,
  CircularProgress,
  Grid,
  Chip,
  Card,
  CardContent,
  CardActions,
  Button,
  Divider,
} from '@mui/material';
import {
  Dns as DnsIcon,
  Description as DescriptionIcon,
  LocationOn as LocationIcon,
  Link as LinkIcon,
} from '@mui/icons-material';

interface MutationEntry {
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

// Function to determine chip color based on clinical significance
const getSignificanceColor = (significance: string): string => {
  const significanceLower = significance.toLowerCase();
  if (significanceLower.includes('pathogenic')) return 'error';
  if (significanceLower.includes('benign')) return 'success';
  if (significanceLower.includes('uncertain')) return 'warning';
  return 'default';
};

const AnalysisHistory: React.FC = () => {
  const [mutationEntries, setMutationEntries] = useState<MutationEntry[]>([]);
  const [loading, setLoading] = useState<boolean>(true);

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

  return (
    <Container maxWidth="lg" sx={{ py: 4 }}>
      <Paper elevation={4} sx={{ p: { xs: 2, md: 5 }, borderRadius: 4, boxShadow: 6 }}>
        <Typography variant="h4" component="h1" gutterBottom sx={{ mb: 4, fontWeight: 500 }}>
          Analysis History
        </Typography>
        {loading ? (
          <Box sx={{ display: 'flex', justifyContent: 'center', my: 4 }}>
            <CircularProgress />
          </Box>
        ) : mutationEntries.length > 0 ? (
          <Grid container spacing={3}>
            {mutationEntries.map((entry, index) => (
              <Grid item xs={12} md={6} key={index}>
                <Card elevation={3} sx={{ height: '100%', display: 'flex', flexDirection: 'column' }}>
                  <CardContent sx={{ flexGrow: 1 }}>
                    <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'flex-start', mb: 2 }}>
                      <Typography variant="h6" gutterBottom component="div" sx={{ fontWeight: 500 }}>
                        {entry.disease_name}
                      </Typography>
                      <Chip 
                        label={entry.clinical_significance} 
                        color={getSignificanceColor(entry.clinical_significance) as any}
                        size="small"
                      />
                    </Box>
                    
                    <Box sx={{ display: 'flex', alignItems: 'center', mb: 1 }}>
                      <DescriptionIcon fontSize="small" sx={{ mr: 1, color: 'text.secondary' }} />
                      <Typography variant="body2" color="text.secondary" sx={{ fontWeight: 500 }}>
                        File: {entry.file_name}
                      </Typography>
                    </Box>
                    
                    <Box sx={{ display: 'flex', alignItems: 'center', mb: 1 }}>
                      <LocationIcon fontSize="small" sx={{ mr: 1, color: 'text.secondary' }} />
                      <Typography variant="body2" sx={{ fontWeight: 500 }}>
                        Chr: {entry.chromosome}, Pos: {entry.position}, Ref: {entry.reference}, Alt: {entry.alternate}
                      </Typography>
                    </Box>
                    
                    <Box sx={{ display: 'flex', alignItems: 'center', mb: 1 }}>
                      <DnsIcon fontSize="small" sx={{ mr: 1, color: 'text.secondary' }} />
                      <Typography variant="body2" color="text.secondary" sx={{ fontWeight: 500 }}>
                        HGVS: {entry.hgvs_id}
                      </Typography>
                    </Box>
                  </CardContent>
                  
                  <Divider />
                  
                  <CardActions sx={{ justifyContent: 'space-between', px: 2, py: 1.5 }}>
                    <Button 
                      size="small" 
                      startIcon={<LinkIcon />}
                      variant="outlined"
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
                      href={`https://www.ensembl.org/Homo_sapiens/Location/View?db=core;r=${entry.chromosome}:${entry.position}`}
                      target="_blank"
                      rel="noopener noreferrer"
                    >
                      Ensembl
                    </Button>
                  </CardActions>
                </Card>
              </Grid>
            ))}
          </Grid>
        ) : (
          <Box sx={{ textAlign: 'center', py: 4 }}>
            <Typography variant="body1" color="text.secondary" sx={{ fontWeight: 500 }}>
              No mutation entries found.
            </Typography>
          </Box>
        )}
      </Paper>
    </Container>
  );
};

export default AnalysisHistory;