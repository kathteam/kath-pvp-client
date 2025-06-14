import React, { useEffect, useState } from 'react';
import {
  Container,
  Typography,
  Paper,
  Box,
  CircularProgress,
} from '@mui/material';

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
      <Paper elevation={4} sx={{ p: 5, borderRadius: 4, boxShadow: 6 }}>
        <Typography variant="h4" component="h1" gutterBottom>
          Analysis History
        </Typography>
        {loading ? (
          <Box sx={{ display: 'flex', justifyContent: 'center', my: 4 }}>
            <CircularProgress />
          </Box>
        ) : mutationEntries.length > 0 ? (
          <Box>
            {mutationEntries.map((entry, index) => (
              <Paper key={index} elevation={3} sx={{ p: 2, mb: 2 }}>
                <Typography variant="h6" gutterBottom>
                  {entry.clinical_significance} - {entry.disease_name}
                </Typography>
                <Typography variant="body2" color="text.secondary">
                  File Name: {entry.file_name}
                </Typography>
                <Typography variant="body2">
                  Chromosome: {entry.chromosome}, Position: {entry.position}, Ref: {entry.reference}, Alt: {entry.alternate}
                </Typography>
                <Typography variant="body2" color="text.secondary">
                  HGVS ID: {entry.hgvs_id}
                </Typography>
                <Typography variant="body2" color="primary">
                  <a
                    href={`https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr${entry.chromosome}%3A${entry.position}`}
                    target="_blank"
                    rel="noopener noreferrer"
                  >
                    View Genome
                  </a>
                </Typography>
              </Paper>
            ))}
          </Box>
        ) : (
          <Typography variant="body1" color="text.secondary">
            No mutation entries found.
          </Typography>
        )}
      </Paper>
    </Container>
  );
};

export default AnalysisHistory;