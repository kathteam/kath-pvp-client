import React, { useEffect, useState } from 'react';
import {
  Container,
  Typography,
  Paper,
  Box,
  CircularProgress,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
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
          <TableContainer component={Paper}>
            <Table>
              <TableHead>
                <TableRow>
                  <TableCell>File Name</TableCell>
                  <TableCell>Clinical Significance</TableCell>
                  <TableCell>Disease Name</TableCell>
                  <TableCell>Synonyms</TableCell>
                  <TableCell>Chromosome</TableCell>
                  <TableCell>Position</TableCell>
                  <TableCell>Reference</TableCell>
                  <TableCell>Alternate</TableCell>
                  <TableCell>HGVS ID</TableCell>
                </TableRow>
              </TableHead>
              <TableBody>
                {mutationEntries.map((entry, index) => (
                  <TableRow key={index}>
                    <TableCell>{entry.file_name}</TableCell>
                    <TableCell>{entry.clinical_significance}</TableCell>
                    <TableCell>{entry.disease_name}</TableCell>
                    <TableCell>{entry.synonyms}</TableCell>
                    <TableCell>{entry.chromosome}</TableCell>
                    <TableCell>{entry.position}</TableCell>
                    <TableCell>{entry.reference}</TableCell>
                    <TableCell>{entry.alternate}</TableCell>
                    <TableCell>{entry.hgvs_id}</TableCell>
                  </TableRow>
                ))}
              </TableBody>
            </Table>
          </TableContainer>
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