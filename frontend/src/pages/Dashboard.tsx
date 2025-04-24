import { JSX, useEffect, useRef, useState } from 'react';
import {
  Box,
  Button,
  Container,
  List,
  ListItem,
  ListItemButton,
  Typography,
  useTheme,
} from '@mui/material';
import { Folder as FolderIcon, Article as FileIcon, ExpandMore as ExpandedIcon, ChevronRight as HiddenIcon } from '@mui/icons-material';
import { CodeLine } from '@/components';
import React from 'react';

interface File {
  name: string;
  type: 'fasta' | 'vcf' | 'pdf';
  path: string;
  content: string;
}

interface Directory {
  name: string;
  files: File[];
}

const mockDirectories: Directory[] = [
  {
    name: 'sample1',
    files: [
      {
        name: 'data.fasta',
        type: 'fasta',
        path: '/sample1/data.fasta',
        content: '>seq1 Homo sapiens chromosome 1 genomic scaffold\nATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC\nGATCGTAGCTAGCTAGCATCGATCGTAGCTAGCTAGCGTACGTAGCTAGCTAGCTAGCTAGCTA\n\n>seq2 Homo sapiens mitochondrial DNA complete genome\nTTAGCGATCGTAGCTACGTAGCTAGCGTACGTAGCTAGCATCGTAGCTAGCTAGCTAGCTAGCG\nATCGTACGTAGCATCGTAGCATCGTAGCTAGCGTACGTAGCATCGTAGCTAGCTAGCTAGCTAG\n\n>seq3 Mus musculus chromosome 2 genomic scaffold\nGGATCGTACGATCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\nTTCGATCGTAGCTAGCATCGTAGCTAGCTAGCTAGCGTACGTAGCATCGTAGCATCGTAGCTAG\n'
      },
      {
        name: 'mutations.vcf',
        type: 'vcf',
        path: '/sample1/mutation.vcf',
        content: '##fileformat=VCFv4.2\n##source=MockVariantCaller\n##reference=GRCh38\n##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n1\t1234567\trs11111111\tG\tA\t50.3\tPASS\tDP=40\tGT\t0/1\n1\t2345678\t.\tT\tC\t20.1\tq10\tDP=22	GT\t0/0\n2\t3456789\trs22222222\tC\tT\t99.9\tPASS\tDP=55\tGT\t1/1\n3\t4567890\trs33333333\tA\tG\t45.7\tPASS\tDP=30\tGT\t0/1\n'
      },
      {
        name: 'diseases.vcf',
        type: 'vcf',
        path: '/sample1/diseases.vcf',
        content: '##fileformat=VCFv4.2\n##source=MockVariantCaller\n##reference=GRCh38\n##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n##INFO=<ID=CLNDN,Number=.,Type=String,Description="ClinVar Disease Name">\n##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical Significance">\n##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">\n##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n1\t1234567\trs11111111\tG\tA\t50.3\tPASS\tDP=40;GENE=CFTR;CLNDN=Cystic_fibrosis;CLNSIG=Pathogenic;AF=0.01\tGT\t0/1\n1\t2345678\t.\tT\tC\t20.1\tq10\tDP=22;GENE=BRCA1;CLNDN=Breast_cancer;CLNSIG=Likely_pathogenic;AF=0.005\tGT\t0/0\n2\t3456789\trs22222222\tC\tT\t99.9\tPASS\tDP=55;GENE=TP53;CLNDN=Li-Fraumeni_syndrome;CLNSIG=Pathogenic;AF=0.0001\tGT\t1/1\n3\t4567890\trs33333333\tA\tG\t45.7\tPASS\tDP=30;GENE=MYH7;CLNDN=Cardiomyopathy;CLNSIG=Uncertain_significance;AF=0.02\tGT\t0/1\n'
      },
      {
        name: 'report.pdf',
        type: 'pdf',
        path: '/sample1/report.pdf',
        content: 'This is a simulated PDF document. Actual content is binary.'
      }
    ]
  },
  {
    name: 'sample2',
    files: [
      {
        name: 'data.fasta',
        type: 'fasta',
        path: '/sample2/data.fasta',
        content: '>ecoli_K12_complete_genome\nATGACCATGATTACGCCAAGCTATTCAGATGAGCCTAGGTATTACGGTGTCGTTGTTGAACTGGA\nAGTGGCACAGAGTTATCGGATCCAGATGAAGTTCGTGAGCGGATAACAATTTCACACAGGAAAC\n\n>saccharomyces_cerevisiae_chrIV\nTTCTGACCTGACTCGAGTCCGTAATTCCTGGCCTGACTTGGTGAAACTGACCTGAGTGTTGGTTG\nTGCATAGGCGAGGTTGAAATGGTCTTGGCATGGCTGTTAGTGCTAGTGTTTTGTTAGCGTGAGA\n\n>arabidopsis_thaliana_chr1_region\nGGCGTCTGCTTATGGTGGATGTGAACGTCTGTTTCCTGTCCTGAGTCTGGTAGCTGATGGTGGAT\nCTTTCATCAGGTGATGGTGATCTTGTAGCTGATGGTGGTGATCTAGTGGTGATGGTGATGGTGGA\n\n>human_chr7_partial\nCAGTGGTGATGGAGTGTGTGTGCTTCTTCTTGGAGGAGTGGTGAGTGGTGGATGTGTGTGTGTGT\nGTGAGTGGTGGAGTGGAGGAGGTGGTGGTGAGTGGTGAGTGGTGGAGTGGTGGAGTGGTGAGTGG\n'
      }
    ]
  }
];

export default function Dashboard(): JSX.Element {
  const theme = useTheme();

  const scrollRefLeft = useRef<HTMLDivElement>(null);
  const scrollRefRight = useRef<HTMLDivElement>(null);

  const [expandedDirectories, setExpandedDirectories] = useState<Set<string>>(new Set());
  const [selectedDirectory, setSelectedDirectory] = useState<Directory | null>(null);
  const [selectedFile, setSelectedFile] = useState<File | null>(null);
  const [nextFile, setNextFile] = useState<File | null>(null);

  const handleToggleDirectory = (directoryName: string) => {
    setExpandedDirectories((prev) => {
      const newExpandedDirectories = new Set(prev);
      if (newExpandedDirectories.has(directoryName)) {
        newExpandedDirectories.delete(directoryName);
      } else {
        newExpandedDirectories.add(directoryName);
      }
      return newExpandedDirectories;
    });
  };

  const handleFileSelection = (file: File) => {
    setSelectedFile(file);
    setSelectedDirectory(
      mockDirectories.find((dir) =>
        dir.files.some((f) => f.path === file.path)
      ) || null
    );
  };

  useEffect(() => {
    const getNextFile = (file: File): File | null => {
      if (!selectedDirectory) return null;
  
      const currentIndex = selectedDirectory.files.findIndex((f) => f.path === file.path);
      return selectedDirectory.files[currentIndex + 1] || null;
    };

    if (selectedFile)
      setNextFile(() => getNextFile(selectedFile));
    else setNextFile(null);

    if (scrollRefLeft.current) {
      scrollRefLeft.current.scrollTop = 0;
    }

    if (scrollRefRight.current) {
      scrollRefRight.current.scrollTop = 0;
    }
  }, [selectedFile, selectedDirectory]);

  return (
    <Container
      maxWidth={false}
      style={{
        padding: 0,
        display: 'flex',
        flexDirection: 'row',
        height: '100%',
        overflow: 'hidden',
        backgroundColor: theme.palette.background.paper,
      }}
    >
      {/* Left Panel: File Tree */}
      <Box
        sx={{
          flex: 1,
          maxWidth: '270px',
          borderRight: `1px solid ${theme.palette.divider}`,
          overflowY: 'auto',
        }}
      >
        <Typography
          variant="h6"
          sx={{
            p: 2,
            borderBottom: `1px solid ${theme.palette.divider}`
          }}
        >
          Directories & Files
        </Typography>
        <List sx={{ padding: 0 }}>
          {mockDirectories.map((directory) => (
            <React.Fragment key={directory.name}>
              <ListItem
                disablePadding
                sx={{
                  backgroundColor: directory === selectedDirectory
                    ? theme.palette.background.default
                    : 'transparent',
                }}
              >
                <ListItemButton 
                  disableGutters 
                  onClick={() => handleToggleDirectory(directory.name)}
                >
                  {expandedDirectories.has(directory.name) 
                    ? <ExpandedIcon sx={{ color: theme.palette.text.secondary }}/>
                    : <HiddenIcon sx={{ color: theme.palette.text.primary }}/>}
                  <FolderIcon
                    sx={{
                      marginRight: 1,
                      color: theme.palette.primary.dark
                    }}
                  />
                  <Typography
                    variant="subtitle1"
                    fontWeight="bold"
                    sx={{
                      wordBreak: 'break-all'
                    }}
                  >
                    {directory.name}
                  </Typography>
                </ListItemButton>
              </ListItem>

              {/* Show files only if the directory is expanded */}
              {expandedDirectories.has(directory.name) && directory.files.map((file) => (
                <ListItem
                  key={file.path}
                  disablePadding
                  sx={{
                    backgroundColor: file === selectedFile 
                      ? theme.palette.action.selected
                      : directory === selectedDirectory
                        ? theme.palette.background.default
                        : 'transparent'
                  }}
                >
                  <ListItemButton 
                    disableGutters
                    onClick={() => handleFileSelection(file)}
                  >
                    <FileIcon
                      sx={{
                        marginLeft: 5,
                        marginRight: 1,
                        color: theme.palette.text.secondary
                      }}
                    />
                    <Typography
                      sx={{
                        wordBreak: 'break-all',
                      }}
                    >
                      {file.name}
                    </Typography>
                  </ListItemButton>
                </ListItem>
              ))}
            </React.Fragment>
          ))}
        </List>
      </Box>

      <Box
        sx={{
          flex: 1,
          display: 'flex',
          flexDirection: 'row',
          overflow: 'hidden',
        }}
      >
        {/* Center Panel: Selected File Content */}
        <Box
          sx={{
            flex: 1,
            borderRight: `1px solid ${theme.palette.divider}`,
            display: 'flex',
            flexDirection: 'column',
            bgcolor: theme.palette.background.default,
          }}
        >
          <Typography
            variant="h6"
            sx={{
              p: 2,
              borderBottom: `1px solid ${theme.palette.divider}`
            }}
          >
            {selectedFile ? selectedFile.name : 'Select a File'}
          </Typography> 
          <Box
            ref={scrollRefLeft} 
            sx={{
              borderRadius: 0,
              boxShadow: 'none',
              flex: 1,
              whiteSpace: 'pre-wrap',
              overflowY: 'auto',
              bgcolor: 'none',
            }}>
            {selectedFile ? (
              (() => {
                return selectedFile.content
                  .split('\n')
                  .map((line, index) => (
                    <CodeLine 
                      key={index} lineNumber={index + 1}
                      content={line}
                      bgcolor={theme.palette.background.default}
                    />
                  ));
              })()
            ) : (
              'No file selected'
            )}
          </Box>
        </Box>

        {/* Right Panel: Next File Content */}
        <Box
          sx={{
            flex: 1,
            display: 'flex',
            flexDirection: 'column',
            bgcolor: theme.palette.background.default,
          }}
        >
          <Typography
            variant="h6"
            sx={{
              p: 2,
              borderBottom: `1px solid ${theme.palette.divider}`
            }}
          >
            {nextFile ? nextFile.name : selectedFile ? 'No File Available' : 'Select a File to see Next File'}
          </Typography> 
          <Box
            ref={scrollRefRight}
            sx={{
              boxShadow: 'none',
              borderRadius: 0,
              flex: 1,
              whiteSpace: 'pre-wrap',
              overflowY: 'auto'
            }}
          >
            {nextFile ? (
              (() => {
                return nextFile.content
                  .split('\n')
                  .map((line, index) => (
                    <CodeLine 
                      key={index}
                      lineNumber={index + 1}
                      content={line}
                      bgcolor={theme.palette.background.default}
                    />
                  ));
              })()
            ) : selectedFile ? (
              <Box 
                sx={{
                  alignContent: 'center',
                  textAlign: 'center'
                }}
              >
                <Button
                  variant="contained"
                  onClick={() => alert('Generate next file')}
                >
                  Generate Next File
                </Button>
              </Box>
            ) : (
              <Typography>Select a file to see what comes next</Typography>
            )}
          </Box>
        </Box>
      </Box>
    </Container>
  );
}
