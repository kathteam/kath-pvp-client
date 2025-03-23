"""
FASTA Downloader Module for bioinformatics data retrieval.

This module provides functions to download FASTA files from various bioinformatics
databases including NCBI, Ensembl, and ClinVar.
"""
import os
import sys
import time
import gzip
import shutil
import logging
import subprocess
from pathlib import Path
from typing import Optional, List, Dict, Tuple


# Third-party imports
import requests
from Bio import Entrez
from dotenv import load_dotenv
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def initialize_directories():
    """Initialize directory structure for the application."""
    # Set up file paths
    base_dir = os.path.join(os.path.expanduser("~"), ".KATH")
    service_dir = os.path.join(base_dir, "environment")
    fasta_dir = os.path.join(base_dir, "shared", "data", "fasta_files")
    blast_dir = os.path.join(base_dir, "shared", "data", "blast_results")
    uploads_dir = os.path.join(fasta_dir, "uploads")
    ref_dir = os.path.join(fasta_dir, "reference")
    samples_dir = os.path.join(fasta_dir, "samples")

    # Create directories if they don't exist
    for dir_path in [service_dir, service_dir, fasta_dir, blast_dir, uploads_dir, ref_dir, samples_dir]:
        os.makedirs(dir_path, exist_ok=True)
    
    # Load environment variables
    env_file = os.path.join(service_dir, ".env")
    if os.path.exists(env_file):
        load_dotenv(env_file)
    else:
        logger.warning("No .env file found. Generating a sample .env file")
        # Create a sample .env file
        with open(env_file, "w") as f:
            f.write("NCBI_API_KEY=\n")
            f.write("NCBI_API_EMAIL=\n")
    
    return {
        "base_dir": base_dir,
        "service_dir": service_dir, 
        "fasta_dir": fasta_dir,
        "blast_dir": blast_dir,
        "uploads_dir": uploads_dir,
        "ref_dir": ref_dir,
        "samples_dir": samples_dir
    }

# Default timeout for HTTP requests (in seconds)
DEFAULT_TIMEOUT = 60


def blast_cmdline(query_fasta_path: str, reference_genome_path: str, output_dir: Optional[Path] = None) -> str:
    """
    Run BLAST using direct command-line execution of BLAST+ tools.
    
    Args:
        query_fasta_path: Path to query FASTA file
        reference_genome_path: Path to reference genome
        output_dir: Directory to save results
        
    Returns:
        Path to BLAST results file
    """
    # Extract filename from query path
    query_filename = Path(query_fasta_path).stem
    reference_filename = Path(reference_genome_path).stem
    output_file = os.path.join(output_dir, f"{query_filename}_vs_{reference_filename}_blast.xml")
    
    logger.info(f"Running BLAST+ alignment for {query_filename} against {reference_filename}...")
    
    # Check if database exists, if not create it
    db_path = os.path.join(Path(reference_genome_path).parent, reference_filename)
    db_files = list(Path(reference_genome_path).parent.glob(f"{reference_filename}.n*"))
    
    if not db_files:
        logger.info(f"Creating BLAST database from {reference_genome_path}...")
        makeblastdb_cmd = [
            "makeblastdb",
            "-in", str(reference_genome_path),
            "-dbtype", "nucl",
            "-parse_seqids",
            "-out", str(db_path)
        ]
        
            # Redirect output to devnull to suppress console output
        with open(os.devnull, 'w') as devnull:
            process = subprocess.run(
                makeblastdb_cmd,
                stdout=devnull,
                stderr=subprocess.PIPE,
                text=True,
                check=False
            )
        
        process = subprocess.run(
            makeblastdb_cmd,
            capture_output=True,
            text=True,
            check=False
        )
        
        if process.returncode != 0:
            logger.error(f"Failed to create BLAST database: {process.stderr}")
            raise RuntimeError(f"makeblastdb failed with exit code {process.returncode}")
        
        logger.info("BLAST database created successfully")
    
    # Run BLAST alignment
    blastn_cmd = [
        "blastn",
        "-query", str(query_fasta_path),
        "-db", str(db_path),
        "-out", str(output_file),
        "-outfmt", "5",  # XML output format
        "-evalue", "1e-5",  # Less stringent e-value for faster results
        "-word_size", "28",  # Larger word size speeds up search significantly
        "-max_target_seqs", "5",
        "-num_threads", str(os.cpu_count() or 1),
        "-task", "megablast"  # Use fast mode
    ]
    
    logger.info(f"Executing BLAST command: {' '.join(blastn_cmd)}")
    
    process = subprocess.run(
        blastn_cmd,
        capture_output=True,
        text=True,
        check=False
    )
    
    if process.returncode != 0:
        logger.error(f"BLAST alignment failed: {process.stderr}")
        raise RuntimeError(f"blastn failed with exit code {process.returncode}")
    
    logger.info("BLAST alignment completed successfully")
    
    # Parse results using BioPython's NCBIXML parser
    parse_blast_results(output_file)
    
    return str(output_file)

def extract_chromosome(subject_id: str) -> str:
    """
    Extract chromosome information from the subject ID.
    
    Args:
        subject_id: The subject ID string from BLAST results
        
    Returns:
        Chromosome identifier or "Unknown"
    """
    # Common formats: 
    # - NC_000001.11 (Chromosome 1)
    # - chr1, chrX, chrY, etc.
    # - 1, 2, 3, etc. (just number)
    
    import re
    # Try to find chromosome in format "chr1" or "chrX"
    chr_match = re.search(r'chr([0-9XYMxym]+)', subject_id)
    if chr_match:
        return chr_match.group(1).upper()
    
    # Try to find NC_0000XX format (NCBI RefSeq accessions)
    nc_match = re.search(r'NC_0000(\d{2})', subject_id)
    if nc_match:
        num = int(nc_match.group(1))
        if 1 <= num <= 22:
            return str(num)
        elif num == 23:
            return 'X'
        elif num == 24:
            return 'Y'
    
    # Try to extract just a number if it's at the beginning or isolated
    num_match = re.search(r'\b([0-9]+|[XYxy])\b', subject_id)
    if num_match:
        return num_match.group(1).upper()
    
    return "Unknown"

def parse_blast_results(result_file: str):
    """Parse and log BLAST results from XML file."""
    logger.info(f"Parsing BLAST results from {result_file}")
    
    with open(result_file) as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    # Extract subject/reference ID
                    subject_id = alignment.title.split()[1] if len(alignment.title.split()) > 1 else alignment.title
                    
                    # Extract chromosome information
                    chromosome = extract_chromosome(subject_id)
                    
                    logger.info(f"Sequence: {alignment.title}")
                    logger.info(f"  Chromosome: {chromosome}")
                    logger.info(f"  Length: {alignment.length}")
                    logger.info(f"  E-value: {hsp.expect}")
                    logger.info(f"  Score: {hsp.score}")
                    logger.info(f"  Identities: {hsp.identities}/{hsp.align_length} "
                               f"({hsp.identities/hsp.align_length*100:.1f}%)")
                    
                    # Report mutations (mismatches)
                    mutations = []
                    for i, (q, s) in enumerate(zip(hsp.query, hsp.sbjct)):
                        if q != s and q != '-' and s != '-':
                            position = hsp.sbjct_start + i
                            mutations.append({
                                'chromosome': chromosome,
                                'position': position,
                                'reference': s,
                                'query': q
                            })
                    
                    if mutations:
                        logger.info(f"  Found {len(mutations)} mutations")
                        for mut in mutations[:5]:  # Show first 5 mutations
                            logger.info(f"    Chr {mut['chromosome']}, Position {mut['position']}: {mut['reference']} -> {mut['query']}")
                        
                        if len(mutations) > 5:
                            logger.info(f"    ... and {len(mutations) - 5} more mutations")


# Example usage
if __name__ == "__main__":
    load_dotenv()  # Make sure environment variables are loaded
    PRIVATE_API = os.environ.get("NCBI_API_KEY")
    PRIVATE_EMAIL = os.environ.get("NCBI_API_EMAIL")
    
    dir = initialize_directories()
    reference_path = dir["ref_dir"]
    
    # Use command line BLAST+
    try:
        # Check if BLAST+ is installed
        subprocess.run(["blastn", "-version"], capture_output=True, check=True)
        logger.info("BLAST+ tools are available on this system")
        
        # First check if there's a reference genome FASTA file
        ref_files = [f for f in os.listdir(reference_path) if f.endswith('.fasta') or f.endswith('.fa')]
        if not ref_files:
            logger.error(f"No reference genome FASTA files found in {reference_path}")
            sys.exit(1)
        
        # Use the first reference file found
        reference_fasta = os.path.join(reference_path, ref_files[0])
        logger.info(f"Using reference genome: {reference_fasta}")
        
        user_uploads_folder = Path(dir["uploads_dir"])
        print(user_uploads_folder)
        if user_uploads_folder.exists():
            sample_file_name = os.listdir(user_uploads_folder)
            
            for file in sample_file_name:
                if file.endswith(".fasta"):
                    sample_path = os.path.join(user_uploads_folder, file)
                    print(sample_path)
            
            # # Run BLAST using command line
                result_file = blast_cmdline(str(sample_path), str(reference_fasta), str(dir["blast_dir"]))
                logger.info(f"BLAST analysis complete. Results saved to {result_file}")
        else:
            logger.error(f"Sample directory does not exist: {user_uploads_folder}")
            
    except subprocess.CalledProcessError:
        logger.error("BLAST+ tools not found. Please install BLAST+ (sudo apt install ncbi-blast+)")
    except Exception as e:
        logger.exception(f"Error running BLAST: {e}")
        