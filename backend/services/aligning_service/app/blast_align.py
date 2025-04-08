"""
FASTA Downloader Module for bioinformatics data retrieval.

This module provides functions to download FASTA files from various bioinformatics
databases including NCBI, Ensembl, and ClinVar.
"""

import os
import subprocess
import sys
from pathlib import Path
from typing import Optional

from Bio.Blast import NCBIXML

backend_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../"))
if backend_dir not in sys.path:
    sys.path.append(backend_dir)

from shared.constants import (
    PROGRAM_STORAGE_DIR_SHARED_BLAST,
    PROGRAM_STORAGE_DIR_SHARED_DATA_FASTA,
)
from utils.logger import get_logger


logger = get_logger(__name__)


# Default timeout for HTTP requests (in seconds)
DEFAULT_TIMEOUT = 60


def blast_cmdline(
    query_fasta_path: str, reference_genome_path: str, output_dir: Optional[Path] = None
) -> str:
    """
    Run BLAST using direct command-line execution of BLAST+ tools.

    Args:
        query_fasta_path: Path to query FASTA file
        reference_genome_path: Path to reference genome
        output_dir: Directory to save results

    Returns:
        Path to BLAST results file
    """
    query_filename = Path(query_fasta_path).stem
    reference_filename = Path(reference_genome_path).stem
    output_file = os.path.join(
        output_dir, f"{query_filename}_vs_{reference_filename}_blast.xml"
    )

    logger.info(
        f"Running BLAST+ alignment for {query_filename} against {reference_filename}..."
    )

    db_path = os.path.join(Path(reference_genome_path).parent, reference_filename)
    db_files = list(Path(reference_genome_path).parent.glob(f"{reference_filename}.n*"))

    if not db_files:
        makeblastdb_cmd = [
            "makeblastdb",
            "-in",
            str(reference_genome_path),
            "-dbtype",
            "nucl",
            "-parse_seqids",
            "-out",
            str(db_path),
        ]

        process = subprocess.run(
            makeblastdb_cmd, capture_output=True, text=True, check=False
        )

        if process.returncode != 0:
            logger.error(f"Failed to create BLAST database: {process.stderr}")
            raise RuntimeError(
                f"makeblastdb failed with exit code {process.returncode}"
            )

    # Run BLAST alignment with stricter parameters
    blastn_cmd = [
        "blastn",
        "-query",
        str(query_fasta_path),
        "-db",
        str(db_path),
        "-out",
        str(output_file),
        "-outfmt",
        "5",  # XML output format
        "-evalue",
        "1e-10",  # More stringent e-value for higher quality matches
        "-word_size",
        "28",  # Larger word size speeds up search significantly
        "-max_target_seqs",
        "3",  # Limit to best matches only
        "-perc_identity",
        "95",  # Higher identity requirement (95%)
        "-num_threads",
        str(os.cpu_count() or 1),
        "-task",
        "blastn",  # Use fast mode
    ]

    process = subprocess.run(blastn_cmd, capture_output=True, text=True, check=False)

    if process.returncode != 0:
        logger.error(f"BLAST alignment failed: {process.stderr}")
        raise RuntimeError(f"blastn failed with exit code {process.returncode}")

    logger.info("BLAST alignment completed successfully")

    mutations = parse_blast_results(output_file, min_score=100, min_identity=95, max_mutations=10)

    mutations_file = output_file.replace(".xml", "_significant_mutations.json")
    if mutations:
        import json
        with open(mutations_file, 'w') as f:
            json.dump(mutations, f, indent=2)
        logger.info(f"Saved {len(mutations)} significant mutations to {mutations_file}")

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

    chr_match = re.search(r"chr([0-9XYMxym]+)", subject_id)
    if chr_match:
        return chr_match.group(1).upper()

    # Try to find NC_0000XX format (NCBI RefSeq accessions)
    nc_match = re.search(r"NC_0000(\d{2})", subject_id)
    if nc_match:
        num = int(nc_match.group(1))
        if 1 <= num <= 22:
            return str(num)
        elif num == 23:
            return "X"
        elif num == 24:
            return "Y"

    num_match = re.search(r"\b([0-9]+|[XYxy])\b", subject_id)
    if num_match:
        return num_match.group(1).upper()

    return "Unknown"


def parse_blast_results(result_file: str, min_score=100, min_identity=95, max_mutations=10):
    """
    Parse and log BLAST results from XML file with filtering.
    
    Args:
        result_file: Path to BLAST XML output file
        min_score: Minimum bit score to consider an alignment significant
        min_identity: Minimum percent identity to consider an alignment (0-100)
        max_mutations: Maximum number of mutations to report per alignment
    """
    significant_mutations = []
    
    with open(result_file) as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    percent_identity = (hsp.identities / hsp.align_length) * 100
                    if hsp.bits < min_score or percent_identity < min_identity:
                        continue
                    
                    subject_id = (
                        alignment.title.split()[1]
                        if len(alignment.title.split()) > 1
                        else alignment.title
                    )

                    chromosome = extract_chromosome(subject_id)

                    mutations = []
                    for i, (q, s) in enumerate(zip(hsp.query, hsp.sbjct)):
                        if q != s and q != "-" and s != "-":
                            position = hsp.sbjct_start + i
                            
                            rel_position = i / hsp.align_length
                            
                            if rel_position < 0.05 or rel_position > 0.90:
                                continue
                                
                            mutations.append(
                                {
                                    "chromosome": chromosome,
                                    "position": position,
                                    "reference": s,
                                    "query": q,
                                    "quality": hsp.bits,
                                    "percent_identity": percent_identity
                                }
                            )
                    
                    if 0 < len(mutations) <= max_mutations:
                        significant_mutations.extend(mutations)
                        logger.info(
                            f"Found {len(mutations)} significant mutations in alignment "
                            f"with {percent_identity:.1f}% identity and bit score {hsp.bits:.1f}"
                        )
                    elif len(mutations) > max_mutations:
                        logger.info(
                            f"Skipped {len(mutations)} mutations in low-quality alignment "
                            f"(exceeds threshold of {max_mutations})"
                        )
    
    if significant_mutations:
        logger.info(f"Total significant mutations found: {len(significant_mutations)}")
    else:
        logger.info("No significant mutations found after filtering")
    
    return significant_mutations


def perform_blast_aligning(query_fasta_path: str = None) -> str:
    """
    Perform BLAST analysis using the provided query and reference genome.

    Args:
        query_fasta_path (str, optional): Path to the specific FASTA file to use as query.
            If None, will attempt to use files from the uploads directory.

    Returns:
        Path to BLAST results file
    """

    reference_path = os.path.join(PROGRAM_STORAGE_DIR_SHARED_DATA_FASTA, "reference")
    blast_dir = PROGRAM_STORAGE_DIR_SHARED_BLAST

    try:
        # Check if BLAST+ is installed
        subprocess.run(["blastn", "-version"], capture_output=True, check=True)
        logger.info("BLAST+ tools are available on this system")

        # First check if there's a reference genome FASTA file
        ref_files = [
            f
            for f in os.listdir(reference_path)
            if f.endswith(".fasta") or f.endswith(".fa")
        ]
        if not ref_files:
            logger.error(f"No reference genome FASTA files found in {reference_path}")
            sys.exit(1)
            return None

        # Use the first reference file found
        reference_fasta = os.path.join(reference_path, ref_files[0])

        sample_path = None

        # If a specific FASTA path is provided, use it
        if query_fasta_path:
            if os.path.exists(query_fasta_path) and (
                query_fasta_path.endswith(".fasta") or query_fasta_path.endswith(".fa")
            ):
                sample_path = query_fasta_path
                logger.info(f"Using provided query FASTA file: {sample_path}")
            else:
                logger.error(
                    f"Provided file does not exist or is not a FASTA file: {query_fasta_path}"
                )
                return None
        # Otherwise, try to find a file in the uploads directory
        else:
            user_uploads_folder = Path(dir["uploads_dir"])
            logger.info(f"Checking uploads folder: {user_uploads_folder}")

            if user_uploads_folder.exists():
                sample_files = [
                    f
                    for f in os.listdir(user_uploads_folder)
                    if f.endswith(".fasta") or f.endswith(".fa")
                ]

                if sample_files:
                    sample_path = os.path.join(user_uploads_folder, sample_files[0])
                    logger.info(f"Using query file from uploads: {sample_path}")
                else:
                    logger.error(
                        f"No FASTA files found in uploads directory: {user_uploads_folder}"
                    )
                    return None
            else:
                logger.error(f"Sample directory does not exist: {user_uploads_folder}")
                return None

        # Run BLAST using command line
        result_file = blast_cmdline(
            str(sample_path), str(reference_fasta), str(blast_dir)
        )
        logger.info(f"BLAST analysis complete. Results saved to {result_file}")
        print(f"BLAST analysis complete. Results saved to {result_file}")
        return result_file

    except subprocess.CalledProcessError:
        logger.error(
            "BLAST+ tools not found. Please install BLAST+ (sudo apt install ncbi-blast+)"
        )
    except Exception as e:
        logger.exception(f"Error running BLAST: {e}")
        return None


# Example usage
if __name__ == "__main__":
    # Example with specific path
    # perform_blast_aligning("/path/to/your/specific/file.fasta")
    disease_path = "C:/Users/Kajus/.kath/shared/data/fasta_files/samples/diseases/alzheimer/unknown_transcript1_1519314819.fasta"

    # Example with default behavior (using uploads directory)
    perform_blast_aligning(disease_path)
