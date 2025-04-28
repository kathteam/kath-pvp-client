"""
FASTA Downloader Module for bioinformatics data retrieval.

This module provides functions to download FASTA files from various bioinformatics
databases including NCBI, Ensembl, and ClinVar. It also includes functions
for running BLAST alignments and parsing results to identify sequence variations.
"""

import os
import re
import subprocess
import json
from pathlib import Path
from typing import Optional, List, Dict, Any
from Bio.Blast import NCBIXML

from src.shared.constants import (
    PROGRAM_STORAGE_DIR_SHARED_BLAST,
    PROGRAM_STORAGE_DIR_SHARED_DATA_FASTA,
)

PROGRAM_STORAGE_DIR_SHARED_BLAST.mkdir(parents=True, exist_ok=True)
Path(PROGRAM_STORAGE_DIR_SHARED_DATA_FASTA, "reference").mkdir(
    parents=True, exist_ok=True
)
Path(PROGRAM_STORAGE_DIR_SHARED_BLAST, "uploads").mkdir(parents=True, exist_ok=True)

import logging

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def blast_cmdline(
    query_fasta_path: str, reference_genome_path: str, output_dir: Path
) -> str:
    """
    Run BLAST using direct command-line execution of BLAST+ tools.

    Args:
        query_fasta_path: Path to query FASTA file.
        reference_genome_path: Path to reference genome FASTA file.
        output_dir: Directory to save BLAST results.

    Returns:
        Path to BLAST results file (XML format).

    Raises:
        RuntimeError: If makeblastdb or blastn commands fail.
        FileNotFoundError: If input files do not exist.
    """
    query_path = Path(query_fasta_path)
    ref_path = Path(reference_genome_path)
    output_dir = Path(output_dir)

    if not query_path.is_file():
        raise FileNotFoundError(f"Query FASTA file not found: {query_fasta_path}")
    if not ref_path.is_file():
        raise FileNotFoundError(
            f"Reference genome file not found: {reference_genome_path}"
        )

    output_dir.mkdir(parents=True, exist_ok=True)

    query_filename = query_path.stem
    reference_filename = ref_path.stem
    output_file = output_dir / f"{query_filename}_vs_{reference_filename}_blast.xml"
    db_path = ref_path.parent / reference_filename

    logger.info(
        f"Running BLAST+ alignment for {query_filename} against {reference_filename}..."
    )
    db_files = list(ref_path.parent.glob(f"{reference_filename}.n*"))
    if not db_files:
        logger.info(f"Creating BLAST database for {reference_filename}...")
        makeblastdb_cmd = [
            "makeblastdb",
            "-in",
            str(ref_path),
            "-dbtype",
            "nucl",
            "-parse_seqids",  # Important for easier parsing of standard identifiers
            "-out",
            str(db_path),
            "-title",
            reference_filename,  # Add a title to the database
        ]
        try:
            process = subprocess.run(
                makeblastdb_cmd,
                capture_output=True,
                text=True,
                check=True,
                encoding="utf-8",
            )
            logger.info("BLAST database created successfully.")
            logger.debug(f"makeblastdb stdout:\n{process.stdout}")
            logger.debug(f"makeblastdb stderr:\n{process.stderr}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to create BLAST database: {e.stderr}")
            raise RuntimeError(
                f"makeblastdb failed with exit code {e.returncode}"
            ) from e
        except FileNotFoundError:
            logger.error(
                "`makeblastdb` command not found. Is BLAST+ installed and in PATH?"
            )
            raise
    else:
        logger.info(f"Using existing BLAST database at {db_path}")

    blastn_cmd = [
        "blastn",
        "-query",
        str(query_path),
        "-db",
        str(db_path),
        "-out",
        str(output_file),
        "-outfmt",
        "5",  # XML output format
        "-evalue",
        "1e-10",  # Keep stringent E-value for high confidence hits
        "-word_size",
        "28",  # Good for highly similar sequences (like human vs human)
        "-perc_identity",
        "90",  # Keep high, but consider lowering slightly (e.g., 90) if needed
        "-num_threads",
        str(os.cpu_count() or 1),  # Use available cores
        # "-task", "megablast", # Often default/preferred for intra-species, slightly faster
        "-task",
        "blastn",  # Using user's original choice, also very suitable
        # "-gapopen", "5", "-gapextend", "2", # Example non-default penalties
    ]

    logger.info(f"Executing blastn command: {' '.join(blastn_cmd)}")
    try:
        process = subprocess.run(
            blastn_cmd, capture_output=True, text=True, check=True, encoding="utf-8"
        )
        logger.info("BLAST alignment completed successfully.")
        logger.debug(f"blastn stdout:\n{process.stdout}")
        logger.debug(f"blastn stderr:\n{process.stderr}")
    except subprocess.CalledProcessError as e:
        logger.error(f"BLAST alignment failed: {e.stderr}")
        # Log stdout as well, might contain info even on failure
        logger.error(f"BLAST stdout: {e.stdout}")
        raise RuntimeError(f"blastn failed with exit code {e.returncode}") from e
    except FileNotFoundError:
        logger.error("`blastn` command not found. Is BLAST+ installed and in PATH?")
        raise

    return str(output_file)


def extract_chromosome(subject_id: str) -> str:
    """
    Extract chromosome information from the subject ID (Hit_id or Hit_def).
    Handles common formats like 'chr1', 'NC_000001.11', '1'.

    Args:
        subject_id: The subject ID string from BLAST results.

    Returns:
        Chromosome identifier (e.g., "1", "X", "MT") or "Unknown".
    """
    chr_match = re.search(r"chr([0-9]{1,2}|[XYMxym])", subject_id, re.IGNORECASE)
    if chr_match:
        chrom = chr_match.group(1).upper()
        return "MT" if chrom == "M" else chrom

    # 2. Check for RefSeq chromosome format (NC_0000XX.version)
    nc_match = re.search(r"NC_0000(\d{2})\.\d+", subject_id)
    if nc_match:
        num = int(nc_match.group(1))
        if 1 <= num <= 22:
            return str(num)
        if num == 23:
            return "X"
        if num == 24:
            return "Y"
        if subject_id.startswith("NC_012920"):
            return "MT"
        return "Unknown_NC"

    if re.search(r"^(NT_|NW_|NZ_)", subject_id):
        return "Contig/Scaffold"

    num_match = re.search(r"\b([0-9]{1,2}|[XYMxy])\b", subject_id)
    if num_match:
        chrom = num_match.group(1).upper()
        return "MT" if chrom == "M" else chrom

    ucsc_match = re.search(r"(chr[0-9]{1,2}|[XYMxym])[_.]", subject_id, re.IGNORECASE)
    if ucsc_match:
        chrom = ucsc_match.group(1).upper().replace("CHR", "")
        return "MT" if chrom == "M" else chrom

    logger.warning(f"Could not determine chromosome from subject ID: {subject_id}")
    return "Unknown"


def parse_blast_results(
    result_file: str,
    min_hsp_score: float = 50.0,
    min_identity_perc: float = 95.0,
    min_hsp_length: int = 50,
    evalue_threshold: float = 1e-10,
) -> List[Dict[str, Any]]:
    """
    Parse BLAST XML results, identify variations (substitutions, insertions, deletions),
    and filter based on alignment quality metrics.

    Args:
        result_file: Path to BLAST XML output file.
        min_hsp_score: Minimum bit score for an HSP to be considered.
        min_identity_perc: Minimum percent identity for an HSP (0-100).
        min_hsp_length: Minimum alignment length for an HSP.
        evalue_threshold: Maximum E-value for an HSP.

    Returns:
        A list of dictionaries, each describing a detected variation.
    """
    variations = []
    logger.info(f"Parsing BLAST results from: {result_file}")

    try:
        with open(result_file, "r", encoding="utf-8") as result_handle:
            blast_records = NCBIXML.parse(result_handle)

            query_count = 0
            total_hsps_processed = 0
            total_variations_found = 0

            for blast_record in blast_records:
                query_count += 1
                logger.debug(
                    f"Processing Query: {blast_record.query} ({blast_record.query_length} bp)"
                )
                if not blast_record.alignments:
                    logger.debug(
                        f"  No significant alignments found for query {blast_record.query}"
                    )
                    continue

                for alignment in blast_record.alignments:
                    subject_id = alignment.hit_id or alignment.title
                    chromosome = extract_chromosome(subject_id)
                    logger.debug(
                        f"  Alignment Hit: {subject_id} (Length: {alignment.length} bp) -> Chromosome: {chromosome}"
                    )

                    hsp_index = 0
                    for hsp in alignment.hsps:
                        hsp_index += 1
                        total_hsps_processed += 1
                        percent_identity = (hsp.identities / hsp.align_length) * 100

                        # --- Quality Filtering ---
                        if (
                            hsp.expect > evalue_threshold
                            or hsp.bits < min_hsp_score
                            or percent_identity < min_identity_perc
                            or hsp.align_length < min_hsp_length
                        ):
                            logger.debug(
                                f"    Skipping HSP {hsp_index} (Score: {hsp.bits:.1f}, "
                                f"E-value: {hsp.expect:.2g}, Ident: {percent_identity:.1f}%, "
                                f"Len: {hsp.align_length}) due to quality filters."
                            )
                            continue

                        logger.debug(
                            f"    Processing HSP {hsp_index} (Score: {hsp.bits:.1f}, "
                            f"E-value: {hsp.expect:.2g}, Ident: {percent_identity:.1f}%, "
                            f"Len: {hsp.align_length}, Strand: {hsp.strand})"
                        )

                        # --- Detailed Variation Extraction ---
                        query_seq = hsp.query
                        match_seq = hsp.match
                        sbjct_seq = hsp.sbjct

                        query_pos = hsp.query_start - 1  # 0-based index
                        sbjct_pos = hsp.sbjct_start - 1  # 0-based index

                        hsp_variations = 0
                        for i in range(hsp.align_length):
                            q_base = query_seq[i]
                            s_base = sbjct_seq[i]
                            m_char = match_seq[i]

                            q_pos_current = -1
                            if q_base != "-":
                                query_pos += 1
                                q_pos_current = query_pos

                            s_pos_current = -1
                            if s_base != "-":
                                sbjct_pos += 1
                                s_pos_current = sbjct_pos

                            # Identify variation type
                            variation_type = None
                            ref_allele = None
                            alt_allele = None

                            if q_base != s_base:
                                if q_base == "-":
                                    variation_type = "insertion"
                                    ref_allele = "-"
                                    alt_allele = s_base
                                elif s_base == "-":
                                    variation_type = "deletion"
                                    ref_allele = q_base
                                    alt_allele = "-"
                                else:  # Substitution (Mismatch)
                                    variation_type = "substitution"
                                    ref_allele = s_base  # Reference base from subject
                                    alt_allele = q_base  # Alternative base from query

                            if variation_type:
                                hsp_variations += 1
                                variation_details = {
                                    "query_id": blast_record.query,
                                    "subject_id": subject_id,
                                    "chromosome": chromosome,
                                    "position": (
                                        s_pos_current + 1
                                        if s_pos_current != -1
                                        else None
                                    ),
                                    "variation_type": variation_type,
                                    "reference_allele": ref_allele,
                                    "query_allele": alt_allele,
                                    "query_position": (
                                        q_pos_current + 1
                                        if q_pos_current != -1
                                        else None
                                    ),
                                    "hsp_score": hsp.bits,
                                    "hsp_evalue": hsp.expect,
                                    "hsp_identity": percent_identity,
                                    "hsp_align_length": hsp.align_length,
                                    "hsp_query_start": hsp.query_start,
                                    "hsp_subject_start": hsp.sbjct_start,
                                    "hsp_strand": hsp.strand,
                                    "hsp_gaps": hsp.gaps,
                                }
                                variations.append(variation_details)

                        if hsp_variations > 0:
                            logger.debug(
                                f"      Found {hsp_variations} variations in this HSP."
                            )
                        total_variations_found += hsp_variations

        logger.info(f"Processed {query_count} queries and {total_hsps_processed} HSPs.")
        if variations:
            logger.info(
                f"Total variations found meeting quality criteria: {len(variations)}"
            )
        else:
            logger.info("No variations found meeting quality criteria after filtering.")

    except FileNotFoundError:
        logger.error(f"BLAST result file not found: {result_file}")
        return []
    except Exception as e:
        logger.exception(f"Error parsing BLAST XML file {result_file}: {e}")
        return []

    return variations


# ... (keep existing imports and constants) ...
# ... (keep blast_cmdline, extract_chromosome, parse_blast_results functions) ...

def _run_blast_and_parse_single(
    actual_query_path: Path,
    reference_fasta: Path,
    blast_output_dir: Path,
) -> Optional[str]:
    """
    Helper function to run BLAST, parse results, and save JSON for a single query file.

    Args:
        actual_query_path: Path to the query FASTA file.
        reference_fasta: Path to the reference FASTA file.
        blast_output_dir: Directory to store BLAST XML and JSON results.

    Returns:
        Path to the generated JSON file, or None if an error occurred.
    """
    logger.info(f"--- Processing Query File: {actual_query_path.name} ---")
    query_stem = actual_query_path.stem
    reference_stem = reference_fasta.stem
    # Define unique output filenames based on query and reference
    results_json_filename = f"{query_stem}_vs_{reference_stem}_variations.json"
    json_output_path = blast_output_dir / results_json_filename

    # --- Run BLAST ---
    try:
        blast_result_xml = blast_cmdline(
            query_fasta_path=str(actual_query_path),
            reference_genome_path=str(reference_fasta),
            output_dir=blast_output_dir, # blast_cmdline already creates unique XML name
        )

    except (RuntimeError, FileNotFoundError) as e:
        logger.error(
            f"BLAST command execution failed for {actual_query_path.name}: {e}"
        )
        return None # Indicate failure for this file
    except Exception as e:
        logger.exception(
            f"An unexpected error occurred during BLAST execution for {actual_query_path.name}: {e}"
        )
        return None # Indicate failure for this file

    # --- Parse Results ---
    try:
        variations = parse_blast_results(blast_result_xml)

        # --- Save Variations to JSON ---
        # Always save a JSON file, even if empty, to indicate processing occurred
        with open(json_output_path, "w", encoding="utf-8") as f:
            json.dump(variations if variations else [], f, indent=2)

        if variations:
            logger.info(
                f"Saved {len(variations)} identified variations to: {json_output_path}"
            )
        else:
            logger.warning(
                f"No variations met the filtering criteria for {actual_query_path.name}."
            )
            print(
                f"Analysis complete for {actual_query_path.name}. No significant variations found. Results file: {json_output_path}"
            )
        return str(json_output_path) # Return path whether variations were found or not

    except Exception as e:
        logger.exception(
            f"An unexpected error occurred during BLAST parsing or saving results for {actual_query_path.name}: {e}"
        )
        return None # Indicate failure for this file


def process_single_fasta(
    query_fasta_path: str,
    reference_fasta_dir: Path = Path(PROGRAM_STORAGE_DIR_SHARED_DATA_FASTA) / "reference",
    blast_output_dir: Path = Path(PROGRAM_STORAGE_DIR_SHARED_BLAST),
) -> Optional[str]:
    """
    Orchestrates BLAST alignment and variation parsing for a single specified query FASTA file.

    Args:
        query_fasta_path: Path to the specific query FASTA file.
        reference_fasta_dir: Directory containing the reference genome FASTA file(s).
        blast_output_dir: Directory to store BLAST XML output and the final JSON results.

    Returns:
        Path to the JSON file containing identified variations, or None if an error occurs
        during setup or processing.
    """
    logger.info(f"Starting single file BLAST analysis for: {query_fasta_path}")

    # --- Sanity Checks and Setup ---
    try:
        subprocess.run(["blastn", "-version"], capture_output=True, text=True, check=True)
        subprocess.run(["makeblastdb", "-version"], capture_output=True, text=True, check=True)
        logger.info("BLAST+ tools found.")
    except (subprocess.CalledProcessError, FileNotFoundError):
        logger.error("BLAST+ commands not found in system PATH.")
        return None

    if not reference_fasta_dir.is_dir():
        logger.error(f"Reference FASTA directory not found: {reference_fasta_dir}")
        return None

    ref_files = list(reference_fasta_dir.glob("*.fasta")) + list(reference_fasta_dir.glob("*.fa"))
    if not ref_files:
        logger.error(f"No reference genome FASTA files found in {reference_fasta_dir}")
        return None
    reference_fasta = ref_files[0]
    logger.info(f"Using reference genome: {reference_fasta}")

    # --- Validate Query File ---
    query_path_obj = Path(query_fasta_path)
    if not query_path_obj.is_file() or query_path_obj.suffix.lower() not in [".fasta", ".fa", ".fna"]:
        logger.error(f"Provided query file is not a valid FASTA file or does not exist: {query_fasta_path}")
        return None

    # Ensure output directory exists
    blast_output_dir.mkdir(parents=True, exist_ok=True)

    # --- Process the Single File ---
    result_path = _run_blast_and_parse_single(
        actual_query_path=query_path_obj,
        reference_fasta=reference_fasta,
        blast_output_dir=blast_output_dir,
    )

    if result_path:
        logger.info(f"Single file workflow finished successfully for {query_fasta_path}.")
    else:
        logger.error(f"Single file workflow failed for {query_fasta_path}.")

    return result_path


def process_fasta_directory(
    uploads_dir: Path = Path(PROGRAM_STORAGE_DIR_SHARED_BLAST) / "uploads",
    reference_fasta_dir: Path = Path(PROGRAM_STORAGE_DIR_SHARED_DATA_FASTA) / "reference",
    blast_output_dir: Path = Path(PROGRAM_STORAGE_DIR_SHARED_BLAST),
) -> List[str]:
    """
    Orchestrates BLAST alignment and variation parsing for all FASTA files in a directory.

    Args:
        uploads_dir: Directory to search for query FASTA files (.fa, .fasta, .fna).
        reference_fasta_dir: Directory containing the reference genome FASTA file(s).
        blast_output_dir: Directory to store BLAST XML output and the final JSON results.

    Returns:
        A list of paths to the JSON files containing identified variations for each
        successfully processed query file. Returns an empty list if no files were
        processed or critical setup errors occurred.
    """
    logger.info(f"Starting directory BLAST analysis for: {uploads_dir}")
    processed_json_files = []

    # --- Sanity Checks and Setup ---
    try:
        subprocess.run(["blastn", "-version"], capture_output=True, text=True, check=True)
        subprocess.run(["makeblastdb", "-version"], capture_output=True, text=True, check=True)
        logger.info("BLAST+ tools found.")
    except (subprocess.CalledProcessError, FileNotFoundError):
        logger.error("BLAST+ commands not found in system PATH.")
        return []

    if not reference_fasta_dir.is_dir():
        logger.error(f"Reference FASTA directory not found: {reference_fasta_dir}")
        return []

    ref_files = list(reference_fasta_dir.glob("*.fasta")) + list(reference_fasta_dir.glob("*.fa"))
    if not ref_files:
        logger.error(f"No reference genome FASTA files found in {reference_fasta_dir}")
        return []
    reference_fasta = ref_files[0]
    logger.info(f"Using reference genome: {reference_fasta}")

    # --- Find Query Files ---
    if not uploads_dir.is_dir():
        logger.error(f"Uploads directory does not exist: {uploads_dir}")
        return []

    query_files_to_process = (
        list(uploads_dir.glob("*.fasta"))
        + list(uploads_dir.glob("*.fa"))
        + list(uploads_dir.glob("*.fna"))
    )
    if not query_files_to_process:
        logger.error(f"No FASTA files (.fa, .fasta, .fna) found in uploads directory: {uploads_dir}")
        return []
    logger.info(f"Found {len(query_files_to_process)} query files in uploads directory.")

    # Ensure output directory exists
    blast_output_dir.mkdir(parents=True, exist_ok=True)

    # --- Process Each Query File ---
    for actual_query_path in query_files_to_process:
        result_path = _run_blast_and_parse_single(
            actual_query_path=actual_query_path,
            reference_fasta=reference_fasta,
            blast_output_dir=blast_output_dir,
        )
        if result_path:
            processed_json_files.append(result_path)
        # If result_path is None, the error was already logged in the helper function

    logger.info(f"Directory workflow finished. Successfully processed {len(processed_json_files)} out of {len(query_files_to_process)} files.")
    return processed_json_files


# Example usage:
if __name__ == "__main__":
    # dummy_query_path = Path(PROGRAM_STORAGE_DIR_SHARED_BLAST) / "uploads" / "single_test.fasta"
    # # if not dummy_query_path.exists():
    # #     dummy_query_path.parent.mkdir(parents=True, exist_ok=True)
    # #     with open(dummy_query_path, "w") as f:
    # #         f.write(">TestQuery1\nACGTACGT\n") # Replace with actual sequence if needed
    # #     print(f"Created dummy query file: {dummy_query_path}")

    # print(f"\n--- Running analysis for a single specified file: {dummy_query_path} ---")
    # # Make sure the dummy file exists and reference exists before running
    # if dummy_query_path.exists():
    #      single_result = process_single_fasta(query_fasta_path=str(dummy_query_path))
    #      if single_result:
    #          print(f"Single file processing complete. Result JSON: {single_result}")
    #      else:
    #          print(f"Single file processing failed for {dummy_query_path}.")
    # else:
    #      print(f"Skipping single file test, query file not found: {dummy_query_path}")


    # --- Example for processing a directory ---
    print("\n--- Running analysis using file(s) from uploads directory ---")
    # Ensure upload dir exists and contains FASTA files
    # Ensure reference genome exists in ./fasta_data/reference/

    result_json_list = process_fasta_directory() # Uses default uploads_dir
    if result_json_list:
        print(f"Directory processing complete. Result JSON files generated:")
        print(f"Directory processing complete. Result JSON files generated:")
        for json_path in result_json_list:
            print(f"- {json_path}")
            from find import process_variants
            process_variants(json_path) # Process each JSON file for variants

    else:
            "Directory analysis failed, found no query files, or encountered errors during processing."

# Remove the old perform_blast_alignment_and_find_variations function
# del perform_blast_alignment_and_find_variations