import os
import time
import logging
import re
import csv
from pathlib import Path
from typing import Optional, List, Dict, Tuple, Union

# Third-party imports
from dotenv import load_dotenv
from Bio.Blast import NCBIXML

# Set up logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
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

    # Create directories if they don"t exist
    for dir_path in [
        service_dir,
        service_dir,
        fasta_dir,
        blast_dir,
        uploads_dir,
        ref_dir,
        samples_dir,
    ]:
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
        "samples_dir": samples_dir,
    }


# Default timeout for HTTP requests (in seconds)
DEFAULT_TIMEOUT = 60


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

    # Try to find chromosome in format "chr1" or "chrX"
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

    # Try to extract just a number if it"s at the beginning or isolated
    num_match = re.search(r"\b([0-9]+|[XYxy])\b", subject_id)
    if num_match:
        return num_match.group(1).upper()

    return "Unknown"


def parse_blast_results(
    result_file: Union[str, Path], return_mutations: bool = False
) -> Optional[Tuple[List[Dict], List[Dict]]]:
    """
    Parse and log BLAST results from XML file.

    Args:
        result_file: Path to the BLAST XML results file
        return_mutations: Whether to return mutation data

    Returns:
        Tuple of (mutations, alignment_summary) if return_mutations is True, otherwise None
    """
    result_file = Path(result_file)
    logger.info(f"Parsing BLAST results from {result_file}")

    if not result_file.exists():
        logger.error(f"Result file not found: {result_file}")
        return ([], []) if return_mutations else None

    all_mutations = []
    alignment_summary = []

    try:
        with open(result_file) as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            for i, blast_record in enumerate(blast_records):
                query_id = blast_record.query.split()[0]
                query_length = blast_record.query_length

                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        # Extract subject/reference ID
                        subject_id = (
                            alignment.title.split()[1]
                            if len(alignment.title.split()) > 1
                            else alignment.title
                        )

                        # Extract chromosome information
                        chromosome = extract_chromosome(subject_id)

                        # Calculate alignment statistics
                        identity_pct = hsp.identities / hsp.align_length * 100
                        coverage_pct = (hsp.align_length / query_length) * 100

                        alignment_info = {
                            "query_id": query_id,
                            "subject_id": subject_id,
                            "chromosome": chromosome,
                            "identity_pct": round(identity_pct, 2),
                            "coverage_pct": round(coverage_pct, 2),
                            "align_length": hsp.align_length,
                            "e_value": hsp.expect,
                            "score": hsp.score,
                            "query_start": hsp.query_start,
                            "query_end": hsp.query_end,
                            "subject_start": hsp.sbjct_start,
                            "subject_end": hsp.sbjct_end,
                        }
                        alignment_summary.append(alignment_info)

                        # Extract mutations
                        mutations = extract_mutations_from_hsp(
                            hsp, chromosome, query_id, subject_id
                        )
                        all_mutations.extend(mutations)

                        if mutations:
                            log_mutation_summary(mutations)

        # Always save mutations to CSV with proper CSV format
        csv_file = result_file.with_suffix(".csv")
        save_mutations_to_csv(all_mutations, csv_file)

        if return_mutations:
            return all_mutations, alignment_summary

        return None

    except Exception as e:
        logger.error(f"Error parsing BLAST results: {e}")
        # Create empty CSV with proper structure on error
        try:
            csv_file = result_file.with_suffix(".csv")
            with open(csv_file, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["error"])
                writer.writerow([f"Error parsing BLAST results: {e}"])
        except Exception:
            pass

        if return_mutations:
            return [], []
        return None


def extract_mutations_from_hsp(hsp, chromosome: str, query_id: str, subject_id: str) -> List[Dict]:
    """
    Extract mutations from a high-scoring segment pair.

    Args:
        hsp: BioPython HSP object
        chromosome: Chromosome identifier
        query_id: Query sequence ID
        subject_id: Subject sequence ID

    Returns:
        List of mutation dictionaries
    """
    mutations = []
    position_offset = 0  # Track position changes due to gaps

    for i, (q, s) in enumerate(zip(hsp.query, hsp.sbjct)):
        ref_position = hsp.sbjct_start + i - position_offset

        # Handle mismatches (SNPs)
        if q != s and q != "-" and s != "-":
            mutation = {
                "type": "SNP",
                "chromosome": chromosome,
                "position": ref_position,
                "reference": s,
                "query": q,
                "context": hsp.sbjct[max(0, i - 5) : i]
                + "["
                + s
                + "]"
                + hsp.sbjct[i + 1 : min(i + 6, len(hsp.sbjct))],
                "query_id": query_id,
                "subject_id": subject_id,
                "alignment_start": hsp.sbjct_start,
                "alignment_end": hsp.sbjct_end,
                "e_value": hsp.expect,
            }
            mutations.append(mutation)

        # Handle deletions in query (gaps in query)
        elif q == "-" and s != "-":
            position_offset += 1
            mutation = {
                "type": "DEL",
                "chromosome": chromosome,
                "position": ref_position,
                "reference": s,
                "query": "-",
                "context": hsp.sbjct[max(0, i - 5) : i]
                + "["
                + s
                + "]"
                + hsp.sbjct[i + 1 : min(i + 6, len(hsp.sbjct))],
                "query_id": query_id,
                "subject_id": subject_id,
                "alignment_start": hsp.sbjct_start,
                "alignment_end": hsp.sbjct_end,
                "e_value": hsp.expect,
            }
            mutations.append(mutation)

        # Handle insertions in query (gaps in reference)
        elif q != "-" and s == "-":
            mutation = {
                "type": "INS",
                "chromosome": chromosome,
                "position": ref_position,
                "reference": "-",
                "query": q,
                "context": hsp.sbjct[max(0, i - 5) : i]
                + "[-]"
                + hsp.sbjct[i + 1 : min(i + 6, len(hsp.sbjct))],
                "query_id": query_id,
                "subject_id": subject_id,
                "alignment_start": hsp.sbjct_start,
                "alignment_end": hsp.sbjct_end,
                "e_value": hsp.expect,
            }
            mutations.append(mutation)

    return mutations


def log_mutation_summary(mutations: List[Dict]) -> None:
    """
    Log a summary of the mutations found.

    Args:
        mutations: List of mutation dictionaries
    """
    logger.info(f"  Found {len(mutations)} mutations")
    mutation_types = {}
    for mut in mutations:
        mutation_types[mut["type"]] = mutation_types.get(mut["type"], 0) + 1

    for mut_type, count in mutation_types.items():
        logger.info(f"    {mut_type}: {count}")

    # Display sample mutations
    for mut in mutations[:5]:
        logger.info(
            f"    {mut['type']} at Chr {mut['chromosome']} position {mut['position']}: "
            f"{mut['reference']} -> {mut['query']} (context: {mut['context']})"
        )

    if len(mutations) > 5:
        logger.info(f"    ... and {len(mutations) - 5} more mutations")


def save_mutations_to_csv(mutations: List[Dict], csv_file: Union[str, Path]) -> None:
    """
    Save mutations to a CSV file.

    Args:
        mutations: List of mutation dictionaries
        csv_file: Path to output CSV file
    """
    try:
        with open(csv_file, "w", newline="") as f:
            if not mutations:
                writer = csv.writer(f)
                writer.writerow(["status"])
                writer.writerow(["No mutations found"])
                return

            writer = csv.DictWriter(f, fieldnames=mutations[0].keys())
            writer.writeheader()
            writer.writerows(mutations)
        logger.info(f"Mutations data saved to {csv_file}")
    except Exception as e:
        logger.error(f"Failed to save mutations to CSV: {e}")


def analyze_blast_xml():
    """
    Main function that processes all BLAST results and outputs CSV-formatted mutation data
    into a single consolidated file.
    """
    # Initialize directories
    dir = initialize_directories()
    BLAST_DIR = dir["blast_dir"]

    # Create output directory for consolidated results
    consolidated_dir = os.path.join(BLAST_DIR, "consolidated")
    os.makedirs(consolidated_dir, exist_ok=True)

    # Get current timestamp for unique filename
    timestamp = time.strftime("%Y%m%d-%H%M%S")
    consolidated_csv = os.path.join(consolidated_dir, f"all_mutations_{timestamp}.csv")

    files = os.listdir(BLAST_DIR)

    if not files:
        logger.error("No BLAST result files found in directory")
        # Create empty CSV indicating no files found
        with open(consolidated_csv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["error"])
            writer.writerow(["No BLAST result files found in directory"])
        return consolidated_csv

    # Process all XML files found
    xml_files = [f for f in files if f.endswith(".xml")]
    if not xml_files:
        logger.error("No XML BLAST result files found in directory")
        with open(consolidated_csv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["error"])
            writer.writerow(["No XML BLAST result files found in directory"])
        return consolidated_csv

    # Collect all mutations
    all_mutations = []
    all_alignments = []

    for xml_file in xml_files:
        file_path = os.path.join(BLAST_DIR, xml_file)
        file_id = os.path.splitext(xml_file)[0]  # Use filename without extension as ID

        try:
            logger.info(f"Processing {xml_file}...")
            mutations, alignments = parse_blast_results(file_path, return_mutations=True)

            # Add file identifier to each mutation
            for mut in mutations:
                mut["file_id"] = file_id

            # Add file identifier to each alignment
            for align in alignments:
                align["file_id"] = file_id

            all_mutations.extend(mutations)
            all_alignments.extend(alignments)

        except Exception as e:
            logger.error(f"Error processing {xml_file}: {e}")
            # Add an error entry so we know this file had problems
            all_mutations.append(
                {
                    "file_id": file_id,
                    "type": "ERROR",
                    "error": str(e),
                    "chromosome": "N/A",
                    "position": 0,
                    "reference": "N/A",
                    "query": "N/A",
                    "context": "N/A",
                    "query_id": "N/A",
                    "subject_id": "N/A",
                    "alignment_start": 0,
                    "alignment_end": 0,
                    "e_value": 0.0,
                }
            )

    # Save consolidated mutations to CSV
    try:
        with open(consolidated_csv, "w", newline="") as f:
            if not all_mutations:
                writer = csv.writer(f)
                writer.writerow(["status"])
                writer.writerow(["No mutations found in any files"])
            else:
                # Ensure all mutations have the same fields (add missing ones)
                all_fields = set()
                for mut in all_mutations:
                    all_fields.update(mut.keys())

                # Make sure "file_id" is the first field
                fieldnames = ["file_id"]
                for field in all_fields:
                    if field != "file_id":
                        fieldnames.append(field)

                # Fill in missing fields with None
                for mut in all_mutations:
                    for field in all_fields:
                        if field not in mut:
                            mut[field] = None

                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(all_mutations)

        logger.info(f"All mutations saved to {consolidated_csv}")

        # Also save consolidated alignments if there are any
        if all_alignments:
            alignment_csv = os.path.join(consolidated_dir, f"all_alignments_{timestamp}.csv")
            with open(alignment_csv, "w", newline="") as f:
                if not all_alignments:
                    writer = csv.writer(f)
                    writer.writerow(["status"])
                    writer.writerow(["No alignments found"])
                else:
                    # Same as for mutations, normalize fields
                    all_fields = set()
                    for align in all_alignments:
                        all_fields.update(align.keys())

                    fieldnames = ["file_id"]
                    for field in all_fields:
                        if field != "file_id":
                            fieldnames.append(field)

                    for align in all_alignments:
                        for field in all_fields:
                            if field not in align:
                                align[field] = None

                    writer = csv.DictWriter(f, fieldnames=fieldnames)
                    writer.writeheader()
                    writer.writerows(all_alignments)

            logger.info(f"All alignments saved to {alignment_csv}")

        # Create a simple summary CSV with counts per file
        summary_csv = os.path.join(consolidated_dir, f"summary_{timestamp}.csv")
        with open(summary_csv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(
                ["file_id", "total_mutations", "SNP_count", "INS_count", "DEL_count", "ERROR_count"]
            )

            # Group mutations by file_id
            file_mutations = {}
            for mut in all_mutations:
                file_id = mut["file_id"]
                if file_id not in file_mutations:
                    file_mutations[file_id] = []
                file_mutations[file_id].append(mut)

            # Write summary for each file
            for file_id, mutations in file_mutations.items():
                type_counts = {"SNP": 0, "INS": 0, "DEL": 0, "ERROR": 0}
                for mut in mutations:
                    mut_type = mut.get("type", "Unknown")
                    if mut_type in type_counts:
                        type_counts[mut_type] += 1

                writer.writerow(
                    [
                        file_id,
                        len(mutations),
                        type_counts["SNP"],
                        type_counts["INS"],
                        type_counts["DEL"],
                        type_counts["ERROR"],
                    ]
                )

        logger.info(f"Summary saved to {summary_csv}")

        return consolidated_csv

    except Exception as e:
        logger.error(f"Error saving consolidated results: {e}")
        error_csv = os.path.join(consolidated_dir, f"error_{timestamp}.csv")
        with open(error_csv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["error"])
            writer.writerow([f"Error saving consolidated results: {str(e)}"])
        return error_csv


if __name__ == "__main__":
    analyze_blast_xml()
