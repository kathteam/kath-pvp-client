import os
import sys
backend_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../"))
if backend_dir not in sys.path:
    sys.path.append(backend_dir)
import time
import gzip
import shutil
from pathlib import Path
from typing import Optional, List, Tuple
import requests
from Bio import Entrez
from dotenv import load_dotenv

from src.utils.logger import get_logger
from services.utils.script_setup import (
    EnvSetup,
    FolderSetup,
)

EnvSetup()
FolderSetup()

from shared.constants import REF_GENOME, PROGRAM_STORAGE_DIR_SHARED_DATA_FASTA_SAMPLES, PROGRAM_STORAGE_DIR_SHARED_DATA_FASTA_UPLOADS, PROGRAM_STORAGE_DIR_SHARED_DATA_FASTA

logger = get_logger(__name__)

# Default timeout for HTTP requests (in seconds)
DEFAULT_TIMEOUT = 60

def log_download(source: str, identifier: str, file_path: str, fasta_dir: Path):
    """Log a successful download."""
    log_file = os.path.join(fasta_dir, "downloads.log")
    with open(log_file, "a", encoding="utf-8") as f:
        f.write(f"{time.strftime("%Y-%m-%d %H:%M:%S")} | {source} | {identifier} | {file_path}\n")


def reference_genome_exists(ref_dir: Path, version: str = REF_GENOME) -> bool:
    """Check if the reference genome file exists."""
    return any(Path(ref_dir).glob(f"{version}_direct.fasta"))


def download_reference_genome_direct(
    version: str = REF_GENOME, output_dir: Optional[Path] = None
) -> str:
    """Download reference genome directly from NCBI FTP."""
    if output_dir is None:
        output_dir = Path(
            os.path.join(
                os.path.expanduser("~"), ".kath", "shared", "data", "fasta_files", "reference"
            )
        )

    output_file = os.path.join(output_dir, f"{version}_direct.fasta")

    if os.path.exists(output_file):
        logger.info(f"Reference genome {version} already exists at {output_file}")
        return str(output_file)

    # URLs for complete reference genomes
    # pylint: disable=line-too-long
    urls = {
        "GRCh38": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz",
        "GRCh37": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.14_GRCh37.p13/GCA_000001405.14_GRCh37.p13_genomic.fna.gz",
    }
    # pylint: enable=line-too-long

    url = urls.get(version)
    if not url:
        raise ValueError(f"Unsupported genome version: {version}")

    logger.info(f"Downloading and decompressing {version} from {url}...")

    try:
        # Download and decompress in a single operation directly to the output file
        response = requests.get(url, stream=True, timeout=DEFAULT_TIMEOUT * 2)
        response.raise_for_status()

        with gzip.GzipFile(fileobj=response.raw) as f_in:
            with open(output_file, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

        logger.info(f"Reference genome saved to {output_file}")
        return str(output_file)

    except requests.exceptions.RequestException as e:
        logger.error(f"Error downloading reference genome: {e}")
        raise


def handle_gene_download_nucleotide(
    download_info: Tuple[str, str, int, str, Path], fasta_dir: Path
) -> Optional[str]:
    """
    Helper function to download a gene sequence.

    Args:
        download_info: Tuple containing (symbol, seq_id, index, disease_term, output_dir)
        fasta_dir: Directory where FASTA files are stored

    Returns:
        Path to downloaded file or None if download failed
    """
    symbol, seq_id, i, disease_term, output_dir = download_info

    try:
        # Generate a descriptive filename
        filename = f"{symbol}_transcript{i+1}_{seq_id}.fasta"
        file_path = os.path.join(output_dir, filename)

        # Skip if already downloaded
        if os.path.exists(file_path):
            logger.info(f"File already exists: {file_path}")
            return str(file_path)

        logger.info(f"Downloading sequence for nucleotide ID: {seq_id}")

        # Download sequence
        seq_handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")

        # Handle content as both string or bytes
        content = seq_handle.read()
        if isinstance(content, bytes):
            content = content.decode("utf-8")

        with open(file_path, "w", encoding="utf-8") as out_f:
            out_f.write(content)

        seq_handle.close()

        logger.info(f"Downloaded {seq_id} for gene {symbol} to {file_path}")
        log_download("Disease", f"{disease_term}:{symbol}", str(file_path), fasta_dir)

        # Be nice to NCBI by adding a delay
        time.sleep(2)
        return str(file_path)

    except IOError as e:
        logger.error(f"File IO error downloading {seq_id}: {str(e)}")
        return None
    except Exception as e:
        logger.error(f"Failed to download {seq_id}: {str(e)}")
        return None


def get_gene_symbol(gene_id: str) -> str:
    """
    Get the gene symbol for a given gene ID.

    Args:
        gene_id: NCBI Gene ID

    Returns:
        Gene symbol or a default string with the gene ID
    """
    try:
        gene_handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
        gene_record = Entrez.read(gene_handle)
        gene_handle.close()

        gene_symbol = "unknown"
        if gene_record and "Entrezgene" in gene_record[0]:
            gene_info = gene_record[0]["Entrezgene"]
            gene_symbol = gene_info.get("Gene-ref", {}).get("Gene-ref_locus", f"gene_{gene_id}")

        return gene_symbol
    except Exception as e:
        logger.error(f"Error getting gene symbol for {gene_id}: {e}")
        return f"gene_{gene_id}"


def search_nucleotide_by_gene_id(gene_id: str, max_results: int = 3) -> List[str]:
    """
    Search for nucleotide sequences associated with a gene ID.

    Args:
        gene_id: NCBI Gene ID
        max_results: Maximum number of results to return

    Returns:
        List of nucleotide IDs
    """
    search_term = (
        f"{gene_id}[Gene ID] AND refseq[Filter] AND " f"mRNA[Filter] AND Homo sapiens[Organism]"
    )

    try:
        direct_handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=max_results)
        direct_results = Entrez.read(direct_handle)
        direct_handle.close()

        return direct_results.get("IdList", [])
    except Exception as e:
        logger.error(f"Error searching nucleotide database: {e}")
        return []


def search_nucleotide_by_symbol(gene_symbol: str, max_results: int = 3) -> List[str]:
    """
    Search for nucleotide sequences associated with a gene symbol.

    Args:
        gene_symbol: Gene symbol
        max_results: Maximum number of results to return

    Returns:
        List of nucleotide IDs
    """
    symbol_search = (
        f"{gene_symbol}[Gene Name] AND refseq[Filter] AND "
        f"mRNA[Filter] AND Homo sapiens[Organism]"
    )

    try:
        symbol_handle = Entrez.esearch(db="nucleotide", term=symbol_search, retmax=max_results)
        symbol_results = Entrez.read(symbol_handle)
        symbol_handle.close()

        return symbol_results.get("IdList", [])
    except Exception as e:
        logger.error(f"Error searching nucleotide database by symbol: {e}")
        return []


def download_disease_related_genes_clinvar(
    disease_term: str,
    max_results: int = 1,
    output_dir: Optional[Path] = None,
    samples_dir: Path = PROGRAM_STORAGE_DIR_SHARED_DATA_FASTA_SAMPLES,
    fasta_dir: Path = PROGRAM_STORAGE_DIR_SHARED_DATA_FASTA,
) -> List[str]:
    """
    Search for and download genes associated with a specific disease.

    Args:
        disease_term: Disease name or search term (e.g., "breast cancer")
        max_results: Maximum number of gene sequences to download
        output_dir: Directory to save files (defaults to samples dir)
        samples_dir: Base directory for sample files
        fasta_dir: Directory where FASTA files are stored

    Returns:
        List of paths to downloaded files
    """
    print(f"Downloading disease-related genes for {disease_term}...")
    if output_dir is None:
        output_dir = os.path.join(samples_dir, "diseases", disease_term.replace(" ", "_").lower())
        os.makedirs(output_dir, exist_ok=True)

    downloaded_files = []

    try:
        logger.info(f"Searching for genes related to {disease_term}...")

        # First search for gene IDs related to the disease
        search_handle = Entrez.esearch(
            db="gene",
            term=f"((((\"{disease_term}\"[Disease\/Phenotype]) AND Pathogenic) OR likely_pathogenic) OR risk_factor) AND \"homo sapiens\"[Organism]",
            retmax=max_results,
        )
        search_results = Entrez.read(search_handle)
        search_handle.close()

        gene_ids = search_results["IdList"]

        if not gene_ids:
            logger.warning(f"No genes found for term: {disease_term}")
            return downloaded_files

        logger.info(f"Found {len(gene_ids)} genes related to {disease_term}")

        # For each gene, get the nucleotide sequences
        for gene_id in gene_ids:
            try:
                logger.info(f"Processing gene ID: {gene_id}")

                # Direct search for nucleotide sequences by gene ID
                nucleotide_ids = search_nucleotide_by_gene_id(gene_id)

                gene_symbol = "unknown"
                if not nucleotide_ids:
                    logger.info(f"No direct nucleotide sequences found for gene ID: {gene_id}")

                    # Try to get gene symbol
                    gene_symbol = get_gene_symbol(gene_id)
                    logger.info(f"Found gene symbol: {gene_symbol}")

                    # Try searching by gene symbol
                    nucleotide_ids = search_nucleotide_by_symbol(gene_symbol)
                    logger.info(f"Found {len(nucleotide_ids)} nucleotide sequences by gene symbol")
                else:
                    logger.info(f"Found {len(nucleotide_ids)} nucleotide sequences by gene ID")
                    gene_symbol = get_gene_symbol(gene_id)

                # If we still don"t have any nucleotide IDs, skip this gene
                if not nucleotide_ids:
                    logger.warning(f"No nucleotide sequences found for gene ID: {gene_id}")
                    continue

                # Download each nucleotide sequence
                for i, nucleotide_id in enumerate(nucleotide_ids):
                    download_info = (gene_symbol, nucleotide_id, i, disease_term, output_dir)
                    file_path = handle_gene_download_nucleotide(download_info, fasta_dir)
                    if file_path:
                        downloaded_files.append(file_path)

            except Exception as e:
                logger.error(f"Error processing gene {gene_id}: {str(e)}")
                continue

        return downloaded_files

    except Exception as e:
        logger.error(f"Error searching for disease genes: {str(e)}")
        raise
    
if __name__ == "__main__":
    # Example usage
    disease_term = "breast cancer"
    max_results = 100
    downloaded_files = download_disease_related_genes_clinvar(disease_term, max_results=max_results)
    print(f"Downloaded files: {downloaded_files}")