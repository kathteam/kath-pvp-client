from Bio import Entrez
import json
from typing import Dict, Any
import os
import pandas
from dotenv import load_dotenv

# Initialize logging
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger("analysis-service")


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


dirs = initialize_directories()

import pickle
from pathlib import Path


# Add a cache to avoid repeated conversions and API calls
def get_cache_dir():
    """Get the cache directory"""
    cache_dir = Path(dirs["service_dir"]) / "cache"
    cache_dir.mkdir(exist_ok=True)
    return cache_dir


def get_cache_file(cache_type):
    """Get a cache file path"""
    return get_cache_dir() / f"{cache_type}_cache.pkl"


def load_cache(cache_type):
    """Load a cache"""
    cache_file = get_cache_file(cache_type)
    if cache_file.exists():
        try:
            with open(cache_file, "rb") as f:
                return pickle.load(f)
        except Exception as e:
            logger.error(f"Error loading {cache_type} cache: {e}")
            return {}
    return {}


def save_cache(cache_type, cache_data):
    """Save a cache"""
    cache_file = get_cache_file(cache_type)
    with open(cache_file, "wb") as f:
        pickle.dump(cache_data, f)


def convert_coordinates(chromosome, position, from_build="GRCh38", to_build="GRCh37"):
    """
    Convert genomic coordinates between genome builds using UCSC"s liftOver API via PyLiftover.
    If PyLiftover is not installed, falls back to a manual conversion table for common variants.

    Args:
        chromosome: Chromosome (e.g., "1", "X")
        position: Position in the source build
        from_build: Source genome build (default: "GRCh38")
        to_build: Target genome build (default: "GRCh37")

    Returns:
        Tuple of (converted_chromosome, converted_position) or (None, None) if conversion fails
    """
    # Check cache first
    cache_key = f"{chromosome}:{position}:{from_build}:{to_build}"
    conversion_cache = load_cache("liftover")

    if cache_key in conversion_cache:
        return conversion_cache[cache_key]

    try:
        # Try to use PyLiftover if available
        from pyliftover import LiftOver

        # Map build names to PyLiftover format
        build_map = {"GRCh38": "hg38", "GRCh37": "hg19", "hg38": "hg38", "hg19": "hg19"}

        from_build_lo = build_map.get(from_build, from_build)
        to_build_lo = build_map.get(to_build, to_build)

        lo = LiftOver(from_build_lo, to_build_lo)

        # Ensure chromosome format (with or without "chr" prefix)
        chrom = chromosome.replace("chr", "")

        # Convert position to int
        pos = int(position)

        # Perform liftover
        converted = lo.convert_coordinate(f"chr{chrom}", pos)

        if converted and len(converted) > 0:
            # Get the first (most likely) conversion
            result = (converted[0][0].replace("chr", ""), str(converted[0][1]))

            # Save to cache
            conversion_cache[cache_key] = result
            save_cache("liftover", conversion_cache)

            return result
        else:
            # No conversion found
            conversion_cache[cache_key] = (None, None)
            save_cache("liftover", conversion_cache)
            return (None, None)

    except ImportError:
        # PyLiftover not installed, use fallback method
        logger.warning("PyLiftover not installed. Using fallback manual conversion.")

        # Known conversions for test variants
        manual_conversions = {
            # Example: GRCh38 to GRCh37 for BRCA1 variant
            "17:43045629": "17:43057051",  # BRCA1 (hg38) -> (hg19)
            # Add more known conversions here
        }

        key = f"{chromosome}:{position}"
        if key in manual_conversions:
            # Split the stored value into chromosome and position
            conv = manual_conversions[key].split(":")
            result = (conv[0], conv[1])

            # Save to cache
            conversion_cache[cache_key] = result
            save_cache("liftover", conversion_cache)

            return result

        # No conversion found
        logger.warning(
            f"No manual conversion available for {key}. Install PyLiftover for better results."
        )
        conversion_cache[cache_key] = (None, None)
        save_cache("liftover", conversion_cache)
        return (None, None)


def find_disease_from_variant(
    chromosome,
    position,
    ref_allele,
    alt_allele,
    email="your.email@example.com",
    source_build="GRCh38",
):
    """
    Query ClinVar for disease information, converting coordinates if needed.

    Args:
        chromosome: Chromosome (e.g., "1", "X")
        position: Position in source build
        ref_allele: Reference allele
        alt_allele: Alternate allele
        email: Your email for NCBI
        source_build: Source genome build (default: "GRCh38")

    Returns:
        Dictionary with disease information or None
    """
    # Set NCBI email
    Entrez.email = email

    # Use API key if available
    api_key = os.environ.get("NCBI_API_KEY")
    if api_key:
        Entrez.api_key = api_key

    # First try with original coordinates
    result = query_clinvar(chromosome, position, ref_allele, alt_allele)

    # If no result and using GRCh38, try converting to GRCh37
    if not result and source_build == "GRCh38":
        logger.info(f"No results with GRCh38 coordinates. Trying conversion to GRCh37...")

        # Convert coordinates from GRCh38 to GRCh37
        conv_chrom, conv_pos = convert_coordinates(chromosome, position, "GRCh38", "GRCh37")

        if conv_chrom and conv_pos:
            logger.info(
                f"Converted coordinates: chr{chromosome}:{position} (GRCh38) -> chr{conv_chrom}:{conv_pos} (GRCh37)"
            )
            result = query_clinvar(conv_chrom, conv_pos, ref_allele, alt_allele)
        else:
            logger.warning(f"Could not convert coordinates: chr{chromosome}:{position}")

    return result


def query_clinvar(chromosome, position, ref_allele, alt_allele):
    """
    Query ClinVar API for variant information
    """
    # Check cache first
    cache_key = f"{chromosome}:{position}:{ref_allele}>{alt_allele}"
    clinvar_cache = load_cache("clinvar")

    if cache_key in clinvar_cache:
        logger.info(f"Using cached result for {cache_key}")
        return clinvar_cache[cache_key]

    # Try multiple query formats
    queries = [
        # Format 1: Using genomic notation with "g."
        f"chr{chromosome}:g.{position}{ref_allele}>{alt_allele}",
        # Format 2: Basic coordinates
        f"chr{chromosome}:{position}{ref_allele}>{alt_allele}",
        # Format 3: Using standard search fields
        f"{chromosome}[Chromosome] AND {position}[Base Position]",
        # Format 4: Very specific format
        f"{chromosome}[Chromosome] AND {position}[Base Position] AND {ref_allele}[Reference Allele] AND {alt_allele}[Alternate Allele]",
        # Format 5: Simple text search
        f"chr{chromosome}:{position}",
    ]

    id_list = []
    successful_query = None

    for i, query in enumerate(queries):
        try:
            logger.info(f"Trying query format {i+1}: {query}")
            search_handle = Entrez.esearch(db="clinvar", term=query, retmax=5)
            search_results = Entrez.read(search_handle)
            search_handle.close()

            if search_results.get("Count", "0") != "0":
                id_list = search_results.get("IdList", [])
                if id_list:
                    successful_query = query
                    logger.info(f"Success! Found {len(id_list)} results with query: {query}")
                    break
        except Exception as e:
            logger.error(f"Query format {i+1} failed: {e}")

    if not id_list:
        # No results found
        clinvar_cache[cache_key] = None
        save_cache("clinvar", clinvar_cache)
        return None

    try:
        # Get summary
        summary_handle = Entrez.esummary(db="clinvar", id=id_list[0])
        summary_results = Entrez.read(summary_handle)
        summary_handle.close()

        if not summary_results or not summary_results.get("DocumentSummarySet", {}).get(
            "DocumentSummary", []
        ):
            logger.info("No summary information available")
            clinvar_cache[cache_key] = None
            save_cache("clinvar", clinvar_cache)
            return None

        # Extract information from the summary
        variant_summary = summary_results["DocumentSummarySet"]["DocumentSummary"][0]

        # For more detailed information, get the full record
        fetch_handle = Entrez.efetch(db="clinvar", id=id_list[0], rettype="json")
        fetch_data = json.loads(fetch_handle.read())
        fetch_handle.close()

        # Extract clinical significance and disease names
        clinical_significance = "Unknown"
        disease_names = []

        if fetch_data and "result" in fetch_data:
            variant_data = fetch_data["result"]
            variant_id = list(variant_data.keys())[0]
            variant_info = variant_data[variant_id]

            # Get clinical significance
            if "clinical_significance" in variant_info:
                clinical_significance = variant_info["clinical_significance"].get(
                    "description", "Unknown"
                )

            # Get disease names
            if "trait_set" in variant_info:
                for trait in variant_info["trait_set"]:
                    if "trait_name" in trait:
                        disease_names.append(trait["trait_name"])

        # Construct result dictionary
        result = {
            "clinvar_id": id_list[0],
            "title": variant_summary.get("title", "Unknown"),
            "clinical_significance": clinical_significance,
            "variant_id": variant_summary.get("variation_id", "Unknown"),
            "last_updated": variant_summary.get("update_date", "Unknown"),
            "diseases": (
                disease_names if disease_names else ["No specific disease information available"]
            ),
            "chromosome": chromosome,
            "position": position,
            "ref_allele": ref_allele,
            "alt_allele": alt_allele,
            "query_used": successful_query,
        }

        # Save to cache
        clinvar_cache[cache_key] = result
        save_cache("clinvar", clinvar_cache)

        return result

    except Exception as e:
        logger.error(f"Error processing ClinVar results: {e}")
        clinvar_cache[cache_key] = None
        save_cache("clinvar", clinvar_cache)
        return None


def format_disease_output(result: Dict[str, Any]) -> None:
    """Prints the disease information in a readable format"""
    if not result:
        print("No information available")
        return

    print("\n=== ClinVar Variant Information ===")
    print(
        f"Variant: chr{result["chromosome"]}:{result["position"]} {result["ref_allele"]}>{result["alt_allele"]}"
    )
    print(f"ClinVar ID: {result["clinvar_id"]}")
    print(f"Title: {result["title"]}")
    print(f"Clinical Significance: {result["clinical_significance"]}")
    print("\nAssociated Diseases/Conditions:")
    for disease in result["diseases"]:
        print(f"  - {disease}")
    print(f"\nLast Updated: {result["last_updated"]}")
    print("===================================\n")


def find_mutations_file() -> None:
    """
    Find a specific mutations file in the shared data directory.

    Args:
        filename: Name of the mutations file to find
    """

    must_contain = "all_mutations"

    # Check if the filename contains the required string
    mutation_directory = dirs["blast_dir"]
    consolidated_dir = os.path.join(mutation_directory, "consolidated")

    # for now, its just one
    for filename in os.listdir(consolidated_dir):
        if must_contain in filename:
            data = pandas.read_csv(os.path.join(consolidated_dir, filename))
            for index, row in data.iterrows():
                chromosome = row["chromosome"]
                position = row["position"]
                ref_allele = row["reference"]
                alt_allele = row["query"]

                result = find_disease_from_variant(chromosome, position, ref_allele, alt_allele)
                if result:
                    format_disease_output(result)
                else:
                    print("No information found or an error occurred.")


# At the bottom of the file, modify your main block:

if __name__ == "__main__":
    # Install required packages, if not already installed
    try:
        import pyliftover
    except ImportError:
        print("Installing pyliftover for coordinate conversion...")
        import subprocess

        subprocess.check_call(["pip", "install", "pyliftover"])
        print("pyliftover installed successfully!")

    # Examples in GRCh38 coordinates
    examples = [
        # BRCA1 pathogenic variant (Breast/Ovarian Cancer)
        {
            "chrom": "17",
            "pos": "43045629",
            "ref": "A",
            "alt": "G",
            "disease": "Hereditary breast and ovarian cancer syndrome",
        },
        # CFTR variant (Cystic Fibrosis)
        {"chrom": "7", "pos": "117559593", "ref": "G", "alt": "T", "disease": "Cystic fibrosis"},
        # HBB variant (Sickle Cell)
        {"chrom": "11", "pos": "5225673", "ref": "A", "alt": "T", "disease": "Sickle cell anemia"},
    ]

    # Test with the first example
    example = examples[0]
    chromosome = example["chrom"]
    position = example["pos"]
    ref_allele = example["ref"]
    alt_allele = example["alt"]

    print(
        f"Querying ClinVar for variant chr{chromosome}:{position} {ref_allele}>{alt_allele} (GRCh38)"
    )
    print(f"Expected disease: {example["disease"]}")

    # Replace with your email
    email = os.environ.get("NCBI_API_EMAIL", "your.email@example.com")

    # Query with automatic coordinate conversion
    result = find_disease_from_variant(
        chromosome, position, ref_allele, alt_allele, email, source_build="GRCh38"
    )

    if result:
        format_disease_output(result)
    else:
        print("No information found in ClinVar")
        print("\nThis might be because:")
        print("1. The variant doesn't exist in ClinVar")
        print("2. The coordinate conversion failed")
        print("3. The alleles might be different in the reference genome")
