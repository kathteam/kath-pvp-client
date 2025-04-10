import json
import os
import re
import sys
import time
from logging import Logger
from typing import Any, Dict, List, Optional, Union
import pandas as pd
import requests
from Bio import Entrez

from shared.constants import (
    PROGRAM_STORAGE_DIR_SHARED_BLAST,
    PROGRAM_STORAGE_DIR_SHARED_DATA_DISEASES,
)
from utils.logger import get_logger

backend_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../"))
if backend_dir not in sys.path:
    sys.path.append(backend_dir)

os.makedirs(PROGRAM_STORAGE_DIR_SHARED_DATA_DISEASES, exist_ok=True)
os.makedirs(PROGRAM_STORAGE_DIR_SHARED_BLAST, exist_ok=True)

logger: Logger = get_logger(__name__)

MYVARIANT_API_URL = "https://myvariant.info/v1/variant"
MYVARIANT_BATCH_SIZE = 1000
MYVARIANT_REQUEST_TIMEOUT = 120
API_DELAY_SECONDS = 1.0


def track_api_queries(response) -> None:
    """Tracks API queries and their responses."""
    folder = os.path.join(PROGRAM_STORAGE_DIR_SHARED_BLAST, "api_queries")
    os.makedirs(folder, exist_ok=True)
    file_number = len(os.listdir(folder)) + 1
    file_path = os.path.join(folder, f"api_query_{file_number}.json")
    with open(file_path, "w", encoding="utf-8") as f:
        json.dump(response, f, indent=4)


def validate_disease_name(disease_name: str) -> bool:
    """Validates the disease name.

    Args:
        disease_name (str): Disease name.

    Returns:
        bool: True if the disease name is valid, False otherwise.
    """

    invalid_names = [
        "unknown",
        "not specified",
        "unspecified",
        "not provided",
        "not available",
        "not applicable",
        "not determined",
    ]
    if disease_name.lower() in invalid_names:
        return False

    if not isinstance(disease_name, str) or not disease_name.strip():
        return False

    return True


def validate_significance(clinical_significance: str) -> bool:
    """Validates the clinical significance.

    Args:
        clinical_significance (str): Clinical significance.

    Returns:
        bool: True if the clinical significance is valid, False otherwise.
    """
    invalid_names = [
        "unknown",
        "not specified",
        "unspecified",
        "Uncertain significance",
        "association",
    ]
    if clinical_significance.lower() in invalid_names:
        logger.error(f"Invalid clinical significance: {clinical_significance}")
        return False
    return True


def validate_variant(
    chromosome: Union[int, str], position: int, reference: str, alternate: str
) -> bool:
    """Checks if the variant data components are plausible."""

    if not re.match(r"^(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$", str(chromosome), re.IGNORECASE):
        logger.debug(f"Invalid chromosome format: {chromosome}")
        return False
    # Allow nucleotides and potentially '-' for indel representation in REF/ALT
    if not re.match(r"^[ACGTN-]+$", reference, re.IGNORECASE) or not reference:
        logger.debug(f"Invalid reference allele format: {reference}")
        return False
    if not re.match(r"^[ACGTN-]+$", alternate, re.IGNORECASE) or not alternate:
        logger.debug(f"Invalid alternate allele format: {alternate}")
        return False
    if not isinstance(position, int) or position <= 0:
        logger.debug(f"Invalid position: {position}")
        return False
    return True


def query_into_myvariant(hgvs_list: List[str]) -> Optional[Dict[str, Any]]:
    """Queries MyVariant.info in batches using POST requests."""
    if not hgvs_list:
        return []

    # Add more fields if needed (e.g., dbsnp)
    url = f"{MYVARIANT_API_URL}?assembly=hg38&fields=clinvar"
    data_payload = {"ids": hgvs_list}
    headers = {"Content-Type": "application/json"}

    logger.debug(f"POSTing batch of {len(hgvs_list)} variants to MyVariant.info")

    try:
        response = requests.post(
            url, headers=headers, json=data_payload, timeout=MYVARIANT_REQUEST_TIMEOUT
        )
        response.raise_for_status()
        results = response.json()
        if not isinstance(results, list):
            logger.error(f"MyVariant.info batch response was not a list: {type(results)}")
            return None
        logger.debug(f"Received {len(results)} results from batch.")
        return results
    except requests.exceptions.HTTPError as e:
        if e.response.status_code != 404:
            logger.error(f"HTTP error during MyVariant.info batch query: {e}")
        else:
            logger.warning(
                "MyVariant.info batch query returned 404 "
                + "(likely all variants not found or issue with request format)."
            )
        return None
    except requests.exceptions.RequestException as e:
        logger.error(f"Network error during MyVariant.info batch query: {e}")
        return None
    except json.JSONDecodeError:
        logger.error("Error decoding JSON response from MyVariant.info batch query.")
        try:
            logger.error(f"Response Text: {response.text[:500]}...")  # Log beginning of text
        except NameError:
            logger.error("Response object not available for logging text.")
        return None


def format_variant_hgvs(
    chromosome: Union[int, str], position: int, reference: str, alternate: str
) -> str:
    """Formats variant components into an HGVS string for MyVariant.info."""
    chrom_str = str(chromosome).upper().replace("CHR", "")
    if reference == "-":
        return f"chr{chrom_str}:g.{position}_{position+1}ins{alternate}"
    if alternate == "-":
        if len(reference) == 1:
            return f"chr{chrom_str}:g.{position}del"
        return f"chr{chrom_str}:g.{position}_{position + len(reference) - 1}del"
    # SNP / MNP
    return f"chr{chrom_str}:g.{position}{reference}>{alternate}"


def extract_clinvar_diseases(api_result: Dict) -> List[Dict[str, Any]]:

    results = []

    track_api_queries(api_result)

    for entry_batch in api_result:

        not_found = entry_batch.get("notfound", False)
        if not_found:
            continue

        json_data = entry_batch.get("clinvar", {})
        if not json_data:
            logger.warning("No ClinVar data found")
            return []

        chromosome = json_data.get("chrom", "")
        position = json_data.get("hg38", "").get("start", "")
        reference = json_data.get("reference", "")
        alternate = json_data.get("alt", "")

        for entry in json_data["rcv"]:
            try:
                clinical_significance = entry.get("clinical_significance", "")
                if not validate_significance(clinical_significance):
                    continue

                conditions = entry.get("conditions", [])
                if not conditions:
                    continue

                disease_name = conditions.get("name", "")
                # synonims = conditions.get("synonyms", [])
                if not validate_disease_name(disease_name):
                    continue

                disease_info = {
                    "clinical_significance": clinical_significance,
                    "disease_name": disease_name,
                    # "synonyms": synonims,
                    "chromosome": chromosome,
                    "position": position,
                    "reference": reference,
                    "alternate": alternate,
                }

                results.append(disease_info)

            except KeyError as e:
                logger.error(f"Key error in ClinVar data: {e}")
                continue

        print(results)
        return results


def shrink_and_count_duplicates(df: pd.DataFrame) -> pd.DataFrame:
    """Shrinks the DataFrame to only unique rows and counts duplicates.

    Args:
        df (pd.DataFrame): Input DataFrame.

    Returns:
        pd.DataFrame: DataFrame with unique rows and a count of duplicates.
    """
    df = df.groupby(list(df.columns)).size().reset_index(name="count")
    df = df.remove_duplicates()
    df["count"] = df["count"].astype(int)
    return df


def save_to_csv(data: List[Dict[str, Any]], output_path: str):
    """Saves results to a CSV file.

    Args:
        data (List[Dict[str, Any]]): List of dictionaries with data. JSON.
        output_path (str): Output file path.
    """
    try:
        df = pd.DataFrame(data)
        df.to_csv(output_path, index=False)
        logger.info(f"Saved results to {output_path}")
    except Exception as e:
        logger.error(f"Error saving results to CSV: {e}")


def load_variants_from_csv(csv_path: str) -> pd.DataFrame:
    """Loads variants from a CSV file.

    Args:
        csv_path (str): Path to the CSV file.

    Returns:
        pd.DataFrame: DataFrame with variant data.
    """
    try:
        data = pd.read_csv(csv_path)
        data = data.dropna(subset=["chromosome", "position", "reference", "query"])
        data = data.drop_duplicates(subset=["chromosome", "position", "reference", "query"])

        data["chromosome"] = data["chromosome"].astype(int)
        data["position"] = data["position"].astype(int)

        # Rename query to alternate for consistency
        if "query" in data.columns and "alternate" not in data.columns:
            data = data.rename(columns={"query": "alternate"})

        required_columns = ["chromosome", "position", "reference", "alternate"]
        if not all(col in data.columns for col in required_columns):
            logger.error(f"CSV is missing required columns. Required: {required_columns}")
            return pd.DataFrame()

        logger.info(f"Loaded {len(data)} variants from {csv_path}")
        return data[required_columns]
    except Exception as e:
        logger.error(f"Error loading variants from CSV: {e}")
        return pd.DataFrame()


def process_variants(variants_file=str, batch_size=MYVARIANT_BATCH_SIZE):
    """Main function to process variants in batches.

    Args:
        batch_size (int): Number of variants to process in each batch
    """
    variants = load_variants_from_csv(variants_file)

    if variants.empty:
        logger.error("No valid variants to process")
        return []

    all_disease_data = []
    total_variants = len(variants)
    logger.info(f"Processing {total_variants} variants in batches of {batch_size}")

    # Process variants in batches
    for batch_start in range(0, total_variants, batch_size):
        batch_end = min(batch_start + batch_size, total_variants)
        batch_variants = variants.iloc[batch_start:batch_end]

        logger.info(
            f"Processing batch {batch_start//batch_size + 1}: variants "
            + f"{batch_start+1}-{batch_end} of {total_variants}"
        )

        valid_hgvs = []

        for index, variant in batch_variants.iterrows():
            chromosome = variant["chromosome"]
            position = variant["position"]
            reference = variant["reference"]
            alternate = variant["alternate"]

            if not validate_variant(chromosome, position, reference, alternate):
                logger.warning(f"Invalid variant data at index {index}, skipping")
                continue

            hgvs_id = format_variant_hgvs(chromosome, position, reference, alternate)
            valid_hgvs.append(hgvs_id)

        if not valid_hgvs:
            logger.warning(f"No valid variants in batch {batch_start//batch_size + 1}, skipping")
            continue

        logger.info(f"Querying MyVariant.info with {len(valid_hgvs)} valid variants")
        data = query_into_myvariant(valid_hgvs)

        if data:
            diseases = extract_clinvar_diseases(data)
            if diseases:
                all_disease_data.extend(diseases)
                logger.info(f"Found {len(diseases)} disease associations in batch")
            else:
                logger.info("No disease associations found in batch")

        if batch_end < total_variants:
            logger.debug(f"Sleeping for {API_DELAY_SECONDS} seconds before next batch")
            time.sleep(API_DELAY_SECONDS)

    if all_disease_data:
        timestamp = time.strftime("%Y%m%d-%H%M%S")
        output_file = os.path.join(
            PROGRAM_STORAGE_DIR_SHARED_DATA_DISEASES,
            f"disease_associations_{timestamp}.csv",
        )
        save_to_csv(all_disease_data, output_file)
        logger.info(f"Found total of {len(all_disease_data)} disease associations")
    else:
        logger.warning("No disease associations found in any batch")

    return all_disease_data


def main():
    # Configure NCBI Entrez
    Entrez.email = "kajeliukasc@gmail.com"
    Entrez.api_key = "7ca5ef526507701d64f16a090124cbc4aa08"

    variants_file = os.path.join(
        PROGRAM_STORAGE_DIR_SHARED_BLAST,
        "consolidated_results_unknown_transcript2_2246031136_vs_GRCh38_direct_blast_dir",
        "all_mutations_20250405-185834.csv",
    )
    process_variants(variants_file)


if __name__ == "__main__":
    main()
