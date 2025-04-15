import json
import os
import re
import sys

sys.path.append(
    os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
)
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
from services.utils.hgvs_parser import parse_hgvs, format_to_hgvs, is_valid_hgvs


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

    if not re.match(
        r"^(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$", str(chromosome), re.IGNORECASE
    ):
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
            logger.error(
                f"MyVariant.info batch response was not a list: {type(results)}"
            )
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
            logger.error(
                f"Response Text: {response.text[:500]}..."
            )  # Log beginning of text
        except NameError:
            logger.error("Response object not available for logging text.")
        return None


def extract_clinvar_diseases(api_result: List[Dict]) -> List[Dict[str, Any]]:
    """Extracts disease information from ClinVar data in MyVariant.info response.
    
    Args:
        api_result (List[Dict]): The response from MyVariant.info API
        
    Returns:
        List[Dict[str, Any]]: List of dictionaries containing disease information
    """
    results = []

    if not api_result:
        logger.warning("No API result to process")
        return results

    track_api_queries(api_result)

    for entry_batch in api_result:
        not_found = entry_batch.get("notfound", False)
        if not_found:
            continue

        # Get variant information from HGVS
        hgvs_id = entry_batch.get("_id", "")
        parsed_hgvs = parse_hgvs(hgvs_id)
        
        chromosome = None
        position = None
        reference = None
        alternate = None
        
        if parsed_hgvs:
            chromosome = parsed_hgvs.get('chromosome')
            position = parsed_hgvs.get('position')
            reference = parsed_hgvs.get('reference')
            alternate = parsed_hgvs.get('alternate')
        else:
            # Legacy fallback method if HGVS parsing fails
            chromosome = hgvs_id.split(":g.")[0].replace("chr", "") if ":g." in hgvs_id else None
            
            if ">" in hgvs_id:  # SNP/MNP format
                position_part = hgvs_id.split(":g.")[1]
                position_match = re.search(r'(\d+)', position_part)
                position = int(position_match.group(1)) if position_match else None
                
                alleles_match = re.search(r'(\d+)([ACGT]+)>([ACGT]+)', position_part)
                if alleles_match:
                    reference = alleles_match.group(2)
                    alternate = alleles_match.group(3)
        
        json_data = entry_batch.get("clinvar", {})
        if not json_data:
            continue
            
        # Try to fill in missing data from clinvar fields if we couldn't get it from HGVS
        if position is None and "hg38" in json_data:
            position = json_data.get("hg38", {}).get("start")
        if reference is None:
            reference = json_data.get("ref", "")
        if alternate is None:
            alternate = json_data.get("alt", "")
                     
        # Process RCV records (Reference ClinVar records)
        if "rcv" in json_data:
            rcv_data = json_data["rcv"]
            # Handle both list and single object formats for rcv
            if isinstance(rcv_data, list):
                rcv_list = rcv_data
            else:
                # If it's a single record, wrap it in a list for unified processing
                rcv_list = [rcv_data]
                
            for rcv in rcv_list:
                try:
                    clinical_significance = rcv.get("clinical_significance", "")
                    if not validate_significance(clinical_significance):
                        continue
                        
                    # Handle both string and dictionary formats for conditions
                    conditions = rcv.get("conditions", [])
                    if isinstance(conditions, dict):
                        disease_names = [conditions.get("name", "")]
                        synonyms = conditions.get("synonyms", [])
                    elif isinstance(conditions, list):
                        disease_names = [cond.get("name", "") if isinstance(cond, dict) else cond for cond in conditions]
                        synonyms = []
                        for cond in conditions:
                            if isinstance(cond, dict) and "synonyms" in cond:
                                synonyms.extend(cond.get("synonyms", []))
                    else:
                        disease_names = [str(conditions)]
                        synonyms = []
                        
                    # Process each disease name
                    for disease_name in disease_names:
                        if not validate_disease_name(disease_name):
                            continue
                            
                        disease_info = {
                            "clinical_significance": clinical_significance,
                            "disease_name": disease_name,
                            "synonyms": synonyms,
                            "chromosome": chromosome,
                            "position": position,
                            "reference": reference,
                            "alternate": alternate,
                            "hgvs_id": hgvs_id
                        }
                        results.append(disease_info)
                except Exception as e:
                    logger.error(f"Error processing RCV record: {e}")
                    continue
        
        # Also check for "clinvar.disease_names" which might be present in some responses
        elif "disease_names" in json_data:
            for disease_name in json_data["disease_names"]:
                if validate_disease_name(disease_name):
                    # Try to get clinical significance from other fields if available
                    clinical_significance = json_data.get("clinical_significance", "Pathogenic")
                    
                    disease_info = {
                        "clinical_significance": clinical_significance,
                        "disease_name": disease_name,
                        "synonyms": [],
                        "chromosome": chromosome,
                        "position": position,
                        "reference": reference,
                        "alternate": alternate,
                        "hgvs_id": hgvs_id
                    }
                    results.append(disease_info)

    logger.info(f"Extracted {len(results)} disease associations")
    return results


def shrink_and_count_duplicates(df: pd.DataFrame) -> pd.DataFrame:
    """Shrinks the DataFrame to only unique rows and counts duplicates.

    Args:
        df (pd.DataFrame): Input DataFrame.

    Returns:
        pd.DataFrame: DataFrame with unique rows and a count of duplicates.
    """
    df = df.groupby(list(df.columns)).size().reset_index(name="count")
    df = df.drop_duplicates()  # Use drop_duplicates() instead of remove_duplicates()
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


def load_variants_from_json(json_path: str) -> pd.DataFrame:
    """Loads variants from a JSON file.

    Args:
        json_path (str): Path to the JSON file.

    Returns:
        pd.DataFrame: DataFrame with variant data.
    """

    important_subset = ["chromosome", "position", "reference_allele", "query_allele"]
    
    logger.info(f"Loading variants from {json_path}")

    try:
        data = pd.read_json(json_path)

        data = data[important_subset]
        data = data.drop_duplicates(subset=important_subset)

        data["chromosome"] = data["chromosome"].astype(str)
        data["position"] = data["position"].fillna(0).astype(int)

        # Rename columns to match expected format
        if "query_allele" in data.columns:
            data = data.rename(columns={"query_allele": "alternate"})

        if "reference_allele" in data.columns:
            data = data.rename(columns={"reference_allele": "reference"})

        required_columns = ["chromosome", "position", "reference", "alternate"]
        if not all(col in data.columns for col in required_columns):
            logger.error(
                f"JSON is missing required columns. Required: {required_columns}"
            )
            return pd.DataFrame()

        logger.info(f"Loaded {len(data)} variants from {json_path}")
        return data[required_columns]
    except Exception as e:
        logger.error(f"Error loading variants from JSON: {e}")
        return pd.DataFrame()


def process_variants(variants_file=str, batch_size=MYVARIANT_BATCH_SIZE):
    """Main function to process variants in batches.

    Args:
        variants_file (str): Path to the variants JSON file
        batch_size (int): Number of variants to process in each batch
        
    Returns:
        List[Dict[str, Any]]: List of disease associations found
    """
    variants = load_variants_from_json(variants_file)

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

            hgvs_id = format_to_hgvs(chromosome, position, reference, alternate)
            valid_hgvs.append(hgvs_id)

        if not valid_hgvs:
            logger.warning(
                f"No valid variants in batch {batch_start//batch_size + 1}, skipping"
            )
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
        "gene_7157_transcript1_2246031143_vs_GRCh38_direct_variations.json",
    )
    
    file_list = os.listdir(PROGRAM_STORAGE_DIR_SHARED_BLAST)
    
    jsons_only = [
        file for file in file_list if file.endswith(".json")
    ]
    
    # jsons_only.append(os.path.join(PROGRAM_STORAGE_DIR_SHARED_BLAST, "unknown_transcript1_1890272221_vs_GRCh38_direct_variations.json"))
    
    if not jsons_only:
        print("No JSON files found in the directory")
        return
    
    for json_file in jsons_only:
        variants_file = os.path.join(PROGRAM_STORAGE_DIR_SHARED_BLAST, json_file)
        print(f"Processing file: {variants_file}")
        things = process_variants(variants_file)
        if things:
            print(f"Found {len(things)} disease associations in {variants_file}")
        else:
            print(f"No disease associations found in {variants_file}")


if __name__ == "__main__":
    main()
