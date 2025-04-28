import json
import os
from pathlib import Path
import re
import sys
from xml.etree import ElementTree
import concurrent.futures

sys.path.append(
    os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
)
import time
from typing import Any, Dict, List, Optional, Union
import pandas as pd
import requests
from Bio import Entrez

from shared.constants import (
    PROGRAM_STORAGE_DIR_SHARED_BLAST,
    PROGRAM_STORAGE_DIR_SHARED_DATA_DISEASES,
)
from services.utils.script_setup import EnvSetup, FolderSetup
from services.utils.hgvs_parser import (
    parse_hgvs,
    format_to_hgvs,
    validate_hgvs,
    vcf_to_c_hgvs,
)

EnvSetup()  # Sets up the environment variables ENTREZ included
FolderSetup()

from utils.logger import get_logger

logger = get_logger(__name__)


def track_api_queries(response_xml: str) -> None:
    """Tracks API queries and their responses, saving them as XML."""
    print(response_xml)
    folder = os.path.join(PROGRAM_STORAGE_DIR_SHARED_BLAST, "api_queries")
    os.makedirs(folder, exist_ok=True)
    file_number = len(os.listdir(folder)) + 1
    file_path = os.path.join(folder, f"api_query_{file_number}.xml")
    try:
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(response_xml)
        logger.debug(f"Saved API query response to {file_path}")
    except Exception as e:
        logger.error(f"Error saving API query response to XML: {e}")


def batch_fetch_clinvar_records(
    clinvar_ids: List[str],
) -> Dict[str, List[Dict[str, str]]]:
    """Fetch multiple ClinVar records by ID and extract clinical significance."""
    if not clinvar_ids:
        return {}

    results_by_id = {}
    batch_size = 200
    unique_ids = list(set(clinvar_ids))

    for i in range(0, len(unique_ids), batch_size):
        id_batch = unique_ids[i : i + batch_size]
        ids_str = ",".join(id_batch)
        try:
            logger.debug(f"Fetching summaries for {len(id_batch)} ClinVar IDs.")
            fetch_handle = Entrez.esummary(db="clinvar", id=ids_str, rettype="xml")
            xml_data = fetch_handle.read()
            fetch_handle.close()
            time.sleep(0.4)

            root = ElementTree.fromstring(xml_data)
            for docsum in root.findall(".//DocumentSummary"):
                uid = docsum.get("uid")
                if not uid:
                    continue

                record_data = []
                for gremline_class in docsum.iterfind("germline_classification"):
                    significance = gremline_class.findtext("description")
                    trait_name_text = "Unknown"
                    trait_set = gremline_class.find("trait_set")
                    if trait_set is not None:
                        for trait in trait_set.findall("trait"):
                            trait_name = trait.find("trait_name")
                            if trait_name is not None and trait_name.text:
                                trait_name_text = trait_name.text
                                break

                    record_data.append(
                        {
                            "significance": significance or "Not provided",
                            "trait": trait_name_text,
                        }
                    )
                if record_data:
                    results_by_id[uid] = record_data

        except Exception as e:
            logger.error(f"Error fetching/parsing ClinVar summaries for batch: {e}")
        finally:
            if "fetch_handle" in locals() and fetch_handle:
                fetch_handle.close()

    return results_by_id


def batch_search_ncbi_clinvar(hgvs_list: List[str]) -> Dict[str, List[Dict[str, str]]]:
    """Searches ClinVar for a list of HGVS strings and returns mapped results using parallel requests."""
    if not hgvs_list:
        return {}

    hgvs_to_ids = {}
    all_unique_ids = set()
    valid_hgvs_list = [h for h in hgvs_list if validate_hgvs(h)]

    logger.info(f"Searching ClinVar IDs for {len(valid_hgvs_list)} HGVS terms using parallel requests...")

    # Helper function to search for a single HGVS term
    def search_single_hgvs(hgvs: str) -> tuple[str, list[str]]:
        try:
            search_handle = Entrez.esearch(db="clinvar", term=hgvs, retmax=100)
            search_record = Entrez.read(search_handle)
            search_handle.close()
            time.sleep(0.4)  # Respect NCBI rate limits per request

            id_list = search_record.get("IdList", [])
            if id_list:
                logger.debug(f"Found {len(id_list)} IDs for {hgvs}")
                return hgvs, id_list
            else:
                logger.debug(f"No IDs found for {hgvs}")
                return hgvs, []
        except Exception as e:
            logger.error(f"Error searching ClinVar for {hgvs}: {e}")
            return hgvs, []
        finally:
            if 'search_handle' in locals() and search_handle:
                try:
                    search_handle.close()
                except Exception:
                    pass

    max_workers = 5
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_hgvs = {executor.submit(search_single_hgvs, hgvs): hgvs for hgvs in valid_hgvs_list}

        processed_count = 0
        total_tasks = len(future_to_hgvs)
        for future in concurrent.futures.as_completed(future_to_hgvs):
            hgvs, id_list = future.result()
            hgvs_to_ids[hgvs] = id_list
            all_unique_ids.update(id_list)
            processed_count += 1
            if processed_count % 100 == 0:
                logger.info(f"Completed ClinVar search for {processed_count}/{total_tasks} HGVS terms.")

    for hgvs in hgvs_list:
        if hgvs not in hgvs_to_ids:
            hgvs_to_ids[hgvs] = []

    logger.info(
        f"Found {len(all_unique_ids)} unique ClinVar IDs across {len(hgvs_to_ids)} HGVS terms."
    )

    if not all_unique_ids:
        return {hgvs: [] for hgvs in hgvs_list}

    id_to_data = batch_fetch_clinvar_records(list(all_unique_ids))

    results = {hgvs: [] for hgvs in hgvs_list}
    for hgvs, id_list in hgvs_to_ids.items():
        hgvs_results = []
        for clinvar_id in id_list:
            if clinvar_id in id_to_data:
                hgvs_results.extend(id_to_data[clinvar_id])
        unique_hgvs_results = [
            dict(t) for t in {tuple(d.items()) for d in hgvs_results}
        ]
        results[hgvs] = unique_hgvs_results

    return results


def save_to_csv(data: List[Dict[str, Any]], output_path: str):
    """Saves results to a CSV file.

    Args:
        data (List[Dict[str, Any]]): List of dictionaries with data. JSON.
        output_path (str): Output file path.
    """

    print(data)

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

    important_subset = [
        "transcript_id",
        "chromosome",
        "position",
        "reference_allele",
        "query_allele",
    ]

    logger.info(f"Loading variants from {json_path}")

    try:
        data = pd.read_json(json_path)

        regex_pattern = r"(?P<accession>[A-Z]{2}_\d+\.\d+)"
        data["transcript_id"] = data["query_id"].str.extract(regex_pattern)

        data = data[important_subset]
        data = data.drop_duplicates(subset=important_subset)

        data["chromosome"] = data["chromosome"].astype(str)
        data["position"] = data["position"].fillna(0).astype(int)

        if "query_allele" in data.columns:
            data = data.rename(columns={"query_allele": "alternate"})

        if "reference_allele" in data.columns:
            data = data.rename(columns={"reference_allele": "reference"})

        required_columns = [
            "transcript_id",
            "chromosome",
            "position",
            "reference",
            "alternate",
        ]
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


def process_variant(variants_file: str):
    """Processes variants from a single file by generating HGVS notations and querying ClinVar."""
    variants = load_variants_from_json(variants_file)
    if variants.empty:
        logger.warning(f"No variants loaded from {variants_file}. Skipping.")
        return []

    hgvs_list = []
    processed_count = 0
    total_variants = len(variants)

    logger.info(f"Generating HGVS notations for {total_variants} variants...")
    for index, variant in variants.iterrows():
        transcript_id = variant.get("transcript_id", None)
        chromosome = variant.get("chromosome", None)
        position = variant.get("position", None)
        reference = variant.get("reference", None)
        alternate = variant.get("alternate", None)

        if not all([transcript_id, position is not None, reference, alternate]):
            logger.warning(
                f"Skipping variant at index {index} due to missing data: {variant.to_dict()}"
            )
            continue

        try:
            hgvs = vcf_to_c_hgvs(
                transcript_id=transcript_id,
                c_pos=position,
                ref=reference,
                alt=alternate,
            )
            if validate_hgvs(hgvs):
                hgvs_list.append(hgvs)
            else:
                logger.warning(
                    f"Invalid HGVS generated: {hgvs} from variant {index}. Skipping."
                )

        except Exception as e:
            logger.error(f"Error converting variant {index} to HGVS: {e}. Skipping.")

        processed_count += 1
        if processed_count % 1000 == 0:
            logger.info(
                f"Processed {processed_count}/{total_variants} variants for HGVS generation."
            )

    unique_hgvs_list = list(set(hgvs_list))
    logger.info(f"Generated {len(unique_hgvs_list)} unique valid HGVS strings.")

    if not unique_hgvs_list:
        logger.warning("No valid HGVS strings generated. No ClinVar search performed.")
        return []

    hgvs_disease_map = batch_search_ncbi_clinvar(unique_hgvs_list)

    all_results_for_file = []
    logger.info("Saving results...")
    for hgvs, diseases in hgvs_disease_map.items():
        if diseases:
            output_path = os.path.join(
                PROGRAM_STORAGE_DIR_SHARED_DATA_DISEASES, f"{hgvs}.csv"
            )
            save_to_csv(diseases, output_path)
            for disease_info in diseases:
                all_results_for_file.append({"hgvs": hgvs, **disease_info})

    logger.info(f"Finished processing {variants_file}.")
    return all_results_for_file


def proccess_multiple_variants(folder_name: str) -> None:
    """Process multiple variants files in a folder."""
    all_disease_data_aggregated = []

    json_files = [f for f in os.listdir(folder_name) if f.endswith(".json")]
    logger.info(f"Found {len(json_files)} JSON files in {folder_name}.")

    for file in json_files:
        file_path = os.path.join(folder_name, file)
        logger.info(f"--- Processing file: {file} ---")
        disease_data_for_file = process_variant(file_path)
        all_disease_data_aggregated.extend(disease_data_for_file)
        logger.info(f"--- Finished processing file: {file} ---")

    if all_disease_data_aggregated:
        aggregated_output_path = os.path.join(
            PROGRAM_STORAGE_DIR_SHARED_DATA_DISEASES, "aggregated_results.csv"
        )
        logger.info(
            f"Saving aggregated results ({len(all_disease_data_aggregated)} records) to {aggregated_output_path}"
        )
        save_to_csv(all_disease_data_aggregated, aggregated_output_path)
    else:
        logger.info("No disease data found across all processed files.")


if __name__ == "__main__":
    variants_folder = PROGRAM_STORAGE_DIR_SHARED_BLAST
    logger.info(f"Starting variant processing in folder: {variants_folder}")
    proccess_multiple_variants(variants_folder)
    logger.info("All processing finished.")
