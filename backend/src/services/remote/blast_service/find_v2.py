import json
import os
from pathlib import Path
import re
import sys
from xml.etree import ElementTree
import concurrent.futures
import threading

sys.path.append(
    os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
)
import time
from typing import Any, Dict, List, Optional, Union, Set, Tuple
import pandas as pd
import requests
from Bio import Entrez
from functools import lru_cache
import hashlib

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

# Create a cache directory
CACHE_DIR = os.path.join(PROGRAM_STORAGE_DIR_SHARED_BLAST, "cache")
os.makedirs(CACHE_DIR, exist_ok=True)

# Constants for optimization
MAX_WORKERS_HGVS = 20  # More workers for CPU-bound HGVS generation
CLINVAR_API_BATCH_SIZE = 200  # Reasonable batch size for efficiency
CLINVAR_CACHE_FILE = os.path.join(CACHE_DIR, "clinvar_cache.json")
HGVS_BATCH_SIZE = 50  # Number of HGVS terms to query in a single request

# Entrez API settings
ENTREZ_EMAIL = os.environ.get("ENTREZ_EMAIL", "your.email@example.com")
ENTREZ_API_KEY = os.environ.get("ENTREZ_API_KEY", None)
ENTREZ_DELAY = 0.34  # Delay between Entrez requests in seconds (0.34 = ~3 requests/second)
if ENTREZ_API_KEY:
    ENTREZ_DELAY = 0.1  # With API key we can do ~10 requests/second

# Configure Entrez
Entrez.email = ENTREZ_EMAIL
if ENTREZ_API_KEY:
    Entrez.api_key = ENTREZ_API_KEY
    logger.info("Using Entrez API key for requests")

# Initialize cache
clinvar_cache = {}
if os.path.exists(CLINVAR_CACHE_FILE):
    try:
        with open(CLINVAR_CACHE_FILE, 'r') as f:
            clinvar_cache = json.load(f)
        logger.info(f"Loaded {len(clinvar_cache)} cached ClinVar entries")
    except Exception as e:
        logger.error(f"Error loading ClinVar cache: {e}")


# Rate limiter for Entrez API
class EntrezRateLimiter:
    def __init__(self, delay=ENTREZ_DELAY):
        self.delay = delay
        self.last_request = 0
        self.lock = threading.Lock()
    
    def wait(self):
        """Wait if needed to respect the rate limit"""
        with self.lock:
            elapsed = time.time() - self.last_request
            if elapsed < self.delay:
                time.sleep(self.delay - elapsed)
            self.last_request = time.time()

entrez_limiter = EntrezRateLimiter()


def save_clinvar_cache():
    """Save the ClinVar cache to disk"""
    try:
        with open(CLINVAR_CACHE_FILE, 'w') as f:
            json.dump(clinvar_cache, f)
        logger.info(f"Saved {len(clinvar_cache)} entries to ClinVar cache")
    except Exception as e:
        logger.error(f"Error saving ClinVar cache: {e}")


def track_api_queries(response_xml: str) -> None:
    """Tracks API queries and their responses, saving them as XML."""
    # Skip printing the whole XML to console
    folder = os.path.join(PROGRAM_STORAGE_DIR_SHARED_BLAST, "api_queries")
    os.makedirs(folder, exist_ok=True)
    
    # Use hash to avoid duplicate files
    hash_value = hashlib.md5(response_xml.encode()).hexdigest()[:10]
    file_path = os.path.join(folder, f"api_query_{hash_value}.xml")
    
    if not os.path.exists(file_path):
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

    # Filter out IDs that are already in cache
    ids_to_fetch = [id for id in clinvar_ids if id not in clinvar_cache]
    logger.info(f"Fetching {len(ids_to_fetch)} new ClinVar IDs (skipping {len(clinvar_ids) - len(ids_to_fetch)} cached)")

    if ids_to_fetch:
        results_by_id = {}
        unique_ids = list(set(ids_to_fetch))

        for i in range(0, len(unique_ids), CLINVAR_API_BATCH_SIZE):
            id_batch = unique_ids[i : i + CLINVAR_API_BATCH_SIZE]
            ids_str = ",".join(id_batch)
            try:
                logger.debug(f"Fetching summaries for {len(id_batch)} ClinVar IDs.")
                # Wait for rate limit before making request
                entrez_limiter.wait()
                fetch_handle = Entrez.esummary(db="clinvar", id=ids_str, rettype="xml")
                xml_data = fetch_handle.read()
                fetch_handle.close()
                
                # Save the API response
                track_api_queries(xml_data.decode('utf-8') if isinstance(xml_data, bytes) else xml_data)

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
                        
                        # Extract additional clinical information if available
                        review_status = docsum.findtext("review_status", "Not provided")
                        last_updated = docsum.findtext("update_date", "Unknown")
                        
                        # Check for other classifications (somatic, etc.)
                        other_classifications = []
                        for class_elem in docsum.iterfind("somatic_classification"):
                            if class_elem.findtext("description"):
                                other_classifications.append(class_elem.findtext("description"))

                        record_data.append(
                            {
                                "significance": significance or "Not provided",
                                "trait": trait_name_text,
                                "review_status": review_status,
                                "last_updated": last_updated,
                                "other_classifications": "|".join(other_classifications) if other_classifications else "None"
                            }
                        )
                    if record_data:
                        results_by_id[uid] = record_data
                        # Also update the cache
                        clinvar_cache[uid] = record_data

            except Exception as e:
                logger.error(f"Error fetching/parsing ClinVar summaries for batch: {e}")
                # Add exponential backoff on error
                time.sleep(2)
            finally:
                if "fetch_handle" in locals() and fetch_handle:
                    fetch_handle.close()
    
    # Merge newly fetched results with cached results
    merged_results = {**clinvar_cache}
    for id in clinvar_ids:
        if id in merged_results:
            continue
        if id in clinvar_cache:
            merged_results[id] = clinvar_cache[id]
    
    # Save cache after fetching new data
    if ids_to_fetch:
        save_clinvar_cache()

    return {id: clinvar_cache.get(id, []) for id in clinvar_ids}


def batch_search_ncbi_clinvar(hgvs_list: List[str]) -> Dict[str, List[Dict[str, str]]]:
    """Searches ClinVar for a list of HGVS strings using batch requests for efficiency."""
    if not hgvs_list:
        return {}

    # Create a mapping of HGVS to its hash to use as cache key
    hgvs_to_hash = {hgvs: hashlib.md5(hgvs.encode()).hexdigest() for hgvs in hgvs_list}
    
    # Check which HGVS terms are already cached
    hgvs_cache_file = os.path.join(CACHE_DIR, "hgvs_to_ids_cache.json")
    hgvs_to_ids_cache = {}
    if os.path.exists(hgvs_cache_file):
        try:
            with open(hgvs_cache_file, 'r') as f:
                hgvs_to_ids_cache = json.load(f)
            logger.info(f"Loaded {len(hgvs_to_ids_cache)} HGVS terms from cache")
        except Exception as e:
            logger.error(f"Error loading HGVS cache: {e}")
    
    # Identify which HGVS terms need to be searched
    hgvs_to_search = []
    hgvs_to_ids = {}
    
    for hgvs in hgvs_list:
        hgvs_hash = hgvs_to_hash[hgvs]
        if hgvs_hash in hgvs_to_ids_cache:
            hgvs_to_ids[hgvs] = hgvs_to_ids_cache[hgvs_hash]
        else:
            hgvs_to_search.append(hgvs)
    
    logger.info(f"Using {len(hgvs_to_ids)} cached HGVS terms, searching for {len(hgvs_to_search)} new terms")
    
    if hgvs_to_search:
        valid_hgvs_list = [h for h in hgvs_to_search if validate_hgvs(h)]
        
        # Process HGVS terms in batches using OR operators for ClinVar search
        for i in range(0, len(valid_hgvs_list), HGVS_BATCH_SIZE):
            batch = valid_hgvs_list[i:i+HGVS_BATCH_SIZE]
            
            # Create a query with all HGVS terms in the batch using OR
            # Format: "(hgvs1[Base Position]) OR (hgvs2[Base Position]) OR ..."
            query_terms = [f"({hgvs}[Base Position])" for hgvs in batch]
            batch_query = " OR ".join(query_terms)
            
            try:
                # Wait for rate limit before making request
                entrez_limiter.wait()
                logger.debug(f"Searching batch of {len(batch)} HGVS terms")
                
                # Search for all IDs matching any of the HGVS terms in this batch
                search_handle = Entrez.esearch(db="clinvar", term=batch_query, retmax=1000)
                search_results = Entrez.read(search_handle)
                search_handle.close()
                
                # Save the raw XML search response if possible
                if hasattr(search_handle, 'read') and callable(search_handle.read):
                    try:
                        search_handle.seek(0)
                        xml_data = search_handle.read()
                        track_api_queries(xml_data.decode('utf-8') if isinstance(xml_data, bytes) else xml_data)
                    except:
                        pass
                
                id_list = search_results.get("IdList", [])
                
                if not id_list:
                    logger.debug(f"No results found for batch of {len(batch)} HGVS terms")
                    for hgvs in batch:
                        hgvs_to_ids[hgvs] = []
                        hgvs_to_ids_cache[hgvs_to_hash[hgvs]] = []
                    continue
                
                logger.debug(f"Found {len(id_list)} IDs for batch of {len(batch)} HGVS terms")
                
                # For each ID found, get the details to see which HGVS it matches
                entrez_limiter.wait()
                fetch_handle = Entrez.esummary(db="clinvar", id=",".join(id_list))
                results = Entrez.read(fetch_handle)
                
                # Save the raw XML results if possible
                if hasattr(fetch_handle, 'read') and callable(fetch_handle.read):
                    try:
                        fetch_handle.seek(0)
                        xml_data = fetch_handle.read()
                        track_api_queries(xml_data.decode('utf-8') if isinstance(xml_data, bytes) else xml_data)
                    except:
                        pass
                
                fetch_handle.close()
                
                # Map each result back to the correct HGVS terms
                for result in results:
                    clinvar_id = result.get("uid", "")
                    if not clinvar_id:
                        continue
                        
                    # Look for variant information in the summary
                    title = result.get("title", "")
                    
                    # Match each HGVS term to this result
                    for hgvs in batch:
                        # If the HGVS notation is found in the title or variant description
                        # This is a simplified matching - in a production scenario, you might want more robust matching
                        hgvs_simple = hgvs.split(":")[1] if ":" in hgvs else hgvs  # Get part after colon
                        if hgvs_simple in title:
                            if hgvs not in hgvs_to_ids:
                                hgvs_to_ids[hgvs] = []
                            hgvs_to_ids[hgvs].append(clinvar_id)
                            
                            # Update cache
                            if hgvs_to_hash[hgvs] not in hgvs_to_ids_cache:
                                hgvs_to_ids_cache[hgvs_to_hash[hgvs]] = []
                            if clinvar_id not in hgvs_to_ids_cache[hgvs_to_hash[hgvs]]:
                                hgvs_to_ids_cache[hgvs_to_hash[hgvs]].append(clinvar_id)
                
                # For HGVS terms with no matches in this batch, set empty lists
                for hgvs in batch:
                    if hgvs not in hgvs_to_ids:
                        hgvs_to_ids[hgvs] = []
                        hgvs_to_ids_cache[hgvs_to_hash[hgvs]] = []
                
            except Exception as e:
                logger.error(f"Error searching ClinVar for batch: {e}")
                # Set empty results for all HGVS in this batch
                for hgvs in batch:
                    hgvs_to_ids[hgvs] = []
                    hgvs_to_ids_cache[hgvs_to_hash[hgvs]] = []
                time.sleep(2)  # Backoff on error
        
        # Save the updated cache
        try:
            with open(hgvs_cache_file, 'w') as f:
                json.dump(hgvs_to_ids_cache, f)
            logger.info(f"Saved {len(hgvs_to_ids_cache)} HGVS terms to cache")
        except Exception as e:
            logger.error(f"Error saving HGVS cache: {e}")

    # Ensure all HGVS strings in the original list have entries
    for hgvs in hgvs_list:
        if hgvs not in hgvs_to_ids:
            hgvs_to_ids[hgvs] = []

    # Collect all unique IDs that need to be fetched
    all_unique_ids = set()
    for id_list in hgvs_to_ids.values():
        all_unique_ids.update(id_list)
    
    logger.info(f"Found {len(all_unique_ids)} unique ClinVar IDs across {len(hgvs_to_ids)} HGVS terms.")

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
    """Saves results to a CSV file."""
    try:
        df = pd.DataFrame(data)
        
        # Ensure the directory exists
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        # Add metadata to the file
        df['timestamp'] = time.strftime("%Y-%m-%d %H:%M:%S")
        
        # Sort by significance to prioritize pathogenic variants
        if 'significance' in df.columns:
            # Custom sorting order for clinical significance
            significance_order = {
                'Pathogenic': 0,
                'Likely pathogenic': 1,
                'Uncertain significance': 2,
                'Likely benign': 3,
                'Benign': 4,
                'Not provided': 5,
                'Conflicting interpretations of pathogenicity': 6
            }
            # Create a sorting key using the defined order, defaulting to 999 for any undefined values
            df['_sort_key'] = df['significance'].map(lambda x: significance_order.get(x, 999))
            df = df.sort_values('_sort_key').drop('_sort_key', axis=1)
        
        df.to_csv(output_path, index=False)
        logger.info(f"Saved results to {output_path}")
    except Exception as e:
        logger.error(f"Error saving results to CSV: {e}")


def load_variants_from_json(json_path: str) -> pd.DataFrame:
    """Loads variants from a JSON file."""
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


def generate_hgvs(row):
    """Generate HGVS notation for a variant row."""
    try:
        hgvs = vcf_to_c_hgvs(
            transcript_id=row["transcript_id"],
            c_pos=row["position"],
            ref=row["reference"],
            alt=row["alternate"],
        )
        if validate_hgvs(hgvs):
            return hgvs
    except Exception as e:
        logger.debug(f"Error generating HGVS: {e}")
    return None


def process_variant(variants_file: str):
    """Processes variants from a single file by generating HGVS notations and querying ClinVar."""
    variants = load_variants_from_json(variants_file)
    if variants.empty:
        logger.warning(f"No variants loaded from {variants_file}. Skipping.")
        return []

    # Filter out rows with missing data
    valid_variants = variants.dropna(subset=["transcript_id", "position", "reference", "alternate"])
    valid_variants = valid_variants[
        valid_variants.apply(
            lambda row: all([row["transcript_id"], row["position"] > 0, row["reference"], row["alternate"]]),
            axis=1
        )
    ]
    
    if len(valid_variants) < len(variants):
        logger.warning(f"Filtered out {len(variants) - len(valid_variants)} variants with missing data")
    
    if valid_variants.empty:
        logger.warning(f"No valid variants after filtering in {variants_file}")
        return []

    logger.info(f"Generating HGVS notations for {len(valid_variants)} variants...")
    
    # Use parallel processing for HGVS generation
    with concurrent.futures.ProcessPoolExecutor(max_workers=MAX_WORKERS_HGVS) as executor:
        # Process in chunks to avoid memory issues
        chunk_size = 5000
        all_hgvs = []
        
        for i in range(0, len(valid_variants), chunk_size):
            chunk = valid_variants.iloc[i:i+chunk_size]
            futures = [executor.submit(generate_hgvs, row) for _, row in chunk.iterrows()]
            
            chunk_results = []
            for future in concurrent.futures.as_completed(futures):
                result = future.result()
                if result:
                    chunk_results.append(result)
            
            all_hgvs.extend(chunk_results)
            logger.info(f"Generated {len(chunk_results)} HGVS notations from chunk {i//chunk_size + 1}")

    unique_hgvs_list = list(set(all_hgvs))
    logger.info(f"Generated {len(unique_hgvs_list)} unique valid HGVS strings.")

    if not unique_hgvs_list:
        logger.warning("No valid HGVS strings generated. No ClinVar search performed.")
        return []

    hgvs_disease_map = batch_search_ncbi_clinvar(unique_hgvs_list)

    # Streamline saving results - only save the aggregated file
    all_results_for_file = []
    for hgvs, diseases in hgvs_disease_map.items():
        if diseases:
            for disease_info in diseases:
                # Add source file info for traceability
                result = {"hgvs": hgvs, 
                         "source_file": os.path.basename(variants_file),
                         **disease_info}
                all_results_for_file.append(result)

    # Only save individual files if there are less than 100 results
    if len(unique_hgvs_list) < 100:
        logger.info("Saving individual HGVS results...")
        
        # Create a timestamp-based subfolder to organize results
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        results_subfolder = os.path.join(PROGRAM_STORAGE_DIR_SHARED_DATA_DISEASES, 
                                        f"hgvs_results_{timestamp}")
        os.makedirs(results_subfolder, exist_ok=True)
        
        for hgvs, diseases in hgvs_disease_map.items():
            if diseases:
                # Use a safer filename by replacing any invalid characters
                safe_hgvs = re.sub(r'[\\/*?:"<>|]', "_", hgvs)
                output_path = os.path.join(results_subfolder, f"{safe_hgvs}.csv")
                save_to_csv(diseases, output_path)
    
    logger.info(f"Finished processing {variants_file} with {len(all_results_for_file)} disease associations.")
    return all_results_for_file


def load_all_variants(folder_name: str) -> List[Tuple[str, pd.DataFrame]]:
    """Load all variants from all JSON files in the folder."""
    json_files = [f for f in os.listdir(folder_name) if f.endswith(".json")]
    logger.info(f"Found {len(json_files)} JSON files in {folder_name}.")
    
    all_files_variants = []
    
    for file in json_files:
        file_path = os.path.join(folder_name, file)
        variants = load_variants_from_json(file_path)
        if not variants.empty:
            all_files_variants.append((file_path, variants))
    
    return all_files_variants


def proccess_multiple_variants(folder_name: str) -> None:
    """Process multiple variants files in a folder."""
    json_files = [f for f in os.listdir(folder_name) if f.endswith(".json")]
    logger.info(f"Found {len(json_files)} JSON files in {folder_name}.")
    
    all_disease_data_aggregated = []
    
    # Process each file
    for file in json_files:
        file_path = os.path.join(folder_name, file)
        logger.info(f"--- Processing file: {file} ---")
        disease_data_for_file = process_variant(file_path)
        all_disease_data_aggregated.extend(disease_data_for_file)
        logger.info(f"--- Finished processing file: {file} ({len(disease_data_for_file)} results) ---")

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
