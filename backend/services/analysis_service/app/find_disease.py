import pandas as pd
import os
import sys
import requests
import json
import time
from typing import Dict, Any, List, Optional, Union
from Bio import Entrez
import logging
import re
from pathlib import Path

backend_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../"))
if backend_dir not in sys.path:
    sys.path.append(backend_dir)
    
from utils.logger import get_logger
from shared.constants import PROGRAM_STORAGE_DIR_SHARED_DATA_DISEASES

os.makedirs(PROGRAM_STORAGE_DIR_SHARED_DATA_DISEASES, exist_ok=True)

logger = get_logger(__name__)


def validate_variant(chromosome: Union[int, str], position: int, reference: str, alternate: str) -> bool:
    if not re.match(r"^[0-9XYM]+$", str(chromosome)):
        logger.error(f"Invalid chromosome: {chromosome}")
        return False
    
    if not re.match(r"^[ACGT]+$", reference):
        logger.error(f"Invalid reference: {reference}")
        return False
    
    if not re.match(r"^[ACGT]+$", alternate):
        logger.error(f"Invalid alternate: {alternate}")
        return False
    
    if not isinstance(position, int):
        logger.error(f"Invalid position: {position}")
        return False
    
    return True


def query_into_myvariant(variant: str) -> Optional[Dict[str, Any]]:
    url = f"https://myvariant.info/v1/variant/{variant}"
    
    try:
        response = requests.get(url)
        response.raise_for_status()
        result = response.json()
        logger.info(f"Successfully queried MyVariant.info for {variant}")
        return result
    except requests.exceptions.HTTPError as e:
        if e.response.status_code == 404:
            logger.warning(f"Variant {variant} not found in MyVariant.info")
        else:
            logger.error(f"HTTP error querying MyVariant.info: {e}")
        return None
    except requests.exceptions.RequestException as e:
        logger.error(f"Error querying MyVariant.info: {e}")
        return None
    except json.JSONDecodeError as e:
        logger.error(f"Invalid JSON response from MyVariant.info: {e}")
        return None


def format_variant_hgvs(chromosome: str, position: int, reference: str, alternate: str) -> str:
    return f"chr{chromosome}:g.{position}{reference}>{alternate}"


def extract_gene_id(variant: Dict[str, Any]) -> Optional[str]:
    try:
        # Try dbSNP source
        if "dbsnp" in variant and "gene" in variant["dbsnp"]:
            if isinstance(variant["dbsnp"]["gene"], dict) and "geneid" in variant["dbsnp"]["gene"]:
                gene_id = variant["dbsnp"]["gene"]["geneid"]
                logger.info(f"Found gene ID {gene_id} from dbSNP")
                return gene_id
        
        # Try CADD annotation
        if "cadd" in variant and "gene" in variant["cadd"]:
            gene_data = variant["cadd"]["gene"]
            if isinstance(gene_data, dict):
                if "gene_id" in gene_data:
                    gene_id = gene_data["gene_id"]
                    if gene_id.startswith("ENSG"):
                        entrez_id = _convert_ensembl_to_entrez(gene_id)
                        if entrez_id:
                            logger.info(f"Converted Ensembl ID {gene_id} to Entrez ID {entrez_id}")
                            return entrez_id
                    else:
                        logger.info(f"Found gene ID {gene_id} from CADD")
                        return gene_id
                
                # Try gene name from CADD
                if "genename" in gene_data:
                    gene_name = gene_data["genename"]
                    entrez_id = _lookup_gene_id_by_symbol(gene_name)
                    if entrez_id:
                        logger.info(f"Found Entrez ID {entrez_id} for gene name {gene_name} from CADD")
                        return entrez_id
        
        # Try snpEff annotation
        if "snpeff" in variant and "ann" in variant["snpeff"]:
            annotations = variant["snpeff"]["ann"]
            if isinstance(annotations, list) and annotations:
                for ann in annotations:
                    if not isinstance(ann, dict):
                        continue
                    
                    # Try gene_id
                    gene_symbol = ann.get("gene_id")
                    if gene_symbol:
                        entrez_id = _lookup_gene_id_by_symbol(gene_symbol)
                        if entrez_id:
                            logger.info(f"Found Entrez ID {entrez_id} for gene symbol {gene_symbol} from snpEff")
                            return entrez_id
                    
                    # Try genename
                    gene_name = ann.get("genename")
                    if gene_name:
                        entrez_id = _lookup_gene_id_by_symbol(gene_name)
                        if entrez_id:
                            logger.info(f"Found Entrez ID {entrez_id} for gene name {gene_name} from snpEff")
                            return entrez_id
        
        # Try gene field if present
        if "gene" in variant:
            gene_info = variant["gene"]
            if isinstance(gene_info, dict):
                gene_id = gene_info.get("geneid")
                if gene_id:
                    logger.info(f"Found gene ID {gene_id} from gene field")
                    return gene_id
                
                symbol = gene_info.get("symbol")
                if symbol:
                    entrez_id = _lookup_gene_id_by_symbol(symbol)
                    if entrez_id:
                        logger.info(f"Found Entrez ID {entrez_id} for symbol {symbol} from gene field")
                        return entrez_id
            elif isinstance(gene_info, list) and gene_info:
                for gene in gene_info:
                    if isinstance(gene, dict):
                        gene_id = gene.get("geneid")
                        if gene_id:
                            logger.info(f"Found gene ID {gene_id} from gene list")
                            return gene_id
                        
                        symbol = gene.get("symbol")
                        if symbol:
                            entrez_id = _lookup_gene_id_by_symbol(symbol)
                            if entrez_id:
                                logger.info(f"Found Entrez ID {entrez_id} for symbol {symbol} from gene list")
                                return entrez_id
        
        logger.warning(f"Gene ID not found in variant data")
        return None
    except Exception as e:
        logger.error(f"Error extracting gene ID: {e}")
        return None


def _convert_ensembl_to_entrez(ensembl_id: str) -> Optional[str]:
    try:
        handle = Entrez.esearch(db="gene", term=f"{ensembl_id}[Ensembl ID] AND human[Organism]", retmode="xml")
        record = Entrez.read(handle)
        handle.close()
        
        id_list = record.get("IdList", [])
        return id_list[0] if id_list else None
    except Exception as e:
        logger.error(f"Error converting Ensembl ID to Entrez ID: {e}")
        return None


def _lookup_gene_id_by_symbol(gene_symbol: str) -> Optional[str]:
    try:
        handle = Entrez.esearch(db="gene", term=f"{gene_symbol}[Gene Symbol] AND human[Organism]", retmode="xml")
        record = Entrez.read(handle)
        handle.close()
        
        id_list = record.get("IdList", [])
        return id_list[0] if id_list else None
    except Exception as e:
        logger.error(f"Error looking up gene ID by symbol: {e}")
        return None


def fetch_clinvar_ids_by_gene(gene_id: str) -> List[str]:
    try:
        id_handler = Entrez.esearch(
            db="clinvar", 
            term=f"{gene_id}[geneID] AND (pathogenic[Clinical Significance] OR likely_pathogenic[Clinical Significance])",
            retmode="xml", 
            retmax=100
        )
        
        content = Entrez.read(id_handler)
        id_handler.close()
        
        id_list = content.get("IdList", [])
        logger.info(f"Found {len(id_list)} ClinVar entries for gene ID {gene_id}")
        return id_list
    except Exception as e:
        logger.error(f"Error fetching ClinVar IDs for gene {gene_id}: {e}")
        return []


def fetch_disease_data(clinvar_ids: List[str]) -> List[Dict[str, Any]]:
    results = []
    
    for id in clinvar_ids:
        try:
            logger.info(f"Fetching ClinVar record {id}")
            record = Entrez.esummary(db="clinvar", id=id, retmode="xml")
            content = Entrez.read(record)
            record.close()
            
            document_summaries = content["DocumentSummarySet"]["DocumentSummary"]
            
            for doc_summary in document_summaries:
                if "germline_classification" in doc_summary:
                    germline = doc_summary["germline_classification"]
                    significance = germline.get("description", "")
                    
                    if significance.lower() not in ["pathogenic", "likely pathogenic"]:
                        continue
                    
                    if "trait_set" in germline:
                        trait_sets = germline["trait_set"]
                        
                        for trait_set in trait_sets:
                            if "trait_name" in trait_set:
                                disease_name = trait_set["trait_name"]
                                if disease_name.lower() in ["not specified", "not provided"]:
                                    continue
                                
                                # Fix for potential list issue with genes
                                genes_str = ""
                                if "genes" in doc_summary and isinstance(doc_summary["genes"], list):
                                    gene_symbols = []
                                    for gene in doc_summary["genes"]:
                                        if isinstance(gene, dict) and "symbol" in gene:
                                            gene_symbols.append(gene["symbol"])
                                    genes_str = ", ".join(gene_symbols)
                                
                                disease_info = {
                                    "disease_name": disease_name,
                                    "significance": significance,
                                    "variant_id": id,
                                    "last_evaluated": germline.get("last_evaluated", ""),
                                    "review_status": germline.get("review_status", ""),
                                    "genes": genes_str
                                }
                                
                                results.append(disease_info)
                                logger.info(f"Found disease: {disease_name} ({significance})")
        
        except Exception as e:
            logger.error(f"Error processing ClinVar record {id}: {e}")
            continue
    
    return results


def save_to_csv(data: List[Dict[str, Any]], output_path: str):
    try:
        df = pd.DataFrame(data)
        df.to_csv(output_path, index=False)
        logger.info(f"Saved results to {output_path}")
    except Exception as e:
        logger.error(f"Error saving results to CSV: {e}")


def load_variants_from_csv(csv_path: str) -> pd.DataFrame:
    try:
        data = pd.read_csv(csv_path)
        data = data.dropna(subset=["chromosome", "position", "reference", "query"])
        data = data.drop_duplicates(subset=["chromosome", "position", "reference", "query"])
        
        # Ensure position is int
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


def process_variants():
    variants_file = "backend/shared/all_mutations_20250325-215708.csv"
    variants = load_variants_from_csv(variants_file)
    
    if variants.empty:
        logger.error("No valid variants to process")
        return
    
    results_count = 0
    
    for index, variant in variants.iterrows():
        try:
            # Get variant information
            chromosome = variant["chromosome"]
            position = variant["position"]
            reference = variant["reference"]
            alternate = variant["alternate"]
            
            logger.info(f"Processing variant {index+1}/{len(variants)}: chr{chromosome}:{position}{reference}>{alternate}")
            
            # Validate variant data
            if not validate_variant(chromosome, position, reference, alternate):
                logger.warning(f"Invalid variant data, skipping")
                continue
            
            # Format HGVS identifier
            hgvs_id = format_variant_hgvs(chromosome, position, reference, alternate)
            
            # Query MyVariant.info
            variant_data = query_into_myvariant(hgvs_id)
            if not variant_data:
                logger.warning(f"No data found for {hgvs_id}, skipping")
                continue
            
            # Extract gene ID with improved error handling
            gene_id = extract_gene_id(variant_data)
            if not gene_id:
                logger.warning(f"No gene ID found for {hgvs_id}, skipping")
                continue
            
            # Fetch ClinVar IDs for this gene
            clinvar_ids = fetch_clinvar_ids_by_gene(gene_id)
            if not clinvar_ids:
                logger.warning(f"No ClinVar entries found for gene ID {gene_id}, skipping")
                continue
            
            # Get disease information
            diseases = fetch_disease_data(clinvar_ids)
            if not diseases:
                logger.warning(f"No disease information found for gene ID {gene_id}, skipping")
                continue
            
            # Create safe filename
            safe_hgvs = hgvs_id.replace(":", "_").replace(">", "_to_")
            output_file = os.path.join(PROGRAM_STORAGE_DIR_SHARED_DATA_DISEASES, f"diseases_{safe_hgvs}.csv")
            
            # Save results
            save_to_csv(diseases, output_file)
            results_count += 1
        
        except Exception as e:
            logger.error(f"Error processing variant at index {index}: {e}")
            continue
            
        # Add small delay to avoid rate limits
        if index < len(variants) - 1:
            time.sleep(1)
    
    logger.info(f"Processing complete. Found disease associations for {results_count} variants.")


def main():
    # Configure NCBI Entrez
    Entrez.email = "kajeliukasc@gmail.com"
    Entrez.api_key = "7ca5ef526507701d64f16a090124cbc4aa08"
    
    # Start processing
    process_variants()


if __name__ == '__main__':
    main()