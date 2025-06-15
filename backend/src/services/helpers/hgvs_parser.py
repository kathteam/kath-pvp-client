"""
HGVS Parser Utility

This module provides utilities for parsing and validating HGVS strings
for genetic variants.
"""

import re
from typing import Dict, List, Optional, Union, Tuple

from utils import get_logger

logger = get_logger(__name__)

# Regular expression patterns for different HGVS formats
SNV_PATTERN = re.compile(r"chr([0-9XYM]+):g\.(\d+)([ACGT])>([ACGT])")  # Single nucleotide variant
DNV_MNV_PATTERN = re.compile(
    r"chr([0-9XYM]+):g\.(\d+)([ACGT]{2,})>([ACGT]{2,})"
)  # Di/Multi-nucleotide variant
DELETION_PATTERN = re.compile(r"chr([0-9XYM]+):g\.(\d+)(?:_(\d+))?del(?:[ACGT]*)?")  # Deletion
INSERTION_PATTERN = re.compile(r"chr([0-9XYM]+):g\.(\d+)_(\d+)ins([ACGT]+)")  # Insertion
DELINS_PATTERN = re.compile(r"chr([0-9XYM]+):g\.(\d+)(?:_(\d+))?delins([ACGT]+)")  # Deletion-insertion
DUPLICATION_PATTERN = re.compile(
    r"chr([0-9XYM]+):g\.(\d+)(?:_(\d+))?dup(?:[ACGT]*)?"
)  # Duplication
# Handle non-standard formats
GENOMIC_VARIANT_PATTERN = re.compile(
    r"([0-9XYM]+)[:-](\d+)[-_]?(?:(\d+))?[:-]?([ACGT]*)[:-]>?([ACGT]*)"
)


def parse_hgvs(hgvs_id: str) -> Optional[Dict[str, Union[str, int]]]:
    """
    Parse an HGVS identifier into its components

    Args:
        hgvs_id (str): The HGVS identifier string

    Returns:
        Optional[Dict[str, Union[str, int]]]: Dictionary with parsed components or None if parsing failed
    """
    if not hgvs_id:
        return None

    # Normalize the HGVS string
    hgvs_id = hgvs_id.strip()

    # Try to match with different HGVS patterns

    # Single nucleotide variant (SNV)
    match = SNV_PATTERN.match(hgvs_id)
    if match:
        chrom, pos, ref, alt = match.groups()
        return {
            "chromosome": chrom,
            "position": int(pos),
            "reference": ref,
            "alternate": alt,
            "type": "SNV",
            "hgvs": hgvs_id,
            "end_position": int(pos),
        }

    # Di/Multi-nucleotide variant (DNV/MNV)
    match = DNV_MNV_PATTERN.match(hgvs_id)
    if match:
        chrom, pos, ref, alt = match.groups()
        return {
            "chromosome": chrom,
            "position": int(pos),
            "reference": ref,
            "alternate": alt,
            "type": "MNV",
            "hgvs": hgvs_id,
            "end_position": int(pos) + len(ref) - 1,
        }
        
    # Deletion-insertion variant
    match = DELINS_PATTERN.match(hgvs_id)
    if match:
        groups = match.groups()
        chrom = groups[0]
        start = int(groups[1])
        end = int(groups[2]) if groups[2] else start
        inserted = groups[3]

        return {
            "chromosome": chrom,
            "position": start,
            "end_position": end,
            "reference": "",  # Would need sequence context
            "alternate": inserted,
            "type": "DELINS",
            "hgvs": hgvs_id,
        }

    # Deletion variant
    match = DELETION_PATTERN.match(hgvs_id)
    if match:
        groups = match.groups()
        chrom = groups[0]
        start = int(groups[1])
        end = int(groups[2]) if groups[2] else start

        return {
            "chromosome": chrom,
            "position": start,
            "end_position": end,
            "reference": "",  # Would need sequence context to know what was deleted
            "alternate": "-",  # Represent deletion as dash
            "type": "DEL",
            "hgvs": hgvs_id,
        }

    # Insertion variant
    match = INSERTION_PATTERN.match(hgvs_id)
    if match:
        chrom, pos1, pos2, inserted = match.groups()
        return {
            "chromosome": chrom,
            "position": int(pos1),
            "end_position": int(pos2),
            "reference": "-",  # Represent insertion as having no reference
            "alternate": inserted,
            "type": "INS",
            "hgvs": hgvs_id,
        }

    # Duplication variant
    match = DUPLICATION_PATTERN.match(hgvs_id)
    if match:
        groups = match.groups()
        chrom = groups[0]
        start = int(groups[1])
        end = int(groups[2]) if groups[2] else start

        return {
            "chromosome": chrom,
            "position": start,
            "end_position": end,
            "reference": "",  # Would need sequence context
            "alternate": "DUP",  # Placeholder - actual sequence would be determined with reference
            "type": "DUP",
            "hgvs": hgvs_id,
        }

    # Try for non-standard format (like basic genomic coordinates)
    match = GENOMIC_VARIANT_PATTERN.match(hgvs_id)
    if match:
        chrom, start, end, ref, alt = match.groups()
        if not end:
            end = start

        var_type = determine_variant_type(ref, alt)

        return {
            "chromosome": chrom.replace("chr", ""),
            "position": int(start),
            "end_position": int(end),
            "reference": ref if ref else "",
            "alternate": alt if alt else "",
            "type": var_type,
            "hgvs": hgvs_id,
            "normalized": False,
        }

    # Try for simpler pattern if none of the above matched
    if ":g." in hgvs_id:
        try:
            chrom = hgvs_id.split(":g.")[0].replace("chr", "")
            pos_part = hgvs_id.split(":g.")[1]

            # Try to extract position
            pos_match = re.search(r"(\d+)", pos_part)
            if pos_match:
                position = int(pos_match.group(1))
                # If only position is available (no reference or alternate), return None
                if not re.search(r"[ACGT]", pos_part):
                    return None
                return {
                    "chromosome": chrom,
                    "position": position,
                    "end_position": position,
                    "reference": "",
                    "alternate": "",
                    "type": "UNKNOWN",
                    "hgvs": hgvs_id,
                    "normalized": False,
                }
        except Exception as e:
            logger.warning(f"Failed to parse HGVS with basic method: {hgvs_id}, error: {e}")

    logger.warning(f"Could not parse HGVS: {hgvs_id}")
    return None


def determine_variant_type(ref: str, alt: str) -> str:
    """Determine variant type based on ref and alt alleles"""
    if not ref and not alt:
        return "UNKNOWN"
    elif not ref or ref == "-":
        return "INS"
    elif not alt or alt == "-":
        return "DEL"
    elif len(ref) == 1 and len(alt) == 1:
        return "SNV"
    elif len(ref) > 1 and len(alt) > 1 and len(ref) == len(alt):
        return "MNV"
    elif len(ref) > 0 and len(alt) > 0 and len(ref) != len(alt):
        return "DELINS"
    else:
        return "COMPLEX"


def format_to_hgvs(chrom: str, pos: int, ref: str, alt: str) -> str:
    """
    Format variant components into an HGVS string

    Args:
        chrom (str): Chromosome
        pos (int): Position
        ref (str): Reference allele
        alt (str): Alternate allele

    Returns:
        str: HGVS formatted string
    """
    # Format chromosome
    chrom_str = str(chrom).upper().replace("CHR", "")

    # Handle different variant types
    if ref == "-" or not ref:
        return f"chr{chrom_str}:g.{pos}_{pos+1}ins{alt}"
    elif alt == "-" or not alt:
        if len(ref) == 1:
            return f"chr{chrom_str}:g.{pos}del"
        return f"chr{chrom_str}:g.{pos}_{pos + len(ref) - 1}del"
    else:
        # Check if substitution, deletion-insertion, or MNV
        if len(ref) == 1 and len(alt) == 1:
            # SNP
            return f"chr{chrom_str}:g.{pos}{ref}>{alt}"
        elif len(ref) == len(alt):
            # MNV (multiple substitution)
            return f"chr{chrom_str}:g.{pos}{ref}>{alt}"
        else:
            # Complex event (delins)
            if len(ref) == 1:
                return f"chr{chrom_str}:g.{pos}delins{alt}"
            else:
                return f"chr{chrom_str}:g.{pos}_{pos + len(ref) - 1}delins{alt}"


def convert_genome_coordinates_to_hgvs(chrom: str, start: int, end: int, ref: str, alt: str) -> str:
    """
    Convert genomic coordinates to HGVS format

    Args:
        chrom (str): Chromosome number/letter (e.g., '1', 'X')
        start (int): Start position
        end (int): End position
        ref (str): Reference allele
        alt (str): Alternate allele

    Returns:
        str: HGVS formatted variant
    """
    # Format chromosome
    chrom_str = str(chrom).upper().replace("CHR", "")

    # Determine variant type
    if ref == "-" or not ref:
        # Insertion
        return f"chr{chrom_str}:g.{start}_{end}ins{alt}"
    elif alt == "-" or not alt:
        # Deletion
        if start == end:
            return f"chr{chrom_str}:g.{start}del"
        else:
            return f"chr{chrom_str}:g.{start}_{end}del"
    elif len(ref) == 1 and len(alt) == 1 and start == end:
        # SNV
        return f"chr{chrom_str}:g.{start}{ref}>{alt}"
    elif len(ref) > 1 and len(alt) > 1 and len(ref) == len(alt) and (end - start + 1) == len(ref):
        # MNV
        return f"chr{chrom_str}:g.{start}{ref}>{alt}"
    else:
        # Deletion-insertion
        if start == end:
            return f"chr{chrom_str}:g.{start}delins{alt}"
        else:
            return f"chr{chrom_str}:g.{start}_{end}delins{alt}"


def is_valid_hgvs(hgvs_id: str) -> bool:
    """
    Check if a string is a valid HGVS identifier

    Args:
        hgvs_id (str): The string to check

    Returns:
        bool: True if the string is a valid HGVS identifier, False otherwise
    """
    return parse_hgvs(hgvs_id) is not None


def normalize_chromosome(chrom: str) -> str:
    """
    Normalize chromosome notation (handles X, Y, MT/M and numeric chromosomes)

    Args:
        chrom (str): Chromosome identifier

    Returns:
        str: Normalized chromosome identifier
    """
    # Remove 'chr' prefix if present
    chrom = str(chrom).upper().replace("CHR", "")

    # Handle mitochondrial chromosome
    if chrom in ["M", "MT"]:
        return "MT"

    # Handle numeric chromosomes, X and Y
    if chrom in [str(i) for i in range(1, 23)] + ["X", "Y"]:
        return chrom

    # Try to extract number from string
    match = re.search(r"(\d+|X|Y|MT?)", chrom, re.IGNORECASE)
    if match:
        ch = match.group(1).upper()
        if ch == "M":
            return "MT"
        return ch

    return chrom  # Return as is if we can't normalize


def get_variant_description(variant_type: str) -> str:
    """
    Get a human-readable description of a variant type

    Args:
        variant_type (str): Variant type code (SNV, INS, DEL, etc.)

    Returns:
        str: Human-readable description
    """
    descriptions = {
        "SNV": "Single Nucleotide Variant",
        "MNV": "Multi-Nucleotide Variant",
        "INS": "Insertion",
        "DEL": "Deletion",
        "DELINS": "Deletion-Insertion",
        "DUP": "Duplication",
        "COMPLEX": "Complex Rearrangement",
        "UNKNOWN": "Unknown Variant Type",
    }
    return descriptions.get(variant_type, "Variant")


def batch_parse_hgvs(hgvs_list: List[str]) -> List[Dict[str, Union[str, int]]]:
    """
    Parse multiple HGVS identifiers at once

    Args:
        hgvs_list (List[str]): List of HGVS identifiers

    Returns:
        List[Dict[str, Union[str, int]]]: List of parsed HGVS dictionaries
    """
    results = []
    for hgvs in hgvs_list:
        parsed = parse_hgvs(hgvs)
        if parsed:
            results.append(parsed)
        else:
            logger.warning(f"Failed to parse HGVS: {hgvs}")

    return results


def vcf_to_hgvs(chrom: str, pos: int, ref: str, alt: str) -> str:
    """
    Convert VCF-style variant to HGVS format

    Args:
        chrom (str): Chromosome
        pos (int): Position (1-based)
        ref (str): Reference allele
        alt (str): Alternate allele

    Returns:
        str: HGVS formatted variant
    """
    chrom_str = normalize_chromosome(chrom)

    # Handle the different variation types
    if ref == alt:
        # Not a variant
        return None
    elif len(ref) == 1 and len(alt) == 1:
        # SNV
        return f"chr{chrom_str}:g.{pos}{ref}>{alt}"
    elif len(ref) == 0:
        # Insertion in VCF (not standard, but sometimes seen)
        return f"chr{chrom_str}:g.{pos}_{pos}ins{alt}"
    elif len(alt) == 0:
        # Deletion in VCF (not standard, but sometimes seen)
        if len(ref) == 1:
            return f"chr{chrom_str}:g.{pos}del"
        else:
            return f"chr{chrom_str}:g.{pos}_{pos + len(ref) - 1}del"
    elif len(ref) > len(alt):
        # VCF-style deletion or delins
        if alt[0] == ref[0]:  # Shares prefix
            # Delete everything after the shared prefix
            prefix = alt[0]
            deleted_bases = len(ref) - len(alt)
            if len(alt) == 1:
                # Simple deletion
                if deleted_bases == 1:
                    return f"chr{chrom_str}:g.{pos + 1}del"
                else:
                    return f"chr{chrom_str}:g.{pos + 1}_{pos + deleted_bases}del"
            else:
                # Delins
                return f"chr{chrom_str}:g.{pos + 1}_{pos + deleted_bases}delins{alt[1:]}"
        else:
            # Complete replacement
            return f"chr{chrom_str}:g.{pos}_{pos + len(ref) - 1}delins{alt}"
    elif len(alt) > len(ref):
        # VCF-style insertion or delins
        if alt[0] == ref[0]:  # Shares prefix
            # Insert after the shared prefix
            inserted_bases = alt[1:]
            return f"chr{chrom_str}:g.{pos}_{pos + 1}ins{inserted_bases}"
        else:
            # Complete replacement
            return f"chr{chrom_str}:g.{pos}_{pos + len(ref) - 1}delins{alt}"
    else:  # len(alt) == len(ref)
        # MNV
        return f"chr{chrom_str}:g.{pos}{ref}>{alt}"


def vcf_to_c_hgvs(transcript_id: str, c_pos: int, ref: str, alt: str) -> str:
    """
    Convert to coding HGVS (e.g., NM_000314.8:c.395G>T)
    """
    if len(ref) == 1 and len(alt) == 1:
        return f"{transcript_id}:c.{c_pos}{ref}>{alt}"
    elif len(ref) > 1 and len(alt) == 0:
        if len(ref) == 1:
            return f"{transcript_id}:c.{c_pos}del"
        else:
            return f"{transcript_id}:c.{c_pos}_{c_pos + len(ref) - 1}del"
    elif len(ref) == 0 and len(alt) > 0:
        return f"{transcript_id}:c.{c_pos}_{c_pos + 1}ins{alt}"
    elif len(ref) > 1 and len(alt) > 0:
        return f"{transcript_id}:c.{c_pos}_{c_pos + len(ref) - 1}delins{alt}"
    else:
        return f"{transcript_id}:c.{c_pos}{ref}>{alt}"


def validate_hgvs(hgvs_id: str) -> bool:
    """
    Validate an HGVS identifier

    Args:
        hgvs_id (str): The HGVS identifier string

    Returns:
        bool: True if valid, False otherwise
    """

    HGVS_REGEX = r"""
    (?P<accession>[A-Z]{2}_\d+\.\d+)      # Accession, e.g., NC_000023.11
    :
    (?P<type>[gcmnrp])\.                  # Type: g, c, m, n, r, p
    (?P<description>.+)                   # Description, e.g., 32867885G>A or Arg175His
"""

    # Compile the regex with verbose mode
    hgvs_regex = re.compile(HGVS_REGEX, re.VERBOSE)

    match = hgvs_regex.match(hgvs_id)
    if match:
        return True
    else:
        logger.warning(f"Invalid HGVS identifier: {hgvs_id}")
        return False
