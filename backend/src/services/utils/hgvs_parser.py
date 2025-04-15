"""
HGVS Parser Utility

This module provides utilities for parsing and validating HGVS strings
for genetic variants.
"""

import re
import logging
from typing import Dict, Optional, Tuple, Union

logger = logging.getLogger(__name__)

# HGVS pattern matchers
# Reference: https://varnomen.hgvs.org/recommendations/

# Basic HGVS patterns
SNV_PATTERN = re.compile(r'chr([0-9XYM]+):g\.(\d+)([ACGT])>([ACGT])')  # Single nucleotide variant
DNV_MNV_PATTERN = re.compile(r'chr([0-9XYM]+):g\.(\d+)([ACGT]{2,})>([ACGT]{2,})')  # Di/Multi-nucleotide variant
DELETION_PATTERN = re.compile(r'chr([0-9XYM]+):g\.(\d+)(?:_(\d+))?del(?:[ACGT]+)?')  # Deletion
INSERTION_PATTERN = re.compile(r'chr([0-9XYM]+):g\.(\d+)_(\d+)ins([ACGT]+)')  # Insertion
DELINS_PATTERN = re.compile(r'chr([0-9XYM]+):g\.(\d+)(?:_(\d+))?delins([ACGT]+)')  # Deletion-insertion

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
            'chromosome': chrom,
            'position': int(pos),
            'reference': ref,
            'alternate': alt,
            'type': 'SNV',
            'hgvs': hgvs_id
        }
    
    # Di/Multi-nucleotide variant (DNV/MNV)
    match = DNV_MNV_PATTERN.match(hgvs_id)
    if match:
        chrom, pos, ref, alt = match.groups()
        return {
            'chromosome': chrom,
            'position': int(pos),
            'reference': ref,
            'alternate': alt,
            'type': 'MNV',
            'hgvs': hgvs_id
        }
    
    # Deletion variant
    match = DELETION_PATTERN.match(hgvs_id)
    if match:
        groups = match.groups()
        chrom = groups[0]
        start = int(groups[1])
        end = int(groups[2]) if groups[2] else start
        
        return {
            'chromosome': chrom,
            'position': start,
            'end_position': end,
            'reference': '',  # Would need sequence context to know what was deleted
            'alternate': '-',  # Represent deletion as dash
            'type': 'DEL',
            'hgvs': hgvs_id
        }
    
    # Insertion variant
    match = INSERTION_PATTERN.match(hgvs_id)
    if match:
        chrom, pos1, pos2, inserted = match.groups()
        return {
            'chromosome': chrom,
            'position': int(pos1),
            'end_position': int(pos2),
            'reference': '-',  # Represent insertion as having no reference
            'alternate': inserted,
            'type': 'INS',
            'hgvs': hgvs_id
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
            'chromosome': chrom,
            'position': start,
            'end_position': end,
            'reference': '',  # Would need sequence context
            'alternate': inserted,
            'type': 'DELINS',
            'hgvs': hgvs_id
        }
    
    # Try for simpler pattern if none of the above matched
    if ':g.' in hgvs_id:
        try:
            chrom = hgvs_id.split(':g.')[0].replace('chr', '')
            pos_part = hgvs_id.split(':g.')[1]
            
            # Try to extract position 
            pos_match = re.search(r'(\d+)', pos_part)
            if pos_match:
                position = int(pos_match.group(1))
                
                return {
                    'chromosome': chrom,
                    'position': position,
                    'reference': '',
                    'alternate': '',
                    'type': 'UNKNOWN',
                    'hgvs': hgvs_id
                }
        except Exception as e:
            logger.warning(f"Failed to parse HGVS with basic method: {hgvs_id}, error: {e}")
    
    logger.warning(f"Could not parse HGVS: {hgvs_id}")
    return None

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
        # SNP / MNP
        return f"chr{chrom_str}:g.{pos}{ref}>{alt}"

def is_valid_hgvs(hgvs_id: str) -> bool:
    """
    Check if a string is a valid HGVS identifier
    
    Args:
        hgvs_id (str): The HGVS string to validate
        
    Returns:
        bool: True if valid HGVS format, False otherwise
    """
    # Check basic format first
    if not hgvs_id or ':g.' not in hgvs_id:
        return False
        
    # Valid if we can parse it
    return parse_hgvs(hgvs_id) is not None