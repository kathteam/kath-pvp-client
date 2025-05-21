import sqlite3
import argparse
from pathlib import Path
import os
import logging
import json

# Configure logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Always use this database file
BASE_OUTPUT_DIR = os.path.join(os.path.expanduser("~"), ".kath", "shared", "data", "blast_results")
DB_FILENAME = "json.db"
SQLITE_DB_FILE = os.path.join(BASE_OUTPUT_DIR, DB_FILENAME)

def create_database(sqlite_db_file: str):
    """
    Create the SQLite database and the variants table if it doesn't exist.
    """
    conn = sqlite3.connect(sqlite_db_file)
    cursor = conn.cursor()
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS variants (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        file_id INTEGER,
        query_id TEXT,
        subject_id TEXT,
        chromosome TEXT,
        position INTEGER,
        variation_type TEXT,
        reference_allele TEXT,
        query_allele TEXT,
        query_position INTEGER,
        hsp_score REAL,
        hsp_evalue REAL,
        hsp_identity REAL,
        hsp_align_length INTEGER,
        hsp_query_start INTEGER,
        hsp_subject_start INTEGER,
        hsp_strand TEXT,
        hsp_gaps INTEGER
    );
    """)
    conn.commit()
    conn.close()

def get_next_file_id(sqlite_db_file: str) -> int:
    """
    Get the next file_id to use for a new JSON file.
    """
    conn = sqlite3.connect(sqlite_db_file)
    cursor = conn.cursor()
    cursor.execute("SELECT MAX(file_id) FROM variants")
    result = cursor.fetchone()
    conn.close()
    if result and result[0]:
        return result[0] + 1
    else:
        return 1

def json_to_db(json_file_path: str) -> str:
    """
    Parse a JSON file and insert its data into a shared SQLite database.

    Args:
        json_file_path: Path to the JSON file.

    Returns:
        Path to the shared SQLite database.
    """
    # Ensure the output directory exists
    os.makedirs(BASE_OUTPUT_DIR, exist_ok=True)

    # Always use the same database file
    sqlite_db_file = SQLITE_DB_FILE

    # Create the database if needed
    create_database(sqlite_db_file)

    # Determine the file_id for this JSON file
    file_id = get_next_file_id(sqlite_db_file)

    # Load the JSON file
    with open(json_file_path, "r", encoding="utf-8") as f:
        variant_list = json.load(f)

    # Insert data into the database
    conn = sqlite3.connect(sqlite_db_file)
    cursor = conn.cursor()

    for variant in variant_list:
        cursor.execute("""
        INSERT INTO variants (
            file_id, query_id, subject_id, chromosome, position, variation_type,
            reference_allele, query_allele, query_position, hsp_score, hsp_evalue,
            hsp_identity, hsp_align_length, hsp_query_start, hsp_subject_start,
            hsp_strand, hsp_gaps
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (
            file_id,
            variant.get("query_id"),
            variant.get("subject_id"),
            variant.get("chromosome"),
            variant.get("position"),
            variant.get("variation_type"),
            variant.get("reference_allele"),
            variant.get("query_allele"),
            variant.get("query_position"),
            variant.get("hsp_score"),
            variant.get("hsp_evalue"),
            variant.get("hsp_identity"),
            variant.get("hsp_align_length"),
            variant.get("hsp_query_start"),
            variant.get("hsp_subject_start"),
            json.dumps(variant.get("hsp_strand")),  # store as JSON string
            variant.get("hsp_gaps"),
        ))

    conn.commit()
    conn.close()
    logger.info(f"Inserted variant data from file_id {file_id} into the database: {sqlite_db_file}")

    return sqlite_db_file

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse JSON and insert results into a shared SQLite database.")
    parser.add_argument("--json", required=True, help="Path to JSON file")
    args = parser.parse_args()

    db_path = json_to_db(args.json)
    print(f"Database updated at: {db_path}")