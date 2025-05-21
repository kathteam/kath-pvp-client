import sqlite3
import argparse
from pathlib import Path
import os
import logging
from xml.etree import ElementTree as ET

# Configure logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Always use this database file
BASE_OUTPUT_DIR = os.path.join(os.path.expanduser("~"), ".kath", "shared", "data", "blast_results")
DB_FILENAME = "xml.db"
SQLITE_DB_FILE = os.path.join(BASE_OUTPUT_DIR, DB_FILENAME)

def create_database(sqlite_db_file: str):
    """
    Create the SQLite database and the hsps table if it doesn't exist.
    """
    conn = sqlite3.connect(sqlite_db_file)
    cursor = conn.cursor()

    # Add file_id column
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS hsps (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        file_id INTEGER,
        hit_id INTEGER,
        hsp_num INTEGER,
        bit_score REAL,
        score INTEGER,
        evalue REAL,
        query_from INTEGER,
        query_to INTEGER,
        hit_from INTEGER,
        hit_to INTEGER,
        query_frame INTEGER,
        hit_frame INTEGER,
        identity INTEGER,
        positive INTEGER,
        gaps INTEGER,
        align_len INTEGER,
        qseq TEXT,
        hseq TEXT,
        midline TEXT
    );
    """)

    conn.commit()
    conn.close()

def get_next_file_id(sqlite_db_file: str) -> int:
    """
    Get the next file_id to use for a new XML file.
    """
    conn = sqlite3.connect(sqlite_db_file)
    cursor = conn.cursor()
    cursor.execute("SELECT MAX(file_id) FROM hsps")
    result = cursor.fetchone()
    conn.close()
    if result and result[0]:
        return result[0] + 1
    else:
        return 1

def xml_to_db(xml_file_path: str) -> str:
    """
    Parse a BLAST XML file and insert its data into a shared SQLite database.

    Args:
        xml_file_path: Path to the BLAST XML file.

    Returns:
        Path to the shared SQLite database.
    """
    # Ensure the output directory exists
    os.makedirs(BASE_OUTPUT_DIR, exist_ok=True)

    # Always use the same database file
    sqlite_db_file = SQLITE_DB_FILE

    # Create the database if needed
    create_database(sqlite_db_file)

    # Determine the file_id for this XML file
    file_id = get_next_file_id(sqlite_db_file)

    # Parse the XML file
    tree = ET.parse(xml_file_path)
    root = tree.getroot()

    # Insert data into the database
    conn = sqlite3.connect(sqlite_db_file)
    cursor = conn.cursor()

    for hsp in root.findall(".//Hsp"):
        hsp_num = int(hsp.findtext("Hsp_num"))
        bit_score = float(hsp.findtext("Hsp_bit-score"))
        score = int(hsp.findtext("Hsp_score"))
        evalue = float(hsp.findtext("Hsp_evalue"))
        query_from = int(hsp.findtext("Hsp_query-from"))
        query_to = int(hsp.findtext("Hsp_query-to"))
        hit_from = int(hsp.findtext("Hsp_hit-from"))
        hit_to = int(hsp.findtext("Hsp_hit-to"))
        query_frame = int(hsp.findtext("Hsp_query-frame"))
        hit_frame = int(hsp.findtext("Hsp_hit-frame"))
        identity = int(hsp.findtext("Hsp_identity"))
        positive = int(hsp.findtext("Hsp_positive"))
        gaps = int(hsp.findtext("Hsp_gaps"))
        align_len = int(hsp.findtext("Hsp_align-len"))
        qseq = hsp.findtext("Hsp_qseq")
        hseq = hsp.findtext("Hsp_hseq")
        midline = hsp.findtext("Hsp_midline")

        cursor.execute("""
        INSERT INTO hsps (file_id, hit_id, hsp_num, bit_score, score, evalue, query_from, query_to, hit_from, hit_to,
                          query_frame, hit_frame, identity, positive, gaps, align_len, qseq, hseq, midline)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (file_id, None, hsp_num, bit_score, score, evalue, query_from, query_to, hit_from, hit_to,
              query_frame, hit_frame, identity, positive, gaps, align_len, qseq, hseq, midline))

    conn.commit()
    conn.close()
    logger.info(f"Inserted HSP data from file_id {file_id} into the database: {sqlite_db_file}")

    return sqlite_db_file

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse BLAST XML and insert results into a shared SQLite database.")
    parser.add_argument("--xml", required=True, help="Path to BLAST XML file")
    args = parser.parse_args()

    db_path = xml_to_db(args.xml)
    print(f"Database updated at: {db_path}")