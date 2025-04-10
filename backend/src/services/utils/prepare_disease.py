import os
import pandas as pd
from typing import List, Dict, Any


def _read_disease_data(file_path: str) -> pd.DataFrame:
    """
    Read disease data from a CSV file.

    :param file_path: Path to the CSV file containing disease data.
    :return: DataFrame containing the disease data.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")

    return pd.read_csv(file_path)


def _prepare_disease_data(df: pd.DataFrame) -> List[Dict[str, Any]]:
    """
    Prepare disease data for display.

    :param df: DataFrame containing the disease data.
    :return: List of dictionaries with prepared disease data.
    """
    # Drop rows with missing values in 'Disease' and 'Description' columns
    data = df[["clinical_significance", "disease_name"]].dropna()

    json_format = []
    for _, row in data.iterrows():
        json_format.append(
            {
                "clinical_significance": row["clinical_significance"],
                "disease_name": row["disease_name"],
            }
        )

    return json_format


def generate_json_diseases(file_path: str) -> List[Dict[str, Any]]:
    """
    Main function to read and prepare disease data.

    :param file_path: Path to the CSV file containing disease data.
    :return: List of dictionaries with prepared disease data.
    """
    df = _read_disease_data(file_path)
    prepared_data = _prepare_disease_data(df)

    return prepared_data
