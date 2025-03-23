import os
import csv
import xml.etree.ElementTree as ET
from collections import OrderedDict

# Function to extract data from XML-like text
def extract_data(xml_text):
    try:
        root = ET.fromstring(xml_text)
        data = OrderedDict()

        # Extract key fields
        spdi_values = []
        for elem in root.iter():
            if elem.text and elem.text.strip():
                tag = elem.tag.strip().replace("{http://www.ncbi.nlm.nih.gov}", "")
                if tag == 'SPDI':
                    # Split SPDIs by delimiter (like semicolon or space) if there are multiple SPDIs in one tag
                    spdi_values.extend([spdi.strip() for spdi in elem.text.strip().split(' ')])
                elif tag in ['NAME', 'CLINICAL_SIGNIFICANCE']:  # Only keep desired fields
                    data[tag] = elem.text.strip()

        # Prepare multiple entries for each SPDI
        entries = []
        for spdi in spdi_values:
            entry = data.copy()  # Copy other data for each SPDI
            spdi_parts = spdi.split(':')
            if len(spdi_parts) == 4:
                entry['chromosome'] = spdi_parts[0].split('_')[-1].split('.')[0]  # Extract whole number only
                entry['version'] = spdi_parts[0].split('.')[-1]
                entry['position'] = spdi_parts[1]
                entry['deleted_seq'] = spdi_parts[2]
                entry['inserted_seq'] = spdi_parts[3]
            entries.append(entry)

        return entries
    except ET.ParseError:
        return None  # Skip malformed data

# Read and process the data
parsed_data = []
with open('shared/central_data.csv', 'r', encoding='utf-8') as file:
    reader = csv.reader(file)
    for line_number, row in enumerate(reader, start=1):
        combined_text = ' '.join(row)  # Add spaces when fusing cells
        data_entries = extract_data(combined_text)
        if data_entries:
            for entry in data_entries:
                entry['line_number'] = line_number  # Track original line number
                parsed_data.append(entry)

# Desired column order
final_columns = [
    'line_number', 'NAME', 'chromosome', 'version',
    'position', 'deleted_seq', 'inserted_seq', 'CLINICAL_SIGNIFICANCE'
]

# Filter out unwanted fields
filtered_data = [
    {col: entry.get(col, '') for col in final_columns}
    for entry in parsed_data
]

# Write the structured data into a new CSV
with open('shared/formated_gene_data.csv', 'w', newline='', encoding='utf-8') as file:
    writer = csv.DictWriter(file, fieldnames=final_columns)
    writer.writeheader()
    writer.writerows(filtered_data)

print("Data successfully cleaned and saved as 'formated_gene_data.csv'.")

file_path = 'shared/central_data.csv'
os.remove(file_path)