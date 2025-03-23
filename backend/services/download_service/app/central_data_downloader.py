from Bio import Entrez
import time  # For adding delays between requests

# Entrez.email = "tomkun3@ktu.lt"
# Entrez.api_key = "54d8bbd81b07bfd02ad3a1997a881f79da09"

retmax = 500  # (Maximum) Number of results to return

# --- 1. Search for EYS gene variants ---
search_term = "human[Organism]"
handle = Entrez.esearch(db="gene", term=search_term, retmax=retmax)
search_results = Entrez.read(handle)
handle.close()

if search_results["IdList"]:
    gene_ids = search_results["IdList"]
    print(f"Found {len(gene_ids)} Gene IDs for EYS.")
else:
    print("No gene found for EYS.")
    exit(1)

# --- 2. Download variant data in batches ---
batch_size = 50  # Recommended batch size for stability

with open("shared/central_data.csv", "w") as output_file:
    for i in range(0, len(gene_ids), batch_size):
        batch_ids = gene_ids[i:i + batch_size]  # Split into manageable chunks
        print(f"Downloading batch {i // batch_size + 1}...")

        try:
            handle = Entrez.efetch(db="snp", id=",".join(batch_ids), retmode="text", rettype="flat")
            gene_record = handle.read()
            handle.close()

            output_file.write(gene_record)
            #time.sleep(1)  # Delay to avoid overwhelming NCBI servers
        except Exception as e:
            print(f"Error with batch {i // batch_size + 1}: {e}")
