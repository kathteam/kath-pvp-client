from Bio import Entrez

# Entrez.email = "tomkun3@ktu.lt"
# Entrez.api_key = "54d8bbd81b07bfd02ad3a1997a881f79da09"

# --- 1. Search for the EYS gene in the 'gene' database ---
search_term = "EYS[Gene]"
handle = Entrez.esearch(db="gene", term=search_term, retmax=10)
search_results = Entrez.read(handle)
handle.close()

if search_results["IdList"]:
    gene_ids = search_results["IdList"]
    print("Found Gene ID for EYS:", gene_ids[0])
    print("Found Gene ID for EYS:", len(search_results))
else:
    print("No gene found for EYS.")
    exit(1)

# --- 2. Download gene records for each gene ID ---
for gene_id in gene_ids:
    handle = Entrez.efetch(db="gene", id=gene_id, retmode="csv")
    gene_record = handle.read()
    handle.close()

    # --- 3. Save each gene record to a separate file ---
    output_filepath = f"data/EYS_gene_{gene_id}.csv"
    with open(output_filepath, "w") as f:
        f.write(gene_record)

    print(f"Gene data for ID {gene_id} saved to {output_filepath}")
