import urllib.request
import gzip
import shutil
import subprocess
import os

# URLs of the files
fa_gz_url = "https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.6.fa.gz"
bam_url = "https://ftp.ensembl.org/pub/release-112/bamcov/homo_sapiens/genebuild/GRCh38.illumina.blood.1.bam"

# Filenames
fa_gz_file = "Homo_sapiens.GRCh38.dna.chromosome.6.fa.gz"
fa_file = "Homo_sapiens.GRCh38.dna.chromosome.6.fa"
bam_file = "GRCh38.illumina.blood.1.bam"
vcf_gz_file = "chr6_blood.vcf.gz"

# Function to download files
def download_file(url, filename):
    print(f"Downloading {filename}...")
    urllib.request.urlretrieve(url, filename)
    print(f"Downloaded {filename}")

# Function to extract .gz files
def extract_gzip(gz_filename, output_filename):
    print(f"Extracting {gz_filename}...")
    with gzip.open(gz_filename, 'rb') as f_in:
        with open(output_filename, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    print(f"Extracted to {output_filename}")

# Download the files
download_file(fa_gz_url, fa_gz_file)
download_file(bam_url, bam_file)

# Update package list and install necessary tools
subprocess.run(["sudo", "apt", "update"], check=True)
subprocess.run(["sudo", "apt", "install", "bcftools", "samtools", "tabix", "-y"], check=True)

# Extract the FASTA file
extract_gzip(fa_gz_file, fa_file)

# Run bcftools mpileup and call variants
subprocess.run([
    "bcftools", "mpileup", "-f", fa_file, bam_file
], stdout=subprocess.PIPE, check=True)

subprocess.run([
    "bcftools", "call", "-mv", "-Oz", "-o", vcf_gz_file
], check=True)

# Re-extract the FASTA file again as in the original script
extract_gzip(fa_gz_file, fa_file)

