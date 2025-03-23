#!/bin/bash

wget https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.6.fa.gz
wget https://ftp.ensembl.org/pub/release-112/bamcov/homo_sapiens/genebuild/GRCh38.illumina.blood.1.bam
sudo apt update
sudo apt install bcftools samtools tabix -y

gunzip Homo_sapiens.GRCh38.dna.chromosome.6.fa.gz

bcftools mpileup -f Homo_sapiens.GRCh38.dna.chromosome.6.fa GRCh38.illumina.blood.1.bam | bcftools call -mv -Oz -o chr6_blood.vcf.gz

gunzip Homo_sapiens.GRCh38.dna.chromosome.6.fa.gz
