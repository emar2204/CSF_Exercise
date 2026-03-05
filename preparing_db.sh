# Download the human reference genome (GRCh38.p14) from NCBI:
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz

# In the directory where the file is stored, decompress the FASTA file.
# Make sure you have at least ~8 GB of free disk space available.
gunzip -c GCF_000001405.40_GRCh38.p14_genomic.fna.gz > GRCh38.p14.fna

# Once Bowtie2 is installed, build the host reference index.
bowtie2-build --threads 8 GRCh38.p14.fna GRCh38_p14_refseq_bowtie2


# MetaPhlAn and HUMAnN database setup
# 
# The latest MetaPhlAn database release may not be compatible with HUMAnN.
# If needed, use the vJun23 database release instead.

# After installing MetaPhlAn, download and install its database.
metaphlan --install --db_dir db/

# After installing HUMAnN, download and install its database.
# ChocoPhlAn nucleotide database
humann_databases --download chocophlan full "YOUR_INSTALL_LOCATION"

# UniRef90 protein database (DIAMOND format)
humann_databases --download uniref uniref90_diamond "YOUR_INSTALL_LOCATION"
