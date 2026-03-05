# CSF_Test
Bioinformatics exercise
# Shotgun metagenomics exercise - Crohn’s vs Healthy  

This repository contains a small, reproducible pipeline for shotgun metagenomics processing and profiling:
- **QC:** FastQC  
- **Read trimming:** fastp  
- **Host decontamination:** Bowtie2
- **Taxonomic profiling:** MetaPhlAn  
- **Functional profiling:** HUMAnN  
- **Merging results:** `merge_metaphlan_tables.py`, `humann_join_tables`, `humann_renorm_table`
- **Downstream analysis:** R script `taxonomic_analysis.R` and `functional_analysis.R`
> **Exercise samples (as provided in the assignment):**  
> Healthy: `SRR5983438`, `SRR5983451`, `SRR5983484`  
> Crohn: `SRR5983306`, `SRR5983464`, `SRR5983287`

---

## Repo contents
### Folders
- `output`: contains only the small and final files generated during shotgun metagenomics analysis for reproduce results of the assignment.
- `downstream_analysis`: contains the scripts used for analysis and interpretation of the taxonomic and functional profiling results as well as the plots.
  - taxonomic_analysis.R - Downstream analysis of MetaPhlAn merged profiles: taxonomic composition, alpha/beta diversity, differential abundance, genus-level inspection (Faecalibacterium) and plots.
  - functional_analysis.R - Downstream analysis of HUMAnN merged pathway tables: pathway filtering, clustering, differential functional profiling, contributor analysis and plots.


### Files
- `download_samples.sh` - Download samples using SRA.  
- `preparing_db.sh` - Host reference, MetaPhlAn and HUMAnN db preparation.

- `Dockerfile` - Builds an analysis image with all required tools (installed via micromamba).  
- `pipeline.sh` - End-to-end pipeline: From raw files to MetaPhlAn and HUMAnN tables.  


- `.dockerignore` to avoid committing large DB folders and outputs (UPDATE WITH YOUR PATHS).

---

## Requirements

- Docker
- Enough disk space (**at least 200GB**) for databases (MetaPhlAn and HUMAnN DBs can be **approx 100GB** + Reference Genome (GRCh38 **approx 10GB**) + temporary working files
- Paired-end FASTQs named as:
  ```
  <sample>_1.fastq.gz
  <sample>_2.fastq.gz
  ```

---

## 1) Build the Docker image

From the repo root:

```bash
docker build -t metagenomic_exercise:latest .
```

---

## 2) Prepare folders on your host machine

You will mount these into the container:

- `FASTQ_DIR`  - paired-end FASTQs  
- `RESULTS_DIR` - pipeline outputs  
- `DB_DIR`     - MetaPhlAn + HUMAnN databases  
- `HOST_DIR`   - Bowtie2 host index files (`*.bt2`)

Example layout (host machine):

```
project/
  fastq/
    SRR5983306_1.fastq.gz
    SRR5983306_2.fastq.gz
    ...
  results/
  db/
    metaphlan_db/
    humann/
      chocophlan/
      uniref90/
  host/
    GRCh38_p14_refseq_bowtie2.1.bt2
    GRCh38_p14_refseq_bowtie2.2.bt2
    ...
```

---

## 3) Host genome: build the Bowtie2 index

The pipeline expects a Bowtie2 **index prefix** (not a folder).  
Default in `pipeline.sh`:

- `HOST_INDEX=/metagenomic_exercise/data/host/GRCh38_p14_refseq_bowtie2`

### Example: GRCh38.p14

Download the genome FASTA (example source: NCBI FTP) and build the index:

```bash
gunzip -c GCF_000001405.40_GRCh38.p14_genomic.fna.gz > GRCh38.p14.fna

PREFIX=GRCh38_p14_refseq_bowtie2
bowtie2-build --threads 8 GRCh38.p14.fna "$PREFIX"
```

Put the resulting `PREFIX.*.bt2` files inside your host `HOST_DIR/` folder.

---

## 4) Databases

### MetaPhlAn DB
`pipeline.sh` uses:
- `MPA_DB=$DB_DIR/metaphlan_db`
- `MPA_INDEX=mpa_vJan25_CHOCOPhlAnSGB_202503`

Make sure your MetaPhlAn DB folder contains that index, or override `MPA_INDEX`.

### HUMAnN DBs (big)
`pipeline.sh` expects:
- nucleotide DB: `HUMANN_DB_NUC=$DB_DIR/humann/chocophlan`
- protein DB: `HUMANN_DB_PROT=$DB_DIR/humann/uniref90`

You can download HUMAnN databases with:

```bash
INSTALL_LOCATION=/metagenomic_exercise/db/humann

humann_databases --download chocophlan full "$INSTALL_LOCATION"
humann_databases --download uniref uniref90_diamond "$INSTALL_LOCATION"
```

> Tip: keep DBs mounted on the host so you don’t redownload them each run.

---

## 5) Run the pipeline

### Default run (recommended)
This matches the defaults in `pipeline.sh`:

```bash
FASTQ_DIR="/path/to/fastq"
RESULTS_DIR="/path/to/results"
DB_DIR="/path/to/db"
HOST_DIR="/path/to/host"

docker run --rm -it   -v "${FASTQ_DIR}":/metagenomic_exercise/data/fastq   -v "${RESULTS_DIR}":/metagenomic_exercise/results   -v "${DB_DIR}":/metagenomic_exercise/db   -v "${HOST_DIR}":/metagenomic_exercise/data/host   -w /metagenomic_exercise   metagenomic_exercise:latest
```

### Override threads (example)
```bash
docker run --rm -it   -e THREADS=16   -e MPA_THREADS=16   -e HUMANN_THREADS=16   -v "${FASTQ_DIR}":/metagenomic_exercise/data/fastq   -v "${RESULTS_DIR}":/metagenomic_exercise/results   -v "${DB_DIR}":/metagenomic_exercise/db   -v "${HOST_DIR}":/metagenomic_exercise/data/host   -w /metagenomic_exercise   metagenomic_exercise:latest
```

---

## Configuration (environment variables)

Defaults in `pipeline.sh`:

- `IN_DIR=/metagenomic_exercise/data/fastq`  
- `OUT_DIR=/metagenomic_exercise/results`  
- `DB_DIR=/metagenomic_exercise/db`  
- `THREADS=14`  

Host removal:
- `HOST_INDEX=/metagenomic_exercise/data/host/GRCh38_p14_refseq_bowtie2`

MetaPhlAn:
- `MPA_DB=$DB_DIR/metaphlan_db`
- `MPA_INDEX=mpa_vJan25_CHOCOPhlAnSGB_202503`
- `MPA_THREADS=12`

HUMAnN:
- `HUMANN_DB_NUC=$DB_DIR/humann/chocophlan`
- `HUMANN_DB_PROT=$DB_DIR/humann/uniref90`
- `HUMANN_THREADS=$THREADS`


---

## What the pipeline does (step-by-step)

### 0) FastQC 
FastQC call produce HTML QC reports in `results/fastqc/`.

### 1) Trimming (fastp)
For each sample pair:
- PE adapter detection
- min read length: 50
Outputs:
- `results/trimmed/<sample>_1.trim.fastq.gz`
- `results/trimmed/<sample>_2.trim.fastq.gz`
- `results/trimmed/<sample>.fastp.html`
- `results/trimmed/<sample>.fastp.json`

### 2) Host removal (Bowtie2)
Maps trimmed reads to host index, keeps **unmapped pairs**.
Outputs:
- `results/hostremoval/<sample>.clean_1.fastq.gz`
- `results/hostremoval/<sample>.clean_2.fastq.gz`

### 3) Taxonomic profiling (MetaPhlAn)
Runs MetaPhlAn on the cleaned reads.
Outputs:
- `results/metaphlan/<sample>.profile.tsv`
- `results/metaphlan/<sample>.mapout.bz2`

### 4) Functional profiling (HUMAnN)
HUMAnN is run using the MetaPhlAn taxonomic profile (`--taxonomic-profile`) in order to not recomputing them.
Outputs per sample:
- `results/humann/<sample>/` (HUMAnN tables + intermediates)
- logs: `results/logs/<sample>.humann.log`

### 5) Merge MetaPhlAn profiles into one table
- Combine all sample *.profile.tsv files generated by MetaPhlAn into one merged matrix.

### 6) Merge HUMAnN tables across samples
The pipeline collects all `*pathabundance.tsv` and `*pathcoverage.tsv`

Then creates:
- `pathabundance_merged.tsv`
- `pathcoverage_merged.tsv`
- `pathabundance_merged_cpm.tsv`

---

## Output structure

```
results/
  fastqc/         # optional (if enabled)
  trimmed/        # fastp trimmed reads + reports
  hostremoval/    # host-decontaminated reads (.clean_1/.clean_2)
  metaphlan/      # MetaPhlAn profiles + mapout
  humann/         # HUMAnN per-sample outputs + merged_dir/
  logs/           # fastp/bowtie2/metaphlan/humann logs
  tmp/            # temp folder (if used by tools)
```
---

---
# Downstream analysis files

After running the pipeline and generating the merged MetaPhlAn and HUMAnN output tables, the final comparative analyses can be performed with the two R scripts included in the repository:

## taxonomic_analysis.R:
This script is used to analyze the taxonomic profiles generated by MetaPhlAn after merging sample-level outputs.
### Includes:
- relative abundance exploration across samples and groups
- alpha diversity metrics
- beta diversity
- clustering and heatmaps
- differential abundance between Crohn’s disease and healthy controls
- focused inspection of relevant taxa such as Faecalibacterium

## functional_analysis.R
This script is used to analyze the functional pathway profiles generated by HUMAnN after table joining and normalization.
### Includes:
- pathway abundance exploration across samples and groups
- filtering of low-prevalence pathways
- heatmaps and clustering of discriminative pathways
- comparison of pathway abundance between Crohn’s disease and healthy controls
- interpretation of functional shifts in the context of taxonomic changes
- taxon-contribution analysis for selected pathways
---

