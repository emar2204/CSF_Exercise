IN_DIR="${IN_DIR:-/metagenomic_exercise/data/fastq}"
OUT_DIR="${OUT_DIR:-/metagenomic_exercise/results}"
DB_DIR="${DB_DIR:-/metagenomic_exercise/db}"
THREADS="${THREADS:-14}"

HOST_INDEX="${HOST_INDEX:-/metagenomic_exercise/data/host/GRCh38_p14_refseq_bowtie2}"
MPA_DB="${MPA_DB:-$DB_DIR/metaphlan_db}"
MPA_THREADS="${MPA_THREADS:-12}"
MPA_INDEX="mpa_vJan25_CHOCOPhlAnSGB_202503"
#MPA_INDEX="mpa_vJun23_CHOCOPhlAnSGB_202403"
HUMANN_DB_NUC="${HUMANN_DB_NUC:-$DB_DIR/humann/chocophlan}"
HUMANN_DB_PROT="${HUMANN_DB_PROT:-$DB_DIR/humann/uniref90}"
HUMANN_THREADS="${HUMANN_THREADS:-$THREADS}"
mkdir -p "$OUT_DIR"/{fastqc,trimmed,hostremoval,metaphlan,humann,tmp,logs}

echo "Input : $IN_DIR"
echo "Output: $OUT_DIR"
echo "DB    : $DB_DIR"
echo "Host  : $HOST_INDEX"
echo "Threads: $THREADS"
echo

#QUALITY CONTROL - FASTQC
echo "FastQC..."
fastqc -t "$THREADS" "$IN_DIR"/*.fastq.gz -o "$OUT_DIR/fastqc"


for r1 in "$IN_DIR"/*_1.fastq.gz; do
  srr="$(basename "$r1" _1.fastq.gz)"
  r2="$IN_DIR/${srr}_2.fastq.gz"
  echo " - $srr:"

  t1="$OUT_DIR/trimmed/${srr}_1.trim.fastq.gz"
  t2="$OUT_DIR/trimmed/${srr}_2.trim.fastq.gz"

  #READS TRIMMING
  echo "    Trimming reads..."
  fastp \
      -i "$r1" -I "$r2" \
      -o "$t1" -O "$t2" \
      --detect_adapter_for_pe \
      --length_required 50 \
      -w "$THREADS" \
      -h "$OUT_DIR/trimmed/${srr}.fastp.html" \
      -j "$OUT_DIR/trimmed/${srr}.fastp.json"\
      > "$OUT_DIR/logs/${srr}.fastp.log" 2>&1
  echo "    Host removal..."

  #HOST REMOVAL 
  bowtie2 -x "$HOST_INDEX" \
      -1 "$t1" -2 "$t2" \
      -p "$THREADS" \
      --very-sensitive \
      --un-conc-gz "$OUT_DIR/hostremoval/${srr}.clean_%.fastq.gz" \
      -S /dev/null \
       2> "$OUT_DIR/logs/${srr}.bowtie2.log"
 
  c1="$OUT_DIR/hostremoval/${srr}.clean_1.fastq.gz"
  c2="$OUT_DIR/hostremoval/${srr}.clean_2.fastq.gz"

  echo "    Taxonomic profiling..."

  #Taxonomic profiling
  metaphlan "$c1","$c2" \
    --input_type fastq \
    --db_dir "$MPA_DB" \
    --index "$MPA_INDEX" \
    --nproc "$MPA_THREADS" \
    --mapout "$OUT_DIR/metaphlan/${srr}.mapout.bz2" \
    -o "$OUT_DIR/metaphlan/${srr}.profile.tsv" \
    > "$OUT_DIR/logs/${srr}.metaphlan.log" 2>&1
  echo "   execution ended!"

  concat="$OUT_DIR/hostremoval/${srr}.clean.concat.fastq.gz"
  mpa_profile="$OUT_DIR/metaphlan/${srr}.profile.tsv"

  zcat "$c1" "$c2" | gzip -c > "$concat"
  
  #Functional profiling
  humann \
    --input "$concat" \
    --output "$OUT_DIR/humann/${srr}" \
    --threads "$HUMANN_THREADS" \
    --nucleotide-database "$HUMANN_DB_NUC" \
    --protein-database "$HUMANN_DB_PROT" \
    --taxonomic-profile "$mpa_profile" \
    > "$OUT_DIR/logs/${srr}.humann.log" 2>&1
  echo "   execution ended!"
done

mkdir -p "$OUT_DIR/humann/merged_dir"

#Copy to merged folder + utility 
merge_metaphlan_tables.py "$OUT_DIR"/metaphlan/*.profile.tsv "$OUT_DIR/metaphlan/merged_abundance_table.tsv"
find "$OUT_DIR/humann/" -type f -name "*pathabundance.tsv" -exec cp -p {} "$OUT_DIR/humann/merged_dir" \;
find "$OUT_DIR/humann/" -type f -name "*pathcoverage.tsv"  -exec cp -p {} "$OUT_DIR/humann/merged_dir" \;
humann_join_tables --input "$OUT_DIR/humann/merged_dir" --output "$OUT_DIR/humann/merged_dir/pathabundance_merged.tsv" --file_name "pathabundance"
humann_join_tables --input "$OUT_DIR/humann/merged_dir" --output "$OUT_DIR/humann/merged_dir/pathcoverage_merged.tsv"  --file_name "pathcoverage"
humann_renorm_table --input "$OUT_DIR/humann/merged_dir/pathabundance_merged.tsv" --output "$OUT_DIR/humann/merged_dir/pathabundance_merged_cpm.tsv" --units cpm --update-snames

echo
echo "Pipeline completed successfully."
echo
echo "Downstream plots are generated separately by running the two R scripts."
echo "Before running them, place the output files in the downstream_analysis folders as follows:"
echo
echo "  MetaPhlAn merged file  -> downstream_analysis/metaphlan/"
echo "  HUMAnN merged file     -> downstream_analysis/humann/"
echo
echo "Expected files:"
echo "  downstream_analysis/metaphlan/merged_abundance_table.tsv"
echo "  downstream_analysis/humann/pathabundance_merged_cpm.tsv"
echo
echo "Then run the R scripts from inside downstream_analysis (or with matching relative paths)."
