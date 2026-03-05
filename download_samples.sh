#Ensure you have installed and running SRA toolkit (and have free 10GB of storage)
SRRS="SRR5983438 SRR5983451 SRR5983484 SRR5983306 SRR5983464 SRR5983287"
for srr in $SRRS; do
 echo "Downloading file: $srr \n"
 prefetch -O SRA_files "$srr"
 fasterq-dump "YOUR_DIRECTORY/$srr" -O YOUR_DIRECTORY/fastq --split-files -e 8
 gzip YOUR_DIRECTORY/fastq/${srr}*.fastq
done
