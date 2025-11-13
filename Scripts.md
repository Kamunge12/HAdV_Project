
Used fastp long tool to trim the adapters from the long reads.
```
fastplong -i /home/KWTRP/gkamunge/ADENO_ENTERIC/adeno_data/barcode92_barcode92.fastq.gz -h barcode92_barcode92.html -j barcode92_barcode92.json -o barcode92_barcode92_clean.fastq.gz
```


# converting .gz files to .bgz using bgzip to have BGZF format.
```
for f in *.fasta.gz; do gunzip -c "$f" | bgzip > "${f%.gz}.bgz" && rm "$f"; done
```

# samtools indexing of reference genomes with "samtool faidx"
```                                                        
#!/bin/bash
# load samtools/1.18
module load samtools/1.18
# loops through the reference files and does samtools faidx
for f in /home/KWTRP/gkamunge/ADENO_ENTERIC/ref_adeno/*.fasta.bgz;do
    echo " indexing: $f"
    samtools faidx "$f"
    if [ $? -eq 0 ];then
       echo " succesfully indexed:$f"
    else
       echo " failed to index:$f"
    fi
done
```



# gzipping fasta.gz refernce files
```
for f in *.fasta.gz; do gunzip -c "$f" | bgzip > temp && mv temp "$f"; done
```
# samtool faidex reference files
```
#!/bin/bash
# load samtools/1.18
module load samtools/1.18
# loops through the reference files and does samtools faidx
for f in /home/KWTRP/gkamunge/ADENO_ENTERIC/reference_files/faidex/*.fasta.gz;do
    echo " indexing: $f"
    samtools faidx "$f"
    if [ $? -eq 0 ];then
       echo " succesfully indexed:$f"
    else
       echo " failed to index:$f"
    fi
done
```
# bwa indexing refernce files
```
#!/bin/bash
# loops through the reference files and does bwa index
for f in /home/KWTRP/gkamunge/ADENO_ENTERIC/reference_files/*fasta.gz;do
    echo " indexing: $f"
    bwa index "$f"
    if [ $? -eq 0 ];then
       echo " succesfully indexed:$f"
    else
       echo " failed to index:$f"
    fi
done
```
# mapping of the reads to refernce genomes to get bam files
```
#!/bin/bash

#load samtools
module load samtools/1.18

#number of threads
threads=8

#output directoy
mkdir -p /home/KWTRP/gkamunge/ADENO_ENTERIC/mapped

#define paths
output_dir="/home/KWTRP/gkamunge/ADENO_ENTERIC/mapped"
adeno_reads="/home/KWTRP/gkamunge/ADENO_ENTERIC/adeno_data"
ref_dir="/home/KWTRP/gkamunge/ADENO_ENTERIC/reference_files"

#Process FASTQ files
for fastq_file in "$adeno_reads"/*.fastq.gz; do
    #Geting filename without extension
    filename=$(basename "$fastq_file" .fastq.gz)

    #Processing reference genomes
    for ref_file in "$ref_dir"/*.fasta.gz; do
        #Get reference filename without extension
        refname=$(basename "$ref_file" .fasta.gz)

        #Create output name-bam
        output_bam="${filename}_vs_${refname}.bam"

        echo "Aligning $filename to $refname..."

        # Runing alignment
        bwa mem "$ref_file" "$fastq_file" | samtools view -b -o "$output_dir/$output_bam" -

        echo "Saved as: $output_bam"
    done
```
# sorting and indexing bam files
```
module load samtools/1.18

input_dir="/home/KWTRP/gkamunge/ADENO_ENTERIC/mapped"
output_dir="/home/KWTRP/gkamunge/ADENO_ENTERIC/mapped/mapped_sorted"

mkdir -p "$output_dir"

for bam in "$input_dir"/*.bam
do
  name=$(basename "$bam" .bam)
  echo "Sorting: $bam"
  samtools sort -@ 8 "$bam" -o "$output_dir/${name}_sorted.bam"
  if [ $? -eq 0 ]; then
        echo "Indexing: ${name}_sorted.bam"
        samtools index "$output_dir/${name}_sorted.bam"
        echo "Done: ${name}_sorted.bam"
    else
        echo "Error sorting $bam"
        exit 1
    fi
done
```
# checking if the sorted files are corrupted and okay
```
for bam in /home/KWTRP/gkamunge/ADENO_ENTERIC/mapped/mapped_sorted/*.bam; do
    samtools quickcheck -v "$bam" || echo "Corrupted: $bam"
done
```
# Summary of the bam files
```
for bam in /home/KWTRP/gkamunge/ADENO_ENTERIC/mapped/mapped_sorted/*.bam; do     echo " Stats for $bam";     samtools flagstat "$bam" > "${bam%.bam}_flagstat.txt"; done
```
# coverage statistics of the bam files
```
for bam in *.bam; do
    echo "Coverage for $bam"
    samtools depth "$bam" | awk '{sum+=$3} END {print "Average coverage =", sum/NR}' > "${bam%.bam}_coverage.txt"
done
```
# Mapping stats summary-csv
```
#!/bin/bash

output="mapping_summary.csv"
echo "sample,total_reads,mapped_reads,mapped_percent,properly_paired,avg_coverage" > "$output"

for flagstat in *_flagstat.txt; do
    sample=${flagstat%_flagstat.txt}
    coverage_file="${sample}_coverage.txt"

    # Extract data from flagstat
    total_reads=$(grep "in total" "$flagstat" | awk '{print $1}')
    mapped_reads=$(grep "mapped (" "$flagstat" | head -n1 | awk '{print $1}')
    mapped_percent=$(grep "mapped (" "$flagstat" | head -n1 | awk -F'[()%]' '{print $2}')
    properly_paired=$(grep "properly paired" "$flagstat" | awk '{print $1}')

    # Extract average coverage value
    if [ -f "$coverage_file" ]; then
        avg_cov=$(awk '{print $4}' "$coverage_file")
    else
        avg_cov="NA"
    fi

    echo "${sample},${total_reads},${mapped_reads},${mapped_percent},${properly_paired},${avg_cov}" >> "$output"
done

echo " Summary written to: $output"
```
# mapping with minimap2
```
#!/bin/bash

#load samtools
module load samtools/1.18
# module load minimap2  # Uncomment and adjust version if needed

# Number of threads
threads=8

#output directory
mkdir -p /home/KWTRP/gkamunge/ADENO_ENTERIC/mapped_minimap

# Define paths
output_dir="/home/KWTRP/gkamunge/ADENO_ENTERIC/mapped_minimap"
adeno_reads="/home/KWTRP/gkamunge/ADENO_ENTERIC/adeno_data"
ref_dir="/home/KWTRP/gkamunge/ADENO_ENTERIC/mapped_minimap/reference"

#processing the files fastq ones
for fastq_file in "$adeno_reads"/*.fastq.gz; do
    # Getting filename without extension
    filename=$(basename "$fastq_file" .fastq.gz)

    #Processing reference genomes
    for ref_file in "$ref_dir"/*.fasta.gz; do
        # Getting reference filename without extension
        refname=$(basename "$ref_file" .fasta.gz)

        # Creating output name BAM
        output_bam="${filename}_vs_${refname}.bam"

        echo "Aligning $filename to $refname "

        #running alignment
        minimap2 -ax map-ont -t $threads "$ref_file" "$fastq_file" | samtools view -b -o "$output_dir/$output_bam" -
        
        if [ $? -eq 0 ]; then
            echo "Saved as: $output_bam"
        else
            echo "Alignment failed for $filename vs $refname"
        fi
    done
done

echo "Mapping complete."
```
# sorting and indexing bam files after mapping with minimap2
```
#!bin/bash

module load samtools/1.18
#threads
threads=8

input_dir="/home/KWTRP/gkamunge/ADENO_ENTERIC/mapped_minimap"
output_dir="/home/KWTRP/gkamunge/ADENO_ENTERIC/mapped_minimap/sorted"

mkdir -p "$output_dir"

for bam in "$input_dir"/*.bam
do
  name=$(basename "$bam" .bam)
  echo "Sorting: $bam"
  samtools sort -@ 8 "$bam" -o "$output_dir/${name}_sorted.bam"
  if [ $? -eq 0 ]; then
        echo "Indexing: ${name}_sorted.bam"
        samtools index "$output_dir/${name}_sorted.bam"
        echo "Done: ${name}_sorted.bam"
    else
        echo "Error sorting $bam"
        exit 1
    fi
done
```




