
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
```

done
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



