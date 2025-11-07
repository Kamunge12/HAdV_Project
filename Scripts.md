
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

