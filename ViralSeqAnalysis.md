# Viral Sequencing Analysis
This is a general workflow process for viral sequencing analysis - taking ONT fastq files and creating a consensus sequence and an aligned/sorted bam. There are some comments provided for analyzing data from other sequencing platforms.

## Data organization and cleaning
Concatenate .fastq files from sequencing using `cat`
It's often easiest to do this from a csv or txt file that has the barcode name and sample name on the same line.
eg.
```
barcode01,SampleName1
barcode02,SampleName2
```
They can then be concatenated with
```{shell}
while read N1 N2; do cat fastq_pass/${N1}/*fastq.gz > ${N2}.fastq.gz; done < NanoporeBarcodesFormat.txt
```
It's then often best to split any SIV or HIV samples into separate folders for further analysis.

## First alignment
Align the reads to the reference genome. For SIV we typically use [SIVMM239](https://www.ncbi.nlm.nih.gov/nucleotide/M33262) and for HIV we use [HIVHXB2CG](https://www.ncbi.nlm.nih.gov/nucleotide/K03455).
```{shell}
minimap2 -ax lr:hq [PATH]/[REFERENCE].fas [READS].fastq.gz > [NAME].sam
```
This can change slightly if you are using a separate sequencing platform. For example, the following can be used for paired-end reads:
```{shell}
minimap2 -ax sr [PATH]/[REFERENCE].fas [READS_1].fastq.gz [READS_2].fastq.gz > [NAME].sam
```
See the [minimap documentation](https://github.com/lh3/minimap2/blob/master/README.md) for more details.
For short-read sequencing data, use [BWA](https://github.com/lh3/bwa?tab=readme-ov-file)

### Clean up
This will remove reads that don't map properly and/or have a length < 1200b.
```{shell}
samtools view -e 'rlen>1200' -bS -F 4 -h [NAME].sam > [NAME].bam
samtools fastq [NAME].bam > [NAME]_filtered.fastq
```
The `-e 'rlen>1200'` option is a QC filter for length. This can be shortened or removed entirely, depending on what you want the minimum fragment length to be.

## First consensus
Use [medaka](https://github.com/nanoporetech/medaka) to generate a consensus sequence. We will then use this as our reference for further analysis.
```{shell}
medaka_consensus -i [NAME]_filtered.fastq -d [PATH]/[REFERENCE].fas -o medaka -t [NTHREADS] -m r1041_e82_400bps_sup_v5.0.0 -g
```
The r1041_e82_400bps_sup_v5.0.0 model is the one our lab has been using at the time I'm writing this. You can change this based on your nanopore basecalling model. For sequencing data from other platforms, I recommend using [iVar](https://andersen-lab.github.io/ivar/html/manualpage.html) instead of medaka.

## Realign and build final consensus
Use medaka to realign our reads to the new reference, generate a consensus sequence, and create a bam file.
```{shell}
mv medaka/consensus.fasta [NAME]_firstcon.fasta
rm -r medaka
medaka_consensus -i [NAME]_filtered.fastq -d [NAME]_firstcon.fasta -o medaka -t [NTHREADS] -m r1041_e82_400bps_sup_v5.0.0
mv medaka/consensus.fasta [NAME]_Consensus.fasta
mv medaka/calls_to_draft.bam [NAME]_sorted.bam
mv medaka/calls_to_draft.bam.bai [NAME]_sorted.bam.bai
rm -r medaka
```

## All of that as a loop:
For your convenience
```{shell}
for i in *.fastq.gz; do \
  base=$(basename -s .fastq.gz $i); \
  minimap2 -ax lr:hq [PATH]/[REFERENCE].fas $i > ${base}.sam; \
  samtools view -e 'rlen>1200' -bS -F 4 -h ${base}.sam > ${base}.bam; \
  samtools fastq ${base}.bam > ${base}_filtered.fastq; \
  medaka_consensus -i ${base}_filtered.fastq -d [PATH]/[REFERENCE].fas -o medaka -t [NTHREADS] -m r1041_e82_400bps_sup_v5.0.0 -g; \
  mv medaka/consensus.fasta ${base}_firstcon.fasta; \
  rm -r medaka; \
  medaka_consensus -i ${base}_filtered.fastq -d ${base}_firstcon.fasta -o medaka -t [NTHREADS] -m r1041_e82_400bps_sup_v5.0.0; \
  mv medaka/consensus.fasta ${base}_Consensus.fasta; \
  mv medaka/calls_to_draft.bam ${base}_sorted.bam; \
  mv medaka/calls_to_draft.bam.bai ${base}_sorted.bam.bai; \
  rm -r medaka; \
  done
```

## Old method
```{shell}
for i in [###modify_accotdingly###]*.fastq.gz; \
  do base=$(basename -s .fastq.gz $i); \
  minimap2 -ax lr:hq [PATH]/[REFERENCE].fas $i > ${base}.sam; \
  samtools view -bS -F 4 ${base}.sam > ${base}.bam; \
  samtools sort -o ${base}_sorted.bam ${base}.bam; \
  samtools index ${base}_sorted.bam; done
```
