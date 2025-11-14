# Notes for Lacy

## Basecalling

[Download dorado](https://software-docs.nanoporetech.com/dorado/latest/#installation).

Start a GPU-capable slurm instance (either interactive or scripted).

```shell
dorado basecaller --kit-name SQK-NBD114-24 -o basecalled sup pod5
```

If you didn't give dorado enough time on slurm, you can use `--resume-from [basecalled.bam]` to resume your basecalling from where it stopped.

Once basecalling is finished, you have to split the BAM file into one per barcode:

```shell
cd basecalled
dorado demux --no-classify -o demuxed [basecalled.bam]
cd demuxed
```

## Genome assembly

You'll need an evironment with all the needed tools. If you don't already have one, you can use

```shell
mamba create -n assembly minimap2 samtools ivar -c conda-forge -c bioconda
```

Create a text file with all your barcode names (I'll assume it's called barcodes.txt)

```text
barcode01    MPOX_Sample1
barcode02    MPOX_Sample2
```

Then we can assemble consensus sequences and aligned BAM files:

```shell
while read N1 N2; \
  do samtools fastq *${N1}.bam > ${N2}.fastq; \
  minimap2 -ax lr:hq ref.fasta ${N2}.fastq > ${N2}.sam; \
  samtools view -bS -F 0x904 -h ${N2}.sam > ${N2}.bam; \
  samtools sort ${N2}.bam > ${N2}_sorted.bam; \
  samtools mpileup -aa -A -d 0 -Q 0 ${N2}_sorted.bam | ivar consensus -p ${N2}_con1 -q 10 -m 100 -k; \
  rm ${N2}.sam ${N2}.bam ${N2}_sorted.bam; \
  minimap2 -ax lr:hq ${N2}_con1.fa ${N2}.fastq > ${N2}.sam; \
  samtools view -bS -F 0x904 -h ${N2}.sam > ${N2}.bam; \
  samtools sort ${N2}.bam > ${N2}_sorted.bam; \
  samtools index ${N2}_sorted.bam; \
  samtools mpileup -aa -A -d 0 -Q 0 ${N2}_sorted.bam | ivar consensus -p ${N2}_consensus -q 10 -t 0.5 -n N; \
  rm ${N2}.sam ${N2}.bam ${N2}.fastq; \
done < barcodes.txt
```

## Clade assignment

We can then take all consensus sequences and begin clade assignment

```shell
mkdir phylo
cp *_consensus.fa phylo
cd phylo
```

We'll combine all the consensus sequences into one file

```shell
# Concatenate reads
for i in *_consensus.fa; do \
  base=$(basename -s _consensus.fa $i) \
  echo ">"${base} >> combined.fasta \
  tail -n +2 $i >> combined.fasta \
done
```

We'll need the nextclade software in our environment

```shell
mamba install nextclade -y
```

To get started, we need to download the nextclade hMPXV dataset

```shell
nextclade dataset get --name hMPXV --output-dir hMPXV
```

We'll then run the clade assignment and QC

```shell
nextclade run --verbose --input-dataset=hMPXV --output-tsv nextclade.tsv combined.fasta
```

Let's get rid of anything with >10% missing bases. Start by making a list of all your sequences

```shell
cat combined.fasta | grep ">" > sequences.txt
sed -i "s/>//g" sequences.txt
```

Then look at your `nextclade.tsv` and remove any sequences with >10% missing from `sequences.txt`. Once that is done, I have a script that can filter the fasta [on my github](https://github.com/seth-borrowman/Random/blob/main/filterFastaSamples.py), or you can do it manually. Save this as a new fasta called `filtered.fasta` and we'll use it for all downstream analyses.

## Phylogenetics

We'll build a tree using IQ-Tree v3.

We will need to add the tools for phylogenetics to our environment if they're not already there.

```shell
mamba install mafft iqtree=3 treetime -y
pip install clipkit
```

Align our reads based on the reference file (included in this github folder, or in the hMPXV dataset)

```shell
mafft --auto --thread -1 --keeplength --addfragments filtered.fasta MPXV-M5312_HM12_Rivers.fa > aligned_w_ref.fasta
```

We can now remove our reference sequence -- it speeds up alignment, but is not needed downstream -- and mask our alignment.

```shell
# Remove reference genome
awk 'BEGIN {keep=0} /^>/ {keep++} keep>1' aligned_w_ref.fasta > aligned.fasta
# Mask alignment
clipkit -l -o masked.fasta aligned.fasta
```

We will use iqtree to build our ML phylogeny. We'll use 1000 SH-ALRT replicates and restrict model finder to only models that make sense for viral genomes.

```shell
iqtree -s masked.fasta -T AUTO -alrt 1000 --safe -msub viral -cptime 20
```

## Methods

This should generally match what Ramon has done in the past, with a few small deviations, mostly due to software updates:

> *Genome assembly*
>
> Quality filtered nanopore reads were aligned to the reference sequence MPXV-M5312_HM12_Rivers (NC_063383.1) using minimap2 v.2.26-r1175 with the nanopore specific configuration using the ‘map-ont’ option. Samtools v1.17 was used to generate a full genome consensus sequence for each sample with a minimum alignment depth of 100 reads and minimum base quality of 10. With this initial consensus we performed a sample-specific alignment using minimap2 with the same conditions and generated an improved final consensus. Consensus sequences with ≥ 10% missing bases were discarded. For each consensus sequence hMPXV lineage assessment was performed using Nextclade with the hMPXV (all clades) reference dataset.
>
> *Phylogenetic analysis*
>
> All NM full genome consensus sequences were aligned using MAFFT v7.453 software including the NC_063383.1 reference genome. Regions with poor alignment scores were eliminated using the Gblocks tool implemented in seaview v5.0.4. Maximum Likelihood (ML) phylogenies were inferred with IQ-Tree v2.0.5 using its Model Finder function before each analysis to estimate the nucleotide substitution model best-fitted for each dataset by means of Bayesian information criterion (BIC). We assessed the tree topology for each phylogeny both with the Shimodaira–Hasegawa approximate likelihood-ratio test (SH-aLRT) with 1000 replicates.

### Some differences

- Some software versions might be slightly different - use `mamba env export` to see the version number of all installed software packages in the environment. The biggest one is IQ-Tree is now in v3 instead of v2.
- Instead of exporting our alignment to seaview, using their Gblocks function, and returning it to Quest, we used clipkit to do this all in the command line.
- We used the new nanopore specific minimap configuration `lr:hq` (long reads: high quality) instead of the original `map-ont` that assumes noisier reads. This is in line with current recommendations.
