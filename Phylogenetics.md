# A general overview of a common phylogenetics workflow

Alignment -> Trimming -> ML Phylogeny -> Phylodynamics

The starting point of this workflow should be a fasta file with consensus sequences from different samples. If your consensus sequences are in separate files, but the same folder, you can combine them using something like:
```{shell}
# Concatenate reads
for i in *_Consensus.fasta; do \
  base=$(basename -s _Consensus.fasta $i) \
  echo ">"${base} >> combined.fasta \
  tail -n +2 $i >> combined.fasta \
done
```
In this example, if your files are named Sample1_Consensus.fasta, Sample2_Consensus.fasta, etc, you will have returned a file called combined.fasta that looks like
```{text}
>Sample1
ACGT
>Sample2
ACGT
```
This can be used as your input in the following Alignment step.

## Work Environment

Before moving on, we'll need to ensure some software is available for you to use. I prefer to use conda/mamba for this, but installing software locally on you computer works as well. The software that's needed includes

- mafft
- iqtree
- treetime
- clipkit
- fastplong

To install this as a conda environment, you can use

```{shell}
# Create environment
mamba create -n phylogenetics mafft iqtree treetime clipkit fastplong -c bioconda
```

This creates an environment called phylogenetics that can be activated using mamba or conda. If you prefer to use conda instead of mamba, simply replace any occurence of `mamba` with `conda`.

```{shell}
# Activate environment
mamba activate phylogenetics
```

## Quality filter

An optional step is to apply QC filtering and adapter trimming to your reads. Here, we'll use fastplong to do this. If you're using short-reads (eg. Illumina) instead of long reads (eg. ONT), you can use fastp instead - just be sure to install it into your environment `mamba install fastp`

```{shell}
# Quality filter
fastplong -i [input].fasta -o [output].fasta
```

There are some useful options for fastplong as well, inlcuding `--thread`, `--qualified_quality_phred`, and `--html` for an HTML report.

If you are using amplicon-based sequencing, you can specify length minimums and maximums as desired with `--length-required [min-length]` and `--length_limit [max-length]`

If you're unsure of what to use for different settings, it's best to leave them as the default.

## Alignment

Next we can align all our reads in a Multiple Sequence Alignment (MSA), using a reference genome to speed things up.

```{shell}
# Align reads to a referene genome
mafft --auto --thread -1 --keeplength --addfragments [input].fasta [reference].fasta > [output].fasta
```

This allows mafft to choose the alignment method and number of threads.

## Trimming

We now need to remove the reference genome if we don't want it in our phylogenetic tree. This can be done using awk

```{shell}
# Remove reference genome
awk 'BEGIN {keep=0} /^>/ {keep++} keep>1' [input].fasta > [output].fasta
```

Mask the alignment to remove sites that are not useful. To do this, we'll use clipkit, though you could use MEGA, or other software as well.

```{shell}
# Mask alignment
clipkit -l -o [output].fasta [input].fasta
```

## ML Phylogeny

Build a maximum likelihood phylogenetic tree. The current version of iqtree is iqtree3, though the command should work with iqtree2 as well.

```{shell}
# ML Tree
iqtree -s [trimmed input].fasta -T AUTO -alrt 1000 --safe -msub viral -cptime 20
```

This allows iqtree to choose the phylogenetic model, but limits it to viral models. It runs a 1000x bootstrap on an automatically decided number of threads and stores data to allow the user to resume if the job ends prematurely.

If the best model has already been found, it can be specified like `-m GTR+F+I+R5` in place of `-msub viral`

## Phylodynamics

```{shell}
# Phylodynamics
treetime \
  --confidence \
  --covariation \
  --relax 1.0 0.5 \
  --aln [aligned fasta].fasta \
  --tree [ML tree].treefile \
  --dates [metadata].csv \
  --coalescent skyline \
  --branch-length-mode marginal \
  --outdir treetime_out \
  --name-column [column in csv with names] \
  --date-column [column in csv with dates]
```

Other options could be `--greedy-resolve` or `--reroot`

Options that could be added for SARS-CoV-2 to speed up analysis:
```{shell}
--clock-filter 4 --clock-rate 0.0008 --clock-std-dev 0.0004
```

```{shell}
# Phylogeography
treetime mugration \
  --tree treetime_out/timetree.nexus \
  --states [metadata].csv \
  --attribute [country/city/etc]
```
