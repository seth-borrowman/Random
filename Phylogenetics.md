# A general overview of a common phylogenetics workflow

Alignment -> Trimming -> ML Phylogeny -> Phylodynamics

The starting point of this workflow should be a fasta file with consensus sequences from different samples. If your consensus sequences are in separate files, but the same folder, you can combine them using something like:
```{shell}
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

## Alignment

```{shell}
mafft --auto --thread -1 --keeplength --addfragments [input].fasta [reference].fasta > [output].fasta
```

This allows mafft to choose the alignment method and number of threads.

## Trimming

Trim ends of the alignment and any insertions that only contain ambiguous/missing bases using MEGA.

## ML Phylogeny

```{shell}
iqtree2 -s [trimmed input].fas -T AUTO -alrt 1000 --safe -msub viral -cptime 20
```

This allows iqtree2 to choose the phylogenetic model, but limits it to viral models. It runs a 1000x bootstrap on an automatically decided number of threads and stores data to allow the user to resume if the job ends prematurely.

If the best model has already been found, it can be specified like `-m GTR+F+I+R5` in place of `-msub viral`

## Phylodynamics

```{shell}
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

Options to add for SARS-CoV-2:
```{shell}
--clock-filter 4 --clock-rate 0.0008 --clock-std-dev 0.0004
```

```{shell}
treetime mugration \
  --tree treetime_out/timetree.nexus \
  --states [metadata].csv \
  --attribute [country/city/etc]
```
