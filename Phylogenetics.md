# A general overview of a common phylogenetics workflow

Alignment -> Trimming -> ML Phylogeny -> Phylodynamics

## Alignment

```{shell}
mafft --auto --thread -1 --keeplength --addfragments [input].fasta [reference].fasta > [output].fasta
```

This allows mafft to choose the alignment method and number of threads.

## Trimming

Trim ends of the alignment and any insertions that only contain ambiguous/missing bases using MEGA.

## ML Phylogeny

```{shell}
iqtree2 -s [trimmed input].fas -T AUTO -B 1000 --safe -msub viral -nt 8 -cptime 20
```

This allows iqtree2 to choose the phylogenetic model, but limits it to viral models. It runs a 1000x bootstrap on 8 threads and stores data to allow the user to resume if the job ends prematurely.

## Phylodynamics

```{shell}
treetime \
  --confidence \
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

Other options could be --greedy-resolve or --reroot

```{shell}
treetime mugration \
  --tree treetime_out/timetree.nexus \
  --states [metadata].csv \
  --attribute [country/city/etc]
```
