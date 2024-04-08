# GenomeAnnotationWorkshop2024
An introduction to Genome Annotation of non-model organisms

This is the github repo for the the Genome Annotation Workshop The workshop is focused on annotating genomes of non-model organisms using a custom pipeline of multiple tools.
In the workshop different strategies, such as homology-based, *ab initio*, and *de novo* approaches are implemented, using a combination of short and long reads (Iso-Seq) available on NCBI.
As examples, a single chromosome from three different organisms is used as a demonstration.

## Installation
To install this site locally run the following commands:

Clone the repo and cd into it 
```
git clone [git@github.com:griffithlab/rnabio.org.git](https://github.com/francicco/GenomeAnnotationWorkshop2024.git)
```
Install the following software and their dependencies:

[STAR](https://github.com/alexdobin/STAR)

[Samtools](https://github.com/samtools/samtools)

[Bedtools](https://github.com/arq5x/bedtools2)

[Diamond](https://github.com/bbuchfink/diamond)

[Miniprot](https://github.com/lh3/miniprot)

[BRAKER](https://github.com/Gaius-Augustus/BRAKER)

[BUSCO](https://busco.ezlab.org/)

[Stringtie](https://ccb.jhu.edu/software/stringtie/)

[Scallop](https://github.com/Kingsford-Group/scallop)

[IsoQuant](https://github.com/ablab/IsoQuant)

[Trinity](https://github.com/trinityrnaseq/trinityrnaseq)

[TransDecoder](https://github.com/TransDecoder/TransDecoder)

[Portcullis](https://github.com/EI-CoreBioinformatics/portcullis)

[Mikado](https://mikado.readthedocs.io/en/stable/)

## The workshop
The workshop is divided into four section

1. RNAseq mapping on the reference genome (Short-reads & Iso-Seq)
2. Homology and evidence-based prediction of protein coding genes (PCGs) using BRAKER2
3. *Ab Initio* annotation using Short-reads RNAseq & Iso-Seq
4. *De Novo* annotation using Short-reads RNAseq using Trinity
5. Metrics to evaluate annotations and Annotation Consensus using Mikado

## I hope this will be helpful, Have fun!
