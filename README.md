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

[STAR](https://github.com/alexdobin/STAR) | [Minimap2](https://github.com/lh3/minimap2) | [Samtools](https://github.com/samtools/samtools) | [Bedtools](https://github.com/arq5x/bedtools2) | [Diamond](https://github.com/bbuchfink/diamond) | [Miniprot](https://github.com/lh3/miniprot) | [BRAKER](https://github.com/Gaius-Augustus/BRAKER) | [Cufflinks](https://cole-trapnell-lab.github.io/cufflinks/install/) | [BUSCO](https://busco.ezlab.org/) | [compleasm](https://github.com/huangnengCSU/compleasm) | [Stringtie](https://ccb.jhu.edu/software/stringtie/) | [Scallop](https://github.com/Kingsford-Group/scallop) | [IsoQuant](https://github.com/ablab/IsoQuant) | [Trinity](https://github.com/trinityrnaseq/trinityrnaseq) | [TransDecoder](https://github.com/TransDecoder/TransDecoder) | [Portcullis](https://github.com/EI-CoreBioinformatics/portcullis) | [Mikado](https://mikado.readthedocs.io/en/stable/) | [`Miniprot2SplicedNucl.py`](https://github.com/francicco/GenomeAnnotationWorkshop2024/blob/main/Scripts/Miniprot2SplicedNucl.py) | [`IsoQuantGTF2BED12v0.1.py`](https://github.com/francicco/GenomeAnnotationWorkshop2024/blob/main/Scripts/IsoQuantGTF2BED12v0.1.py) | [`Analyze_Diamond_topHit_coverage.R`](https://github.com/francicco/GenomeAnnotationWorkshop2024/blob/main/Scripts/Miniprot2SplicedNucl.py) | [UniProtDB](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz) ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

## The workshop
The workshop is divided into four section

1. [RNAseq mapping on the reference genome (Short-reads & Iso-Seq)](https://github.com/francicco/GenomeAnnotationWorkshop2024/blob/main/1.Mapping/1.MappingStep.md)
2. [Homology and evidence-based prediction of protein coding genes (PCGs) using BRAKER2](https://github.com/francicco/GenomeAnnotationWorkshop2024/blob/main/2.Prediction/BRAKER.md)
3. [*De Novo* annotation using Short-reads RNAseq using Trinity](https://github.com/francicco/GenomeAnnotationWorkshop2024/blob/main/3.DeNovoAnnotation/DeNovoTrinity.md)
4. *Ab Initio* annotation using [Short-reads](https://github.com/francicco/GenomeAnnotationWorkshop2024/blob/main/4.AbInitioAnnotation/1.ShortReadAnnotation.md) RNAseq & [Long-reads RNA-Seq](https://github.com/francicco/GenomeAnnotationWorkshop2024/blob/main/4.AbInitioAnnotation/2.LongReadAnnotation.md)
5. Metrics to evaluate annotations, [Splice-site filtering](https://github.com/francicco/GenomeAnnotationWorkshop2024/blob/main/5.SpliceJunctionFiltering/5.PortcullisRun.md), and [Annotation Consensus using Mikado](https://github.com/francicco/GenomeAnnotationWorkshop2024/blob/main/6.Consensus/ConsensusAnnotationMikado.md)

All data is available at this [link](https://drive.google.com/drive/folders/1IreMRHaOa1kvOomyjoEm8xFw1fmOR-oK?usp=drive_link), but don't forget to [set up your environment](https://github.com/francicco/GenomeAnnotationWorkshop2024/blob/main/0.VariableSetting.md)!!!

## I hope this will be useful, Have fun!
