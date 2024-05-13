# Merge all annotations into one.
There are multiple both empirical and intuitive evidence as to why integrating multiple annotation protocols enhances the outcome of genome annotation.

<img width="1282" alt="AnnPipeLineCicconardiEtAl" src="https://github.com/francicco/GenomeAnnotationWorkshop2024/assets/9006870/3c1b0332-7644-4a95-b367-9ecc1a96fb6a">

1. **Redundancy Reduction:** Each annotation protocol may have its own strengths and weaknesses, leading to both false positives and false negatives. Integrating multiple protocols allows researchers to cross-validate annotations and reduce the likelihood of errors. Consistent annotations across multiple protocols increase confidence in the accuracy of the final annotation.

2. **Resolution of Ambiguities:** Some regions of the genome may be challenging to annotate due to complex genomic structures, repetitive sequences, or ambiguities in the data. Integrating multiple annotation protocols can help resolve these ambiguities by providing complementary information from different experimental approaches or computational methods.

4. **Improved Gene Prediction:** Gene prediction algorithms may differ in their sensitivity and specificity, leading to variations in predicted gene structures. Integrating multiple annotation protocols can improve gene prediction by incorporating evidence from different sources, such as transcriptomic data, homology-based approaches, and ab initio gene prediction algorithms.

5. **Discovery of Novel Features:** Integrating multiple annotation protocols increases the likelihood of discovering novel genomic features that may have been overlooked by individual protocols. For example, integrating RNA-seq data with computational predictions and comparative genomics can lead to the identification of novel genes, alternative splicing events, and non-coding RNAs.

Overall, the integration of multiple annotation protocols maximizes the depth, accuracy, and completeness of genome annotation, providing researchers with a more comprehensive understanding of the genomic landscape.

Several algorithms can achieve this. Among the most popular there are [`MAKER`](https://weatherby.genetics.utah.edu/MAKER/wiki/index.php/MAKER_Tutorial_for_WGS_Assembly_and_Annotation_Winter_School_2018) and `EVidenceModeler` ([`EVM`](https://github.com/EVidenceModeler/EVidenceModeler/wiki).

In this workshop we decided to shows [`Mikado`](https://mikado.readthedocs.io/en/stable/). This tool can defines loci, score transcripts, determine a representative transcript for each locus, and finally return a set of gene models filtered to individual requirements, in the attempt to remove transcripts that are chimeric, fragmented, or with short or disrupted coding sequences. You can read more on the published [article](https://academic.oup.com/gigascience/article/7/8/giy093/5057872).

Briefly, `Mikado` is a three-part protocol: `configure`, `prepare` and `pick`. 

## 1. `Configure` is the first step in it you gather all info together including a [scoring file](https://mikado.readthedocs.io/en/stable/Tutorial/Scoring_tutorial/).
But Let's first create a configuration file (`Mikado.conf`). It is the main file, you can have a better look at the official wiki on [how to create it](https://mikado.readthedocs.io/en/stable/Tutorial/#creating-the-configuration-file-for-mikado).

In few words it is a 3-field tab-delimited (plus 4 optionals) file describing our annotations.

1. The file location and name (if no folder is specified, Mikado will look for each file in the current working directory)
2. An alias associated with the file, which has to be unique
3. A binary flag (True / False) indicating whether the annotation is strand-specific or not

Following the optional fields correspond to:

4. A score associated with that sample. All transcripts associated with the label will have their score corrected by the value on this field. So eg. in this example all Stringtie models will receive an additional point (`1`), and all Trinity models will be penalised by one point(-1). `Isoquant` has instead a bonus of `3`.
5. A binary flag (`True` / `False`) defining whether the sample is a reference or not.
6. A binary flag (`True` / `False`) defining whether to exclude redundant models or not.
7. A binary flag (`True` / `False`) indicating whether `Mikado` `prepare` should strip the `CDS` of faulty models, but otherwise keep their cDNA structure in the final output (`True`) or whether instead it should completely discard such models (`False`).
8. A binary flag (`True` / `False`) instructing `Mikado` about whether the chimera split routine should be skipped for these models (`True`) or if instead it should proceed normally (`False`).

Ok! Now let's proceed by converting into `gff3` the `bed12` files we otained previously
```bash
mikado util convert -if bed12 -of gff3 $SPECIES.RNA.$CHRNAME.isoquant.out.bed > $SPECIES.RNA.$CHRNAME.isoquant.out.gff3
mikado util convert -if bed12 -of gff3 $SPECIES.minimap2.Trinity.bed > $SPECIES.minimap2.Trinity.gff3
```

And let's credte `Mikado.conf`
```bash	
echo -e "$BRAKERGTF\tbr\tTrue\t3" > Mikado.conf
echo -e "$SPECIES.RNA.$CHRNAME.stringtie.out.gtf\tst\tFalse\t1" >> Mikado.conf
echo -e "$SPECIES.RNA.$CHRNAME.stringtie+BRAKER.out.gtf\tstbr\ttFalse\t1" >> Mikado.conf
echo -e "$SPECIES.RNA.$CHRNAME.scallop.out.gtf\tscal\ttFalse\t1" >> Mikado.conf
echo -e "$SPECIES.lonrRNA.$CHRNAME.stringtie.out.gtf\tstLR\tFalse\t1" > Mikado.conf 
echo -e "$SPECIES.RNA.$CHRNAME.isoquant.out.gff3\tisoq\tTrue\t3" >> Mikado.conf
echo -e "$SPECIES.minimap2.Trinity.gff3\tTrin\tFalse\t-1" >> Mikado.conf 
```	

Now we can run `mikado configure`, to note here is the addition of the annotated splice-sites `portcullis.0.5.pass.junctions.bed`, and the proteome `${SWISSPROTDB}.fasta` that will be used by `Mikado` to check for chimeras and ORFs.
Another important feature in `Mikado` is the scoring file a [`yaml` formatted file](https://spacelift.io/blog/yaml), which is key to `Mikado`. The scoring file is a [user-defined configuration file](https://mikado.readthedocs.io/en/stable/Scoring_files/) that is used by `Mikado` to filter the annotations and to reject or accept transcript.
There's no ideal scoring system, and depending on the need or the genomic property can be customized. In the `Mikado` distribution there are some available such as the one used in this example (`insects.yaml`). In the repository of this workshop we also provide a scoring file (`Insects-Mod.yaml`). It is much more complex than a default one. I used it for the annotation of few ant species. You can have a look and experiment with it.
**There's no right or wrong, there is trial and error!**

```bash
mikado configure --list Mikado.conf --scoring insects.yaml --reference $CHR --codon-table 0 \
	--junctions portcullis.0.5.pass.junctions.bed -bt ${SWISSPROTDB}.fasta configuration.yaml
```

The output is a configuration file: `configuration.yaml`, which will contain all the info for `Mikado`.

## 2. The subsequent step involves running `prepare` to create a sorted, non-redundant `GTF` with all the input annotations.
```bash
mikado prepare --json-conf configuration.yaml --procs $THREADS
```

This step produces a `gtf` (`mikado_prepared.gtf`) and a `fasta` file (`mikado_prepared.fasta`), and we can use the fasta to annotate the ORFs within each transcript.
For that we can use `TransDecoder`.
```bash
TransDecoder.LongOrfs -t mikado_prepared.fasta
```

On the predicted ORFs we can run `diamond blastp` to check for homology. This is another source of info that `TransDecoder` will use to assign the best ORF per transcript.
```bash 
diamond blastp --ultra-sensitive --threads $THREADS --db ${SWISSPROTDB} --out mikado_prepared.fasta.transdecoder_dir/longest_orfs.Diamond.outfmt6.out --outfmt 6 \
	--evalue 1e-5 --max-target-seqs 1 --query mikado_prepared.fasta.transdecoder_dir/longest_orfs.pep
```
*Note: `TransDecoder` can also make use of PFAM, using `hmmscan`, but because this is a relatively time-consuming step, which for time constraints we do not run here, we advice its use*

Now let's give to `TransDecoder` all the info we have!
```bash
TransDecoder.Predict -T 1000 -t mikado_prepared.fasta --retain_blastp_hits mikado_prepared.fasta.transdecoder_dir/longest_orfs.Diamond.outfmt6.out -v
```

`Mikado` will also accept a full `blastx` search in a different format in `xml` format (`--outfmt 5`). `Diamond` is relatively quick, so... let's give it a go!
```bash
diamond blastx --ultra-sensitive --threads $THREADS --db ${SWISSPROTDB} --out mikado.diamond.xml --outfmt 5 --evalue 1e-5 --query mikado_prepared.fasta
gzip mikado.diamond.xml
```

## 3. Finally, `pick` will integrate the data with the positional and structural data present in the `GTF` file to select the best transcript models.
But before `pick`, let's covert all into a `SQLite` database...
```bash
mikado serialise -p $THREADS --json-conf configuration.yaml --xml mikado.blast.xml.gz --orfs mikado_prepared.fasta.transdecoder.bed
```

...and run `pick`!
Wait! `pick` has five different modes: `nosplit`, `stringent`, `lenient`, `permissive`, and `split`.
- `nosplit`: keep the transcripts whole.
- `stringent`: split multi-orf transcripts if two consecutive ORFs have both BLAST hits and none of those hits is against the same target.
- `lenient`: split multi-orf transcripts as in stringent, and additionally, also when either of the ORFs lacks a BLAST hit (but not both).
- `permissive`: like lenient, but also split when both ORFs lack BLAST hits.
- `split`: split multi-orf transcripts regardless of what BLAST data is available.

Now *pick* the one you like... or check for differences!
```bash
MODE=split
mikado pick --mode $MODE --prefix $SPECIES -p $THREADS --json-conf configuration.yaml --subloci-out mikado.subloci.out.gff3 --output-dir $SPECIES.Mikado.split
```

## Now that you concluded the pipeline you can have a look at the results
```bash
cd $SPECIES.Mikado.$MODE 
```

Do e bit of cleaning of the resulted `gff3`...
```bash
grep -P -v "Mikado_loci\tsuperlocus" mikado.loci.gff3 | grep -v '###' > tmp
```

Extract the amino-acid sequences:
```bash
gffread tmp -g $CHR -y mikado.loci.aa.fasta
```

And run `BUSCO`...
```bash
busco -f -i mikado.loci.aa.fasta --cpu $THREADS -m prot -l $BUSCODIR/$BUSCODB --out run_mikado.loci.aa.$BUSCODB
cat run_mikado.loci.aa.$BUSCODB/short_summary.*.txt
```

and compare it with our reference annotation, if available.
```bash
mikado compare -r $ANNREF -p mikado.loci.gff3 -o MikadoVsRef 2> MikadoVsRef.errors.log
cat MikadoVsRef.stats
```

`Mikado` has also various utilities such as `stats`. You can use it to check stats in your annotations:
```bash
mikado util stats mikado.loci.gff3
```

## How does your annotation look? Can you improve it?
