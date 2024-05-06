# Merge all annotations into one.
There are multiple both empirical and intuitive evidence as to why integrating multiple annotation protocols enhances the outcome of genome annotation.

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
```bash
mikado prepare --json-conf configuration.yaml --procs 10
```

```bash
TransDecoder.LongOrfs -S -t mikado_prepared.fasta
```

```bash 
diamond blastp --ultra-sensitive --threads 30 --db ${SWISSPROTDB} --out mikado_prepared.fasta.transdecoder_dir/longest_orfs.Diamond.outfmt6.out --outfmt 6 \
	--evalue 1e-5 --max-target-seqs 1 --query mikado_prepared.fasta.transdecoder_dir/longest_orfs.pep
```

```bash
diamond blastx --ultra-sensitive --threads 30 --db ${SWISSPROTDB} --out mikado.diamond.xml --outfmt 5 --evalue 1e-5 --query mikado_prepared.fasta
```

```bash
gzip mikado.diamond.xml
```

```bash
TransDecoder.Predict -T 1000 -t mikado_prepared.fasta --retain_blastp_hits mikado_prepared.fasta.transdecoder_dir/longest_orfs.Diamond.outfmt6.out -v
```

```bash
mikado serialise -p 30 --json-conf configuration.yaml --xml mikado.blast.xml.gz --orfs mikado_prepared.fasta.transdecoder.bed
```

```bash
mikado pick --mode split --prefix $SPECIES -p 20 --json-conf configuration.yaml --subloci-out mikado.subloci.out.gff3 --output-dir $SPECIES.Mikado.split
```

```bash
cd $SPECIES.Mikado.split 
```

```bash
grep -P -v "Mikado_loci\tsuperlocus" mikado.loci.gff3 | grep -v '###' > tmp
```

```bash
gffread tmp -g $CHR -y mikado.loci.aa.fasta
```

```bash
busco -f -i mikado.loci.aa.fasta --cpu 30 -m prot -l $BUSCODIR/$BUSCODB --out run_mikado.loci.aa.$BUSCODB
cat run_mikado.loci.aa.$BUSCODB/short_summary.*.txt
```

```bash
mikado compare -r $ANNREF -p mikado.loci.gff3 -o MikadoVsRef 2> MikadoVsRef.errors.log
cat MikadoVsRef.stats
```

```bash
mikado util stats mikado.loci.gff3
```
