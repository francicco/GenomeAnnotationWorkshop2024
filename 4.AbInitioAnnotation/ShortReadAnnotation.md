# *Ab Initio* transcriptome assembly Stringtie and Scallop
Both [`StringTie`](https://www.nature.com/articles/nbt.3122) and [`Scallop`](https://www.nature.com/articles/nbt.4020) are well-established methodologies for the reconstruction of the *ab initio* transcriptomes from short-reads RNA-seq data.

**StringTie** is based on the construction of a splice graph, in which nodes represent exons and edges (paths through the graph) represent possible splice variants.

<img width="968" alt="StringtieGraph" src="https://github.com/francicco/GenomeAnnotationWorkshop2024/assets/9006870/123439c2-c9a3-4375-9e0f-34772dd394ee">

The algorithm is similar to other transcript assemblers; however, instead of relying solely on a parsimony-based algorithm, which aims to generate the minimal number of transcripts, StringTie also incorporates transcript abundance to assess and filter out erroneous transcripts while simultaneously estimating their expression levels. Initially, StringTie clusters the reads, then constructs a splice graph for each cluster to identify transcripts. Subsequently, for each transcript, it constructs a separate flow network to estimate its expression level, employing a maximum flow algorithm.

**Scallop** also employs a graph-based algorithm, aiming to concurrently minimize discrepancies in read coverage and the number of expressed transcripts. This is achieved through an iterative process of decomposing vertices within the splice graph.

<img width="1019" alt="ScallopGraph" src="https://github.com/francicco/GenomeAnnotationWorkshop2024/assets/9006870/f5b02934-61f1-40e7-9eae-b4082f1f41ad">

Scallop capitalizes on reads spanning more than two exons by encoding them as phasing paths in the splice graph, which the algorithm endeavors to maintain during decomposition. When a vertex is fully covered by phasing paths, Scallop decomposes it by tackling two linear programming instances, optimizing to minimize coverage deviation. In cases where a vertex is not fully covered by phasing paths, Scallop formulates and resolves a subset-sum instance to either decrease transcript count or eliminate a false-positive edge.


## *Ab Initio* pipeline:
In this step of the annotation pipeline, we will therefore employ both methods, which will be merged in the final step of the annotation pipeline. The idea is to capitalize of the strengths of both metodologies.

### 1. We first run `StringTie` on the BAM file previously generate:
```bash
stringtie $SPECIES.RNA.$CHRNAME.Aligned.sortedByCoord.out.bam \
  -o $SPECIES.RNA.$CHRNAME.stringtie.out.gtf -p $THREADS -m 1000 $STRAND
```
*NOTE: because we have a gene prediction annotation we could also use it in the mapping step generating a **new** BAM file and use it which will likely have more reads mapped on spice-sites*

### Now we can convert the resulted `gtf` file produced by `StringTie` into a `bed12` file format...
```bash
mikado util convert -if gtf -of bed12 $SPECIES.RNA.$CHRNAME.stringtie.out.gtf | sed 's/;coding=False//' | sed 's/ID=//' > $SPECIES.RNA.$CHRNAME.stringtie.out.bed
```

### ...that will be used to extract the corresponding nucleotide sequences for each transcript.
```bash
bedtools getfasta -nameOnly -fo $SPECIES.RNA.$CHRNAME.stringtie.out.fasta -fi $DATADIR/$SPECIES.$CHRNAME.fasta -bed $SPECIES.RNA.$CHRNAME.stringtie.out.bed
```

### And we can use `busco` to try and infer how complete the annotation would be.
```bash
busco -f -m transcriptome -i $SPECIES.RNA.$CHRNAME.stringtie.out.fasta -f -o $SPECIES.RNA.$CHRNAME.stringtie.Busco.$BUSCODB -l $BUSCODIR/$BUSCODB --cpu 20
```

### We can also `blastx` the transcripts...
```bash
diamond blastx  --ultra-sensitive --max-target-seqs 1 --threads 20 --query $SPECIES.RNA.$CHRNAME.stringtie.out.fasta --outfmt 6 --db ${SWISSPROTDB} \
	--evalue 1e-5 --out $SPECIES.RNA.$CHRNAME.stringtie.out.outfmt6
```

### ...to check the completeness of the transcripts, just like we did during the *de novo* step.
```bash
$TRINITY_HOME/util/analyze_blastPlus_topHit_coverage.pl $SPECIES.RNA.$CHRNAME.stringtie.out.outfmt6 $SPECIES.RNA.$CHRNAME.stringtie.out.fasta ${SWISSPROTDB}.fasta
grep -v '^#' $SPECIES.RNA.$CHRNAME.stringtie.out.outfmt6.w_pct_hit_length | sed 's/^/Stringtie\t/' > $SPECIES.RNA.$CHRNAME.stringtie.out.outfmt6.w_pct_hit_length.tsv
```

### 2. Another strategy would be add an annotation, *i.e.*: the `BRAKER` one, to guide the transcript reconstruction. 
```bash
stringtie $SPECIES.RNA.$CHRNAME.Aligned.sortedByCoord.out.bam \
    -o $SPECIES.RNA.$CHRNAME.stringtie+BRAKER.out.gtf -p $THREADS -m 1000 $STRAND -G $BRAKERGTF
```

### We can now follow the same steps to check the quality of this strategy.
```bash
mikado util convert -if gtf -of bed12 $SPECIES.RNA.$CHRNAME.stringtie+BRAKER.out.gtf | sed 's/;coding=False//' | sed 's/ID=//' > $SPECIES.RNA.$CHRNAME.stringtie+BRAKER.out.bed
```

```bash
bedtools getfasta -nameOnly -fo $SPECIES.RNA.$CHRNAME.stringtie+BRAKER.out.fasta -fi $DATADIR/$SPECIES.$CHRNAME.fasta -bed $SPECIES.RNA.$CHRNAME.stringtie+BRAKER.out.bed
```

```bash
busco -f -m transcriptome -i $SPECIES.RNA.$CHRNAME.stringtie+BRAKER.out.fasta -f -o $SPECIES.RNA.$CHRNAME.stringtie+BRAKER.Busco.$BUSCODB -l $BUSCODIR/$BUSCODB --cpu 20
```

```bash
diamond blastx  --ultra-sensitive --max-target-seqs 1 --threads 20 --query $SPECIES.RNA.$CHRNAME.stringtie+BRAKER.out.fasta --outfmt 6 --db ${SWISSPROTDB} \
    --evalue 1e-5 --out $SPECIES.RNA.$CHRNAME.stringtie+BRAKER.out.outfmt6
```

```bash
$TRINITY_HOME/util/analyze_blastPlus_topHit_coverage.pl $SPECIES.RNA.$CHRNAME.stringtie+BRAKER.out.outfmt6 $SPECIES.RNA.$CHRNAME.stringtie+BRAKER.out.fasta ${SWISSPROTDB}.fasta
grep -v '^#' $SPECIES.RNA.$CHRNAME.stringtie+BRAKER.out.outfmt6.w_pct_hit_length | sed 's/^/Stringtie+BRAKER\t/' > $SPECIES.RNA.$CHRNAME.stringtie+BRAKER.out.outfmt6.w_pct_hit_length.tsv
```

### 3. We can repeat the same steps with `Scallp`.
```bash
scallop -i $SPECIES.RNA.$CHRNAME.Aligned.sortedByCoord.out.bam -o $SPECIES.RNA.$CHRNAME.scallop.out.gtf --library_type $LIBTYPE > $SPECIES.RNA.$CHRNAME.scallop.log
```

```bash
mikado util convert -if gtf -of bed12 $SPECIES.RNA.$CHRNAME.scallop.out.gtf | sed 's/;coding=False//' | sed 's/ID=//' > $SPECIES.RNA.$CHRNAME.scallop.out.bed
```

```bash
bedtools getfasta -nameOnly -fo $SPECIES.RNA.$CHRNAME.scallop.out.fasta -fi $DATADIR/$SPECIES.$CHRNAME.fasta -bed $SPECIES.RNA.$CHRNAME.scallop.out.bed
```

```bash
busco -f -m transcriptome -i $SPECIES.RNA.$CHRNAME.scallop.out.fasta -f -o $SPECIES.RNA.$CHRNAME.scallop.Busco.$BUSCODB -l $BUSCODIR/$BUSCODB --cpu 20
```

```bash
diamond blastx  --ultra-sensitive --max-target-seqs 1 --threads 20 --query $SPECIES.RNA.$CHRNAME.scallop.out.fasta --outfmt 6 --db ${SWISSPROTDB} \
	--evalue 1e-5 --out $SPECIES.RNA.$CHRNAME.scallop.out.outfmt6
```

```bash
$TRINITY_HOME/util/analyze_blastPlus_topHit_coverage.pl $SPECIES.RNA.$CHRNAME.scallop.out.outfmt6 $SPECIES.RNA.$CHRNAME.scallop.out.fasta ${SWISSPROTDB}.fasta
grep -v '^#' $SPECIES.RNA.$CHRNAME.scallop.out.outfmt6.w_pct_hit_length | sed 's/^/Scallop\t/' > $SPECIES.RNA.$CHRNAME.scallop.out.outfmt6.w_pct_hit_length.tsv
```

### Finally we can concatenate all these stats obtained for each strategy and look at the results.
```bash
head -n1 $DATADIR/$SPECIES.$CHRNAME.proteins.*.outfmt6.w_pct_hit_length > $SPECIES.AbInitioVsRef.w_pct_hit_length.tsv
cat $SPECIES.RNA.$CHRNAME.stringtie.out.outfmt6.w_pct_hit_length.tsv $SPECIES.RNA.$CHRNAME.stringtie+BRAKER.out.outfmt6.w_pct_hit_length.tsv $SPECIES.RNA.$CHRNAME.scallop.out.outfmt6.w_pct_hit_length.tsv >> $SPECIES.AbInitioVsRef.w_pct_hit_length.tsv
```

### We provide you with a `R` script: `Analyze_Diamond_topHit_coverage.R`, located in the `Scripts` folder, to plot the data.
```bash
Rscript Analyze_Diamond_topHit_coverage.R $SPECIES.AbInitioVsRef.w_pct_hit_length.tsv $SPECIES.AbInitioVsRef.w_pct_hit_length.png
```


```bash
mikado compare -r $ANNREF -p $SPECIES.RNA.$CHRNAME.stringtie.out.gtf -o StringtieVsRef
cat StringtieVsRef.stats
```

```bash
mikado compare -r $ANNREF -p $SPECIES.RNA.$CHRNAME.stringtie+BRAKER.out.gtf -o StringtieBrakerVsRef
cat StringtieBrakerVsRef.stats
```

```bash
mikado compare -r $ANNREF -p $SPECIES.RNA.$CHRNAME.scallop.out.gtf -o ScallopVsRef
cat ScallopVsRef.stats
```

```bash
mikado compare -r $BRAKERGTF -p $SPECIES.RNA.$CHRNAME.stringtie.out.gtf -o StringtieVsBRAKER
cat StringtieVsBRAKER.stats
```

```bash
mikado compare -r $BRAKERGTF -p $SPECIES.RNA.$CHRNAME.stringtie+BRAKER.out.gtf -o StringtieBrakerVsBRAKER
cat StringtieBrakerVsBRAKER.stats
```

```bash
mikado compare -r $BRAKERGTF -p $SPECIES.RNA.$CHRNAME.scallop.out.gtf -o ScallopVsBRAKER
cat ScallopVsBRAKER.stats
```
