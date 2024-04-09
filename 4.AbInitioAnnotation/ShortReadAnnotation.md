```bash
stringtie $SPECIES.RNA.$CHRNAME.Aligned.sortedByCoord.out.bam \
  -o $SPECIES.RNA.$CHRNAME.stringtie.out.gtf -p $THREADS -m 1000 $STRAND
```

```bash
mikado util convert -if gtf -of bed12 $SPECIES.RNA.$CHRNAME.stringtie.out.gtf | sed 's/;coding=False//' | sed 's/ID=//' > $SPECIES.RNA.$CHRNAME.stringtie.out.bed
```

```bash
bedtools getfasta -nameOnly -fo $SPECIES.RNA.$CHRNAME.stringtie.out.fasta -fi $DATADIR/$SPECIES.$CHRNAME.fasta -bed $SPECIES.RNA.$CHRNAME.stringtie.out.bed
```

```bash
busco -f -m transcriptome -i $SPECIES.RNA.$CHRNAME.stringtie.out.fasta -f -o $SPECIES.RNA.$CHRNAME.stringtie.Busco.$BUSCODB -l $BUSCODIR/$BUSCODB --cpu 20
```

```bash
diamond blastx  --ultra-sensitive --max-target-seqs 1 --threads 20 --query $SPECIES.RNA.$CHRNAME.stringtie.out.fasta --outfmt 6 --db ${SWISSPROTDB} \
	--evalue 1e-5 --out $SPECIES.RNA.$CHRNAME.stringtie.out.outfmt6
```

```bash
$TRINITY_HOME/util/analyze_blastPlus_topHit_coverage.pl $SPECIES.RNA.$CHRNAME.stringtie.out.outfmt6 $SPECIES.RNA.$CHRNAME.stringtie.out.fasta ${SWISSPROTDB}.fasta
grep -v '^#' $SPECIES.RNA.$CHRNAME.stringtie.out.outfmt6.w_pct_hit_length | sed 's/^/Stringtie\t/' > $SPECIES.RNA.$CHRNAME.stringtie.out.outfmt6.w_pct_hit_length.tsv
```

```bash
stringtie $SPECIES.RNA.$CHRNAME.Aligned.sortedByCoord.out.bam \
    -o $SPECIES.RNA.$CHRNAME.stringtie+BRAKER.out.gtf -p $THREADS -m 1000 $STRAND -G $BRAKERGTF
```

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

```bash
head -n1 $DATADIR/$SPECIES.$CHRNAME.proteins.*.outfmt6.w_pct_hit_length > $SPECIES.AbInitioVsRef.w_pct_hit_length.tsv
cat $SPECIES.RNA.$CHRNAME.stringtie.out.outfmt6.w_pct_hit_length.tsv $SPECIES.RNA.$CHRNAME.stringtie+BRAKER.out.outfmt6.w_pct_hit_length.tsv $SPECIES.RNA.$CHRNAME.scallop.out.outfmt6.w_pct_hit_length.tsv >> $SPECIES.AbInitioVsRef.w_pct_hit_length.tsv
```

```bash
Rscript ~/software/bioscript/Analyze_Diamond_topHit_coverage.R $SPECIES.AbInitioVsRef.w_pct_hit_length.tsv $SPECIES.AbInitioVsRef.w_pct_hit_length.png
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
