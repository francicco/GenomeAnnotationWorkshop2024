```bash
stringtie $SPECIES.lonrRNA.$CHRNAME.sortedByCoord.out.bam -L \
		-o $SPECIES.lonrRNA.$CHRNAME.stringtie.out.gtf -p $THREADS -m 100 $STRAND #-G $BRAKERGTF
```

```bash
mikado util convert -if gtf -of bed12 $SPECIES.lonrRNA.$CHRNAME.stringtie.out.gtf | sed 's/;coding=False//' | sed 's/ID=//' > $SPECIES.lonrRNA.$CHRNAME.stringtie.out.bed
```

```bash
bedtools getfasta -name -fo $SPECIES.lonrRNA.$CHRNAME.stringtie.out.fasta -fi $DATADIR/$SPECIES.$CHRNAME.fasta -bed $SPECIES.lonrRNA.$CHRNAME.stringtie.out.bed
```

```bash
busco -f -i $SPECIES.lonrRNA.$CHRNAME.stringtie.out.fasta --cpu 30 -m transcriptome -l $BUSCODIR/$BUSCODB --out run_$SPECIES.lonrRNA.$CHRNAME.stringtie.out.$BUSCODB
```

```bash
 diamond blastx  --ultra-sensitive --max-target-seqs 1 --threads 20 --query $SPECIES.lonrRNA.$CHRNAME.stringtie.out.fasta --outfmt 6 --db ${SWISSPROTDB} \
		--evalue 1e-5 --out $SPECIES.lonrRNA.$CHRNAME.stringtie.out.outfmt6
```

```bash
$TRINITY_HOME/util/analyze_blastPlus_topHit_coverage.pl $SPECIES.lonrRNA.$CHRNAME.stringtie.out.outfmt6 $SPECIES.lonrRNA.$CHRNAME.stringtie.out.fasta ${SWISSPROTDB}.fasta
grep -v '^#' $SPECIES.lonrRNA.$CHRNAME.stringtie.out.outfmt6.w_pct_hit_length | sed 's/^/Scallop\t/' > $SPECIES.lonrRNA.$CHRNAME.stringtie.out.outfmt6.w_pct_hit_length.tsv
```


```bash
grep -vw 'CDS\|5UTR\|3UTR\|start_codon\|stop_codon' $BRAKERGTF > $SPECIES.BRAKER.isoquant.gff3
```

```bash
isoquant.py -d nanopore --report_novel_unspliced true --stranded none --threads $THREADS --genedb $SPECIES.BRAKER.isoquant.gff3 \
  --bam $SPECIES.lonrRNA.$CHRNAME.sortedByCoord.out.bam --reference $CHR --output $SPECIES.IsoQuant \
  --model_construction_strategy all --check_canonical --prefix ${SPECIES}IsoSeq --data_type $LONGTYPE \
  --check_canonical --splice_correction_strategy all --model_construction_strategy all --stranded none --count_exons
```

```bash
IsoQuantGTF2BED12v0.1.py -g ./$SPECIES.IsoQuant/${SPECIES}IsoSeq/${SPECIES}IsoSeq.transcript_models.gtf > $SPECIES.RNA.$CHRNAME.isoquant.out.bed
```

```bash
bedtools getfasta -name -fo $SPECIES.RNA.$CHRNAME.isoquant.out.fasta -fi $DATADIR/$SPECIES.$CHRNAME.fasta -bed $SPECIES.RNA.$CHRNAME.isoquant.out.bed
```

```bash
busco -m transcriptome -i $SPECIES.RNA.$CHRNAME.isoquant.out.fasta -f -o run_$SPECIES.lonrRNA.$CHRNAME.isoquant.out.$BUSCODB -l $BUSCODIR/$BUSCODB --cpu 20
```

```bash
diamond blastx  --ultra-sensitive --max-target-seqs 1 --threads 20 --query $SPECIES.RNA.$CHRNAME.isoquant.out.fasta --outfmt 6 --db ${SWISSPROTDB} \
		--evalue 1e-5 --out $SPECIES.RNA.$CHRNAME.isoquant.out.outfmt6
```

```bash
$TRINITY_HOME/util/analyze_blastPlus_topHit_coverage.pl $SPECIES.RNA.$CHRNAME.isoquant.out.outfmt6 $SPECIES.RNA.$CHRNAME.isoquant.out.fasta ${SWISSPROTDB}.fasta
grep -v '^#' $SPECIES.RNA.$CHRNAME.isoquant.out.outfmt6.w_pct_hit_length | sed 's/^/IsoQuant\t/' > $SPECIES.RNA.$CHRNAME.isoquant.out.outfmt6.w_pct_hit_length.tsv
```

```bash
head -n1 $DATADIR/$SPECIES.Chr1.proteins.*.outfmt6.w_pct_hit_length > $SPECIES.AbInitioLongReadsVsRef.w_pct_hit_length.tsv
cat $SPECIES.RNA.$CHRNAME.stringtie.out.outfmt6.w_pct_hit_length.tsv $SPECIES.RNA.$CHRNAME.isoquant.out.outfmt6.w_pct_hit_length.tsv >> $SPECIES.AbInitioLongReadsVsRef.w_pct_hit_length.tsv
Rscript ~/software/bioscript/Analyze_Diamond_topHit_coverage.R $SPECIES.AbInitioLongReadsVsRef.w_pct_hit_length.tsv $SPECIES.AbInitioLongReadsVsRef.w_pct_hit_length.png
```