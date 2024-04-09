```bash
SRA1=''
SRA2=''

for SRA in $srSRAS; do
	SRA1="$SRA1 ${DATADIR}/${SRA}_1.extracted.fastq.gz"
	SRA2="$SRA2 ${DATADIR}/${SRA}_2.extracted.fastq.gz"
done
```

```bash
Trinity --max_memory 150G --seqType fq --left $SRA1 --right $SRA2 --CPU $THREADS --output $SPECIES.trinity --full_cleanup --bflyCPU $THREADS"
```

```bash
minimap2 -ax splice:hq -C5 -t $THREADS $CHR $DATADIR/$SPECIES.Trinity.fasta.gz | samtools sort - | \
	samtools view -b -h -@ $THREADS - | bedtools bamtobed -bed12 -i - | awk '{ if ( $10 > 1 ) print $0 }' | awk '{ $4 = $4"."NR; print }' | sed 's/ /\t/g' > $SPECIES.minimap2.Trinity.bed
```

```bash
bedtools getfasta -nameOnly -fo $SPECIES.RNA.$CHRNAME.trinity.out.fasta -fi $DATADIR/$SPECIES.$CHRNAME.fasta -bed $SPECIES.minimap2.Trinity.bed
```

```bash
busco -f -m transcriptome -i $SPECIES.RNA.$CHRNAME.trinity.out.fasta -f -o $SPECIES.RNA.$CHRNAME.trinity.Busco.$BUSCODB -l $BUSCODIR/$BUSCODB --cpu 20
```

```bash
diamond blastx  --ultra-sensitive --max-target-seqs 1 --threads 20 --query $SPECIES.RNA.$CHRNAME.trinity.out.fasta --outfmt 6 --db ${SWISSPROTDB} \
	--evalue 1e-5 --out $SPECIES.RNA.$CHRNAME.trinity.out.outfmt6
```

```bash
$TRINITY_HOME/util/analyze_blastPlus_topHit_coverage.pl $SPECIES.RNA.$CHRNAME.trinity.out.outfmt6 $SPECIES.RNA.$CHRNAME.trinity.out.fasta ${SWISSPROTDB}.fasta
grep -v '^#' $SPECIES.RNA.$CHRNAME.trinity.out.outfmt6.w_pct_hit_length | sed 's/^/Trinity\t/' > $SPECIES.RNA.$CHRNAME.trinity.out.outfmt6.w_pct_hit_length.tsv
```

```bash
mikado compare -r $BRAKERGTF -p $SPECIES.minimap2.Trinity.bed -o TrinityVsBRAKER
cat TrinityVsBRAKER.stats
```

```bash
mikado compare -r $ANNREF -p $SPECIES.minimap2.Trinity.bed -o TrinityVsRef
cat TrinityVsRef.stats
```
