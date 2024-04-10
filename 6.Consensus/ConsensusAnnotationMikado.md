
```bash
mikado util convert -if bed12 -of gff3 $SPECIES.RNA.$CHRNAME.isoquant.out.bed > $SPECIES.RNA.$CHRNAME.isoquant.out.gff3
mikado util convert -if bed12 -of gff3 $SPECIES.minimap2.Trinity.bed > $SPECIES.minimap2.Trinity.gff3
```

```bash	
echo -e "$BRAKERGTF\tbr\tTrue\t3" > Mikado.conf
echo -e "$SPECIES.RNA.$CHRNAME.stringtie.out.gtf\tst\tTrue\t1" >> Mikado.conf
echo -e "$SPECIES.RNA.$CHRNAME.stringtie+BRAKER.out.gtf\tstbr\tTrue\t1" >> Mikado.conf
echo -e "$SPECIES.RNA.$CHRNAME.scallop.out.gtf\tscal\tTrue\t1" >> Mikado.conf
echo -e "$SPECIES.lonrRNA.$CHRNAME.stringtie.out.gtf\tstLR\tFalse\t1" > Mikado.conf 
echo -e "$SPECIES.RNA.$CHRNAME.isoquant.out.gff3\tisoq\tTrue\t3" >> Mikado.conf
echo -e "$SPECIES.minimap2.Trinity.gff3\tTrin\tFalse\t1" >> Mikado.conf 
```	

```bash
mikado configure --list Mikado.conf --scoring insects.yaml --reference $CHR --codon-table 0 \
	--junctions portcullis.0.5.pass.junctions.bed -bt ${SWISSPROTDB}.fasta configuration.yaml
```

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
