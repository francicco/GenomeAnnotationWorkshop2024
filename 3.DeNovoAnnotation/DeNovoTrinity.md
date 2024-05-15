# *De novo* transcriptome reconstruction using Trinity
`Trinity` is a well-established methodology for the reconstruction of the *de novo* transcriptomes from RNA-seq data. It implements many individual de Bruijn graphs, each representing the transcriptional complexity at a given gene or locus, to generate the full-length splicing isoforms from each gene and their paralogous, without the need for a reference genome assembly.

The reason why we are using `Trinity` when we have a genome assembly is because isoform reconstruction is not an easy task, and the complexity of primary mRNA splicing cannot be grasped with a single methodology. Let's for example imagine there are particularly long introns that are not easy to be fully identified with standard short-read mapping. Using `Trinity` for the de novo transcriptome reconstruction, with a subsequent mapping of those reconstructed isoforms, can bring a different type of information not fully exploited by short-reads.

## `Trinity` pipeline:
In this step of the annotation pipeline, we will therefore employ `Trinity` to do exactly that:
1. *De novo* transcriptome reconstruction
2. Mapping of those transcript onto the reference genome to integrate this type of information with other pieces of the puzzle.


In cases like this one where we have multiple libraries, we can run a small `for` loop in `bash` to create variables for the left and right pairs:
```bash
SRA1=''
SRA2=''

for SRA in $srSRAS; do
	SRA1="$SRA1 ${DATADIR}/${SRA}_1.extracted.fastq.gz"
	SRA2="$SRA2 ${DATADIR}/${SRA}_2.extracted.fastq.gz"
done
```

And subsequently run `Trinity`. This step will require some time, I provided you the final result using a standard run.
Note: if you have multiple libraries from multiple tissues and experiment you can combine all of them or run them separately. But in case you have very high sequencing depth you can tell `Trinity`, or use third part software (*e.g.*: [BBmap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/)), to [*in silico* normalise](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Insilico-Normalization) the data, which means removing reads from highly expressed genes, to make the computation more affordable in terms of RAM and CPU load.
```bash
Trinity --max_memory 150G --seqType fq --left $SRA1 --right $SRA2 --CPU $THREADS --output $SPECIES.trinity --full_cleanup --bflyCPU $THREADS
```


### Checking `Trinity` output:
Once you have your reconstructed transcripts you can map them with `minimap2`, making sure to convert the output in a `bed` file. You can use `|` to combine multiple commands:
```bash
minimap2 -ax splice:hq -C5 -t $THREADS $CHR $DATADIR/$SPECIES.Trinity.fasta.gz | samtools sort - | \
	samtools view -b -h -@ $THREADS - | bedtools bamtobed -bed12 -i - | awk '{ if ( $10 > 1 ) print $0 }' | awk '{ $4 = $4"."NR; print }' | sed 's/ /\t/g' > $SPECIES.minimap2.Trinity.bed
```

We can now use the `bed` file or whatever other format (*e.g.*: `gff`) you think is appropriate to extract the nucleotide sequence, not the aa. And check some stats.
```bash
bedtools getfasta -nameOnly -fo $SPECIES.RNA.$CHRNAME.trinity.out.fasta -fi $DATADIR/$SPECIES.$CHRNAME.fasta -bed $SPECIES.minimap2.Trinity.bed
```
Do you know why we're using nt instead of aa? 

Now we can use `BUSCO -m transcriptome` to check the completeness to check how complete is annotation might be. What would you expect?
```bash
busco -f -m transcriptome -i $SPECIES.RNA.$CHRNAME.trinity.out.fasta -f -o $SPECIES.RNA.$CHRNAME.trinity.Busco.$BUSCODB -l $BUSCODIR/$BUSCODB --cpu $THREADS
```

We also want to look at other metrics such as the representation of full-length reconstructed protein-coding genes. This can be done by searching the annotated transcripts against a database of known protein sequences.
you can find more of this explanation [here](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Counting-Full-Length-Trinity-Transcripts). In brief, what we want to do is to check how many transcripts appear to be full-length or nearly full-length compared with a reference proteome/transcriptome. Of course, doing this when closely related species are model organisms is a relatively straightforward procedure since one can easily align the transcripts to the reference transcripts and examine the length coverage. But when your target genome doesn't have close model organisms things can appear different without being necessarily a problem in the assembly itself. The idea of this metric is therefore to determine the number of unique top-matching proteins that align across more than X% of its length.

Since we have nt we can use `diamond bastx` to translate and map to proteins, providing a protein database indexed with `diamond makedb`.
```bash
diamond blastx  --ultra-sensitive --max-target-seqs 1 --threads $THREADS --query $SPECIES.RNA.$CHRNAME.trinity.out.fasta --outfmt 6 --db ${SWISSPROTDB} \
	--evalue 1e-5 --out $SPECIES.RNA.$CHRNAME.trinity.out.outfmt6
```

We can then compute the X% length against the best hits. We can also use the final tsv file to make a plot of the distribution on the frequency of X%.
```bash
analyze_blastPlus_topHit_coverage.pl $SPECIES.RNA.$CHRNAME.trinity.out.outfmt6 $SPECIES.RNA.$CHRNAME.trinity.out.fasta ${SWISSPROTDB}.fasta
grep -v '^#' $SPECIES.RNA.$CHRNAME.trinity.out.outfmt6.w_pct_hit_length | sed 's/^/Trinity\t/' > $SPECIES.RNA.$CHRNAME.trinity.out.outfmt6.w_pct_hit_length.tsv
```

You can also combine this table with the one from BRAKER
```bash
cat braker_utr.aa.out.outfmt6.w_pct_hit_length.tsv > AllAnnotations.aa.out.outfmt6.w_pct_hit_length.tsv
cat $SPECIES.RNA.$CHRNAME.trinity.out.outfmt6.w_pct_hit_length.tsv >> AllAnnotations.aa.out.outfmt6.w_pct_hit_length.tsv
```

And plot them together
```bash
Rscript Analyze_Diamond_topHit_coverage.R AllAnnotations.aa.out.outfmt6.w_pct_hit_length.tsv AllAnnotations.aa.out.outfmt6.w_pct_hit_length.png
```

`Mikado compare` can be used to check differences with the previous `BRAKER` predicion 
```bash
mikado compare -r $BRAKERGTF -p $SPECIES.minimap2.Trinity.bed -o TrinityVsBRAKER
cat TrinityVsBRAKER.stats
```

... and a reference.
```bash
mikado compare -r $ANNREF -p $SPECIES.minimap2.Trinity.bed -o TrinityVsRef
cat TrinityVsRef.stats
```

What can you tell? How is `Trinity` performing?
