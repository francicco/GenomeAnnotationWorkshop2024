# Prediction of Protein Coding Genes using the BRAKER pipeline
`BRAKER` is a fully automated pipeline that generates a training set of gene which is then used to with `GeneMark-ES/ET/EP/ETP` and `AUGUSTUS`. It uses both evidence-based information from both RNA-Seq experiments and/or protein homology.

But it can also used to obtain good gene prediction accuracy in the absence of proteins from closely related species and RNA-Seq data.

In this step, you will leverage RNAseq data mapped in the previous step and protein mapping here to execute (I should say simulate) a BRAKER analysis. The whole pipeline can be quite slow, depending on the amount of data and the size of the genome you want to annotate. Here you will try to execute a run but you won't be able to finish it.
But, don't worry I'll provide a pre-run analysis that you can explore and use for the next steps of the annotation pipeline.

## Beginning of the pipeline:
This section will deal with mapping proteins from closely related species of your choice and format the output to be digested from `BRAKER`.
In this example, you will use proteins from `Swiss prot DB`, and you will map that with `miniprot`.
The mapper is very fast and can provide multiple types of output. Execute `miniprot -h` to check.

For example the parameter `-G` controls the expected size of the longest introns, you need to adjust accordingly to your genome. One way to check the results is to load the output on `IGV` and check in the different maps have too much overlap or not enough coding exons.
```bash
miniprot --aln --gff --trans -t $THREADS --trans -G 50000 $CHR $SWISSPROTDB.fasta > $SPECIES.MiniProt.gff
```

Once you are satisfied with your mapping you can start to convert the mapping into a `gff` format file and use it as input in `BRAKER`.
The steps involve a fist conversion into a nucleotide fasta file, which is then remapped with `minimap2` to obrain a `SAM` file which is than converted into a `psl` and finallu into a `gff`

To extract the nucleotides from a `miniprot` run you can use any tool/script you want. In this case I wrote a small python script `Miniprot2SplicedNucl.py` to do the job.
```bash
Miniprot2SplicedNucl.py -g $SPECIES.MiniProt.gff > $SPECIES.MiniProt.nt.fasta
```

Now run `minimap2`:
```bash
minimap2 -t $THREADS -ax splice:hq $CHR $SPECIES.MiniProt.nt.fasta > $SPECIES.MiniProt.nt.sam
```

Convert the `sam` into a `psl`:
```bash
sam2psl.py -i $SPECIES.MiniProt.nt.sam -s -o $SPECIES.MiniProt.psl
```

Sort it:
```bash
cat $SPECIES.MiniProt.psl | sort -n -k 16,16 | sort -s -k 14,14 > $SPECIES.MiniProt.sorted.psl
```

and finally convert it check it:
```bash
blat2hints.pl --nomult --in=$SPECIES.MiniProt.sorted.psl --out=$SPECIES.MiniProt.hints.gff
wc -l $SPECIES.MiniProt.hints.gff
```


Ok, now that all the `BRAKER` inputs are prepared we are ready to execute `BRAKER`. One important step here is to check whether all environmental path are correctly set. For example `$AUGUSTUS_SCRIPTS_PATH` or `$PYTHON3_PATH`.
You can also manually link these paths to `BRAKER` using specific flags like `--PYTHON3_PATH`.
`BRAKER` has other useful settings now, like `--UTR` option, which allow you to predict also the UTRs, in case you're interested in them. You could also check if the analyses run with or without it produce the same results.

Now run `BRAKER`:
```bash
braker.pl --UTR=on --threads $THREADS \
  --workingdir=. --alternatives-from-evidence=true --crf \
  --nocleanup --species=${SPECIES}_$CHRNAME --softmasking \
  --genome=$CHR --gff3 --hints=$SPECIES.MiniProt.hints.gff \
  --AUGUSTUS_ab_initio -grass --verbosity=4 \
  --bam=$SPECIES.RNA.$CHRNAME.Aligned.sortedByCoord.out.bam 
```

### Now the run has finished you can check the quality of it, for example asking if enough core genes (`BUSCO`) have been predicted and if there's room for improvement.
Running `BUSCO` on the genome AND on the annotation can give you an idea how efficient the annotation method you're using is performing. In theory you want to approach the values found with `BUSCO -m genome`

Now extract amino-acid sequences from the `BRAKER` annotation:
```bash
gffread $BRAKERGTF -g $CHR -y braker_utr.aa.fasta
```

... and run `BUSCO -m protein`:
```bash
busco -f -i braker_utr.aa.fasta --cpu $THREADS -m prot -l $BUSCODIR/$BUSCODB --out run_braker_utr.aa.$BUSCODB 
```

### Let's have a look at the fragmentation of the genes! One way to do it would be to blast our set of proteins to a reference DB and check the distribution of alignment. For this we can you uniprot. If you haven't done it already you can download the fasta file:
```bash
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
```

Move `uniprot_sprot.fasta.gz` to your `$DATADIR` and unpack the `gz` file:
```bash
cp uniprot_sprot.fasta.gz $DATADIR
cd $DATADIR
gunzip uniprot_sprot.fasta.gz
```

Now we have to build the index for `diamond` using the subcommand `makedb`.
```bash
diamond makedb --in uniprot_sprot.fasta --db uniprot_sprot
```

Now we are ready to `blastp` our proteome.
```bash
diamond blastp  --ultra-sensitive --max-target-seqs 1 --threads $THREADS --query braker_utr.aa.fasta --outfmt 6 --db ${SWISSPROTDB} \
	--evalue 1e-5 --out braker_utr.aa.out.outfmt6
```

We can then compute the X% length against the best hits. We can also use the final tsv file to make a plot of the distribution on the frequency of X%.
```bash
analyze_blastPlus_topHit_coverage.pl braker_utr.aa.out.outfmt6 braker_utr.aa.asta ${SWISSPROTDB}.fasta
sed 's/^/Blaker\t/' > braker_utr.aa.out.outfmt6.w_pct_hit_length.tsv
```

In case you already have a reference annotation you can compare it with the `BRAKER` output, to see how the analysis performed.
In normal circumstances, you won't have an annotation. So in this case we're using this step to show you the limitations of each step of the pipeline.

The comparison can be done with a sub-program of `mikado`, `mikado compare` to obtain all the stats. You can also use the comparison to compare the different steps of your pipeline.
```bash
mikado compare -r $ANNREF -p $BRAKERGTF -o BrakerVsRef 2> BrakerVsRef.errors.log
cat BrakerVsRef.stats
```

`Mikado util` instead extracts multiple pieces of information from your annotation. Execute it a have a look.
```bash
mikado util stats $BRAKERGTF
```

Well done, the `BRAKER` step is complete!!!

