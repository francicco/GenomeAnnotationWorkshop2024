#!/usr/bin/env python2

import optparse
import re
from collections import defaultdict

################################# Command line options

desc='Converte a IsoQuant GTF (transcript/exon) into a BED12 file'

parser = optparse.OptionParser(description=desc, version='%prog version 0.1 - 25-09-2016 - Author: FCicconardi')

parser.add_option('-g', '--GTF', dest='gtf', help='GTF file from Exonerate. Mandatory opt.', action='store', metavar='FILE')

(opts, args) = parser.parse_args()

mandatories = ['gtf']
for m in mandatories:
        if not opts.__dict__[m]:
                print "\nWARNING! One or more options not specified\n"
                parser.print_help()
                exit(-1)

############################## Reading files and parametersfrom sys import argv

file=open(opts.gtf, 'r')

OldTrToNew=defaultdict(list)
OldGeneToNew=defaultdict(list)
NewGeneToTr=defaultdict(list)

gtf_trs=defaultdict(list)
gtf_exo=defaultdict(list)
gtf_scf=defaultdict(list)

gtf=file.readline()
g=0
while gtf:
	if gtf.startswith('#'): gtf=file.readline()
	else:
		el=gtf.strip().split('\t')
		scf=el[0]
		if el[2] == 'gene':
			info=el[-1].split(';')
			for i in range(0,len(info)):
				if 'gene_id' in info[i]:
					g+=1
					Ngene='%sG%s' % (scf,g)
					gene=info[i].strip().split(' ')[1].replace('"','')
					OldGeneToNew[gene].append(Ngene)
		elif el[2] == 'transcript':
			if el[6] != '.':
				strand=el[6].strip()
				start=el[3]
				end=el[4]
				info=el[-1].split(';')
				for i in range(0,len(info)):
					if 'transcript_id' in info[i]:
						transid=info[i].strip().split(' ')[1].replace('"','')

				t=1
				Ntransid='%s.%s' % (Ngene,t)

				if Ngene not in NewGeneToTr:
					NewGeneToTr[Ngene].append(Ntransid)
				else:
					while Ntransid in NewGeneToTr[Ngene]:
						t+=1
						Ntransid='%s.%s' % (Ngene,t)
				NewGeneToTr[Ngene].append(Ntransid)

				OldTrToNew[transid].append(Ntransid)
				score='100'
				trInfo='%s\t%s\t%s\t%s\t%s\t%s' % (scf,start,end,strand,transid,score)
				gtf_scf[scf].append(transid)
				gtf_trs[transid].append(trInfo)
		gtf=file.readline()
		

file=open(opts.gtf, 'r')
gtf=file.readline()

while gtf:
	if gtf.startswith('#'): gtf=file.readline()
	else:
		el=gtf.strip().split('\t')
		if el[2] == 'exon':
			if el[6] != '.':
				info=el[-1].replace('-','_')
				tmp=info.split(';')
				strand=el[6]
				start=el[3]
				end=el[4]
				for i in range(0,len(tmp)):
					if 'transcript_id' in tmp[i]:
						transid=tmp[i].strip().split(' ')[1].replace('"','')
				gtf_exo[transid].append(gtf.strip())
		gtf=file.readline()

for scf in sorted(gtf_scf.keys()):
	for i in range(0,len(gtf_scf[scf])):
		lengths=[]
		starts=[]
		transid=gtf_scf[scf][i]
		tr=gtf_trs[transid][0].split('\t')
		trID=tr[4]
		trStart=int(tr[1])
		trEnd=int(tr[2])
		strand=tr[3]
		score=100
		NtrID=OldTrToNew[trID][0]
		transcript_info='%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t255,0,0' % (scf, trStart-1, trEnd, NtrID, score, strand, trStart-1, trEnd)
		for j in range(0,len(gtf_exo[transid])):
			start=int(gtf_exo[transid][j].split('\t')[3])
			end=int(gtf_exo[transid][j].split('\t')[4])
			lengths.append(str(end-start+1))
			intLen=start-trStart
			starts.append(str(intLen))
		nexons=str(len(lengths))
		if strand == '-':
			starts.reverse()
			lengths.reverse()
		print '%s\t%s\t%s\t%s' % (transcript_info,nexons,','.join(lengths),','.join(starts))

