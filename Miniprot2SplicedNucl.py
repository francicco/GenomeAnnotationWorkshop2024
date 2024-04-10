#!/usr/bin/python

import optparse
from collections import defaultdict

################################# Command line options

desc='Converte a miniprot GFF into an Exonerate GFF file'

parser = optparse.OptionParser(description=desc, version='%prog version 0.1 - 25-09-2016 - Author: FCicconardi')

parser.add_option('-g', '--GFF', dest='gff', help='GFF file from Miniprot. Mandatory opt.', action='store', metavar='FILE')

(opts, args) = parser.parse_args()

mandatories = ['gff']
for m in mandatories:
        if not opts.__dict__[m]:
                print "\nWARNING! One or more options not specified\n"
                parser.print_help()
                exit(-1)

############################## Reading files and parametersfrom sys import argv

#print '##gff-version 2'

file=open(opts.gff, 'r')

line=file.readline()

while line:
	if line.startswith('##PAF'):
		ID=line.split()[1]
	elif line.startswith('##ATN'):
		CDS = "".join(filter(str.isupper, line.split()[1])).replace('-','')
	elif line.startswith('##STA'):
		print '>%s\n%s' % (ID,CDS)
	line=file.readline()
