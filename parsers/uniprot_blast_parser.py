#! /usr/bin/env python

from Bio.Blast import NCBIXML
import sys

def hit_gen(infile):
	'''Given blast xml output (--outfmt 5), return generator of queries and 
	top hit titles.'''
	with open(infile) as f:
		records = NCBIXML.parse(f)
		for rec in records:
			try:
				yield rec.query, rec.alignments[0].title, rec.descriptions[0].e
			except IndexError: # no hits below evalue threshold
				yield rec.query, None, None
			
def tsv_line_gen(infile):
	'''Generator of parsed Blast XML output from a Uniprot search. Yields 
	tuples of "query", "Uniprot ID", and long Uniprot title.'''
	count = 0
	for query,top_hit,evalue in hit_gen(infile):
		if top_hit != None:
			descr = top_hit.strip().split(" ")
			uni = descr[1].split("|")
			assert uni[0] == 'sp' or uni[0] == 'tr', \
			"Unknown database: %s" % uni[0]
			uniprot_id = uni[1]
			title = ' '.join(descr[2:])
		else:
			uniprot_id, title = '',''
		yield (query,uniprot_id,title,evalue)
		count +=1
		if count % 100 == 0:
			print str(count)
		
if __name__ == '__main__':
	infile = sys.argv[1]
	print "\t".join(["Query","Uniprot Id", "Title"])
	for line in tsv_line_gen(infile):
		print '\t'.join(line)
