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
				yield rec.query, rec.alignments[0].title
			except IndexError: # no hits below evalue threshold
				yield rec.query, None
			
def tsv_line_gen(infile):
	count = 0
	for query,top_hit in hit_gen(infile):
		if top_hit != None:
			descr = top_hit.strip().split(" ")
			uni = descr[1].split("|")
			assert uni[0] == 'sp' or uni[0] == 'tr', \
			"Unknown database: %s" % uni[0]
			uniprot_id = uni[1]
			name = ' '.join(descr[2:])
		else:
			uniprot_id, name = '',''
		yield "\t".join([query,uniprot_id,name])
		count +=1
		if count % 100 == 0:
			print str(count)
		
if __name__ == '__main__':
	infile = sys.argv[1]
	print "\t".join(["Query","Uniprot Id", "Title"])
	for line in tsv_line_gen(infile):
		print line
