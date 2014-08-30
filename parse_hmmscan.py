#! /usr/bin/env python

## For parsing hmmscan table generated with --tblout option.

## TO DO: find a way to dump information from scan_generator to database

import sys
from Bio import SeqIO

def main_generator(hmmscan_table):
	'''Main line generator. Yields list of field values'''
	with open(hmmscan_table) as f:
		for line in f:
			if not line.startswith("#"):
				line = line.split()
				target_descr = ' '.join(line[18:])
				yield line[:18] + [target_descr]
				
def scan_generator(hmmscan_table):
	'''
	Yield dictionaries corresponding to one query at a time. Format is
	the query as key mapped to a list of dictionaries corresponding to the 
	field values:
	{Query1:
		[{target_name: hmm1, accession: PF0001, evalue: 0}, # First hit
		{target_name: hmm2, accession: PF0002, evalue: .000001}]
	}
	'''
	scan = {}
	for line in main_generator(hmmscan_table):
		query = line[2]
		hit_D = {"target_name": line[0], 
			"accession": line[1], "evalue": float(line[4])}
		if scan == {}:
			scan[query] = [hit_D]
		elif query in scan:
			scan[query].append(hit_D)
		else:
			yield scan
			scan = {}
			scan[query] = [hit_D]
	yield scan
	
def byline_scan_generator(hmmscan_table):
	'''
	Yield dictionaries corresponding to one query at a time. Format is
	the query as key mapped to a list of dictionaries corresponding to the 
	field values:
	
	{query: name, target_name: hmm1, accession: PF0001, evalue: 0}
	'''
	fields = ["target","accession","name","evalue","expected_domains","description"]
	for line in main_generator(hmmscan_table):
		values = [line[0],line[1],line[2],float(line[4]),float(line[10]),line[-1]]
		yield dict(zip(fields,values))

def matches(hmmscan_table,domain,thresh):
	'''Generator function. Returns stream of query names that hit <domain>
	with an evalue below <thresh>'''
	for scan in scan_generator(hmmscan_table):
		assert len(scan) == 1, "Multiple queries in scan: %s" % \
		('\n'.join(scan.keys()))
		query = scan.keys()[0]
		for hit in scan[query]:
			if hit["target_name"] == domain \
			and hit["evalue"] <= float(thresh):
				yield query
				break
		else:
			continue
			
def filter_seqs(infile,hmmscan_table,domain,thresh,format='fasta'):
	'''
	Return generator of SeqRecords that matched <domain> below evalue
	<thresh> in <hmmscan_table>.
	Can be used like so:
	
	>>> filtered_seq_gen = filter_seqs(infile,hmmscan_table,etc...)
	>>> SeqIO.write(filtered_seq_gen, "my_matches.fas", "fasta")
	'''
	match_list = [m for m in matches(hmmscan_table,domain,thresh)]
	records = SeqIO.parse(infile,format)
	total_count = 0
	filtered_count = 0
	for rec in records:
		total_count +=1
		if rec.id in match_list:
			filtered_count +=1
			yield rec
	sys.stderr.write("%i out of %i sequences had domain: %s\n" % \
		(filtered_count,total_count,domain))
