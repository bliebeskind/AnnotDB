#! /usr/bin/env python
import time

import sqlite3 as sql
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from itertools import izip_longest
#import pickle
from trinity_parser import parse_trinity
from uniprot_blast_parser import tsv_line_gen
from parse_hmmscan import byline_scan_generator
from Models import *
	


class TrinityDB:	
	'''
	Class for loading a database with output from a Trinity assembly and
	annotations from Uniprot using Blastp and PFAM using hmmscan.
	'''
	
	def __init__(self,db_name):
		self.con = sql.connect(db_name)
		
	def close(self):
		'''Close database'''
		self.con.close()
	
	def load_trinity(self,infile,min_length=0):
		'''
		Load sequences output by a Trinity assembly into table: "Trinity
		Dictionary fields are:
	
		name: 		The unique name of the sequence given by Trinity (str)
		seq:		String of the full sequence
		orf:		String of the longest ORF
		prot: 		String of the translated longest ORF
		is_canonical:	Boolean - whether this is the longest transcript
		dna_len:	Length of full sequence
		orf_len:	Length of longest ORF
		prot_len:	Length of tranlated longest ORF
		num_var:	Number of variants associated with this comp
	
		Specifying a minimum length will skip sequences that don't have a 
		tranlatable region above the specified number of amino acids.
	
		**Note: only the forward frames are searched for protein sequences."
		'''
		t0 = time.time()
		try:
			self.con.execute('''
				CREATE TABLE
				Trinity
				(name TEXT PRIMARY KEY,
				seq TEXT,
				orf TEXT,
				prot TEXT,
				is_canonical INTEGER,
				dna_len INTEGER,
				orf_len INTEGER,
				prot_len INTEGER,
				num_var INTEGER)
				''')
		except sql.Error as e: # catch if table already exists
			return e # Could make this more sophisticated, user input?
		with self.con:
			self.con.executemany("INSERT INTO Trinity VALUES (?,?,?,?,?,?,?,?,?)", \
				parse_trinity(infile,min_length))
		t1 = time.time()
		print "Time taken: %i" % (t1 - t0)
			
	def get_canonicals(self,cutoff=0,kind='prot'):
		'''Returns an iterator of SeqRecord objects corresponding to the 
		longest sequences for each comp. Kinds can be "seq" for the whole
		sequence; "orf" for the longest open reading frame; and "prot" for
		translated ORFs'''
		L = ["prot","orf","seq"]
		try:
			assert kind in L
		except AssertionError:
			raise Exception("Kind must be 'prot','orf', or 'seq': %s" % kind)
		selection = self.con.execute('''SELECT name, %s FROM Trinity \
			WHERE Trinity.is_canonical==1 AND Trinity.prot_len >= %s''' % (kind,cutoff))
		count = 0
		for name,prot in selection:
			seq = SeqRecord(Seq(prot),id=name,description='')
			count +=1
			if count % 100 == 0:
				print str(count)
			yield seq

	def write_canonicals(self,outfile,cutoff=0,kind='prot',format='fasta'):
		'''Write canonical (longest) transcripts to outfile. Kinds can be 
		"seq" for the whole sequence; "orf" for the longest open reading 
		frame; and "prot" for translated ORFs.'''
		SeqIO.write(self.get_canonicals(cutoff,kind),outfile,format)
		

	def load_uniprot(self,infile):
		'''Load top hits from .xml blast of Uniprot database.'''
		try:
			self.con.execute('''
				CREATE TABLE
				Uniprot
				(name TEXT PRIMARY KEY,
				uniprot_id TEXT,
				title TEXT,
				FOREIGN KEY(name) REFERENCES Trinity(name))
				''')
		except sql.Error as e: # catch if table already exists
			return e # Could make this more sophisticated, user input?
		with self.con:
			self.con.executemany("INSERT INTO Uniprot VALUES (?,?,?)", \
				tsv_line_gen(infile))
				
	def load_pfam(self,infile):
		'''Load PFAM table from hmmscan table (must use --tblout option).'''
		PFAM.create_table()
		with self.database.transaction():
			for row in byline_scan_generator(infile):
				transcript_name = Trinity.get(Trinity.name == str(row['name']))
				PFAM.create(name=transcript_name,
					description=row['description'],
					evalue=row['evalue'],
					expected_domains=row['expected_domains'],
					target=row['target'])
				
	def load_all(self,trinity_infile,uniprot_xml,hmmscan_table,min_length=0):
		'''
		Load trinity fasta file and annotations. Uses BLASTp annotations
		against Uniprot and hmmscan annotations agains PFAM. Blast output
		must be in xml format (outfmt 5), and hmmscan output must be table
		format (--tblout option)
		'''
		print "Loading from Trinity fasta file"
		self.load_trinity(trinity_infile,min_length)
		print "Loading Uniprot BLAST annotations"
		self.load_uniprot(uniprot_xml)
		print "Loading PFAM hmmscan annotations"
		self.load_pfam(hmmscan_table)
