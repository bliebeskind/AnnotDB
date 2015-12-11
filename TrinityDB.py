#! /usr/bin/env python

import time, csv
import sqlite3 as sql
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from trinity_parser import parse_trinity
from uniprot_blast_parser import tsv_line_gen
from parse_hmmscan import byline_scan_generator	

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
		
	def drop(self,table):
		'''Drop a table from the database. PERMANENTLY!'''
		cur = self.con.cursor()
		cur.execute('drop table %s' % table)
		self.con.commit()
		
	def table_names(self):
		'''Show tables'''
		tables = self.con.execute('''
			SELECT NAME FROM sqlite_master WHERE TYPE="table"''').fetchall()
		print "Found %i tables:" % len(tables)
		for t in tables:
			print '  ' + t[0]
			
	def table_info(self):
		'''Show tables and number of rows and columns in each table.'''
		print "\t".join(["Table","Columns","Rows"])
		tables = self.con.execute('''
			SELECT NAME FROM sqlite_master WHERE TYPE="table"''').fetchall()
		table_list = [t[0] for t in tables]
		for table in table_list:
			num_cols = len(self.con.execute('''
				PRAGMA table_info(%s)''' % table).fetchall())
			num_rows = self.con.execute('''
				SELECT Count() FROM %s''' % table).fetchone()[0]
			print "\t".join([table,str(num_cols),str(num_rows)])
	
	def load_trinity(self,infile,min_length=0):
		'''
		Load sequences output by a Trinity assembly into table "Trinity"
		Fields are:
	
		name: 		The unique name of the sequence given by Trinity (TEXT, P. KEY)
		seq:		TEXT of the full sequence
		orf:		TEXT of the longest ORF
		prot: 		TEXT of the translated longest ORF
		is_canonical:	Boolean - whether this is the longest transcript
		dna_len:	Length of full sequence (INTEGER)
		orf_len:	Length of longest ORF (INTEGER)
		prot_len:	Length of tranlated longest ORF (INTEGER)
		num_var:	Number of variants associated with this comp (INTEGER)
	
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

	def load_uniprot(self,infile):
		'''Load top hits from .xml blast of Uniprot database to table "Uniprot"
		Fields are:
		name: 	Trinity name (TEXT, PRIMARY KEY, FOREIGN KEY to Trinity.name)
		uniprot_id: Unique Uniprot identifier (TEXT)
		title: 	Uniprot description of Uniprot.uniprot_id (TEXT)
		'''
		try:
			self.con.execute('''
				CREATE TABLE
				Uniprot
				(name TEXT PRIMARY KEY,
				uniprot_id TEXT,
				title TEXT,
				evalue REAL,
				FOREIGN KEY(name) REFERENCES Trinity(name))
				''')
		except sql.Error as e: # catch if table already exists
			return e # Could make this more sophisticated, user input?
		with self.con:
			self.con.executemany("INSERT INTO Uniprot VALUES (?,?,?,?)", \
				tsv_line_gen(infile))
				
	def load_pfam(self,infile):
		'''
		Load PFAM table from hmmscan table (must use --tblout option).
		Fields are:
		target: 	The PFAM domain hit (TEXT)
		accession:	PFAM accession (TEXT)
		name:		Trinity name (TEXT, FOREIGN KEY to Trinity.name)
		evalue:		(REAL)
		expected_domains: The expected number of PFAM.target domains (REAL)
		description:	Description of PFAM.target (TEXT)
		'''
		try:
			self.con.execute('''
				CREATE TABLE
				PFAM
				(target TEXT,
				accession TEXT,
				name TEXT,
				evalue REAL,
				expected_domains REAL,
				description TEXT,
				FOREIGN KEY(name) REFERENCES Trinity(name))
				''')
		except sql.Error as e: # catch if table already exists
			return e # Could make this more sophisticated, user input?
		with self.con:
			self.con.executemany("INSERT INTO PFAM VALUES (?,?,?,?,?,?)", \
				byline_scan_generator(infile))
				
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
		
	def load_user_table(self,csv_file,columns,table_name,dialect='excel'):
		'''
		Load a user defined table from a csv file. 
		
		Parameters:
		csv_file: Flat text infile. Header must match column dict.
		columns: Dictionary mapping column names to type. Must have column
			"name" which matches "name" column from Trinity table.
			Example: {'name':'text','test1':'integer','test2':'integer'}
		table: The name of your table. Must be one word.
		dialect: Optional parameter for reading infile. See python's csv.reader
			for documentation.
			
		**Note: Does not read xlsx. Please export your Excel file to csv.
		'''
		with open(csv_file,'rb') as f:
			reader = csv.reader(f,dialect)
			# Create sql columns and values to be substituted:
			# cols: ({},{}...FOREIGN KEY(name) REFERENCES Trinity(name)))
			# vals: (?,?,?...)
			header_line = reader.next()
			assert len(header_line) == len(columns), \
			"Number of columns don't match table"
			temp = "(" + ",".join(("{}" for i in header_line))
			cols_string = temp + ",FOREIGN KEY(name) REFERENCES Trinity(name))"
			vals_string = "(" + ",".join(("?" for i in header_line)) + ")"
			
			# Create substitution: (column1 integer, column2 text...)
			func = lambda x: ' '.join((x,columns[x]))
			try: # create column names and types from infile and column dict
				sub_tuple = tuple(map(func,header_line))
			except KeyError, e: # header line and column dict don't match
				raise Exception("Bad column from infile: %s" % e)
			
			# Create table
			try:
				sql_string = "CREATE TABLE {} {}".format(table_name,cols_string)
				self.con.execute(sql_string.format(*sub_tuple))
			except sql.Error as e: # catch if table already exists
				return e # Could make this more sophisticated, user input?
			
			# Populate table
			csv_gen = (row for row in reader)
			with self.con:
				self.con.executemany('''
					INSERT INTO {} VALUES {}
					'''.format(table_name,vals_string), csv_gen)
		
	def _assert_kind(self,kind):
		L = ["prot","orf","seq"]
		try:
			assert kind in L
		except AssertionError:
			raise Exception("Kind must be 'prot','orf', or 'seq': %s" % kind)
		
	def get_canonicals(self,cutoff=0,seq_type='prot'):
		'''Returns an iterator of SeqRecord objects corresponding to the 
		longest sequences for each comp. Kinds can be "seq" for the whole
		sequence; "orf" for the longest open reading frame; and "prot" for
		translated ORFs'''
		self._assert_kind(seq_type)
		selection = self.con.execute('''SELECT name, %s FROM Trinity \
			WHERE Trinity.is_canonical==1 AND Trinity.prot_len >= %s''' % (seq_type,cutoff))
		count = 0
		for name,prot in selection:
			seq = SeqRecord(Seq(prot),id=name,description='')
			count +=1
			if count % 100 == 0:
				print str(count)
			yield seq

	def write_canonicals(self,outfile,cutoff=0,seq_type='prot',format='fasta'):
		'''Write canonical (longest) transcripts to outfile. Kinds can be 
		"seq" for the whole sequence; "orf" for the longest open reading 
		frame; and "prot" for translated ORFs.'''
		SeqIO.write(self.get_canonicals(cutoff,seq_type),outfile,format)
		
	def uniprot_descriptions_from_domain(self,domain,delimiter='\t'):
		'''
		Returns a generator of sequence names and Uniprot descriptions
		for all rows that match a PFAM domain.
		'''
		search = self.con.execute('''
			SELECT t.name,u.title
			FROM Trinity t JOIN PFAM p ON t.name=p.name 
			JOIN Uniprot u ON t.name=u.name
			WHERE p.target=?''', (domain,))
		return (delimiter.join([row[0],row[1]]) for row in search.fetchall())
	
	def to_fasta_from_domain(self,domain,outfile,seq_type='prot',evalue=10):
		'''
		Search for all sequences that match a PFAM domain and write them to
		a fasta file. seq_type must be "prot" "orf" or "seq"
		'''
		self._assert_kind(seq_type)
		c = self.con.cursor()
		search = c.execute('''
			SELECT t.name,{},u.title 
			FROM TRINITY t JOIN PFAM p ON t.name=p.name
			JOIN Uniprot u ON t.name=u.name
			WHERE p.target=?
			AND p.evalue <= ?'''.format("t."+seq_type), (domain,evalue))
		seq_gen = (SeqRecord(Seq(j),id=i,description=k) for i,j,k in search)
		SeqIO.write(seq_gen,outfile,'fasta')
		
	def to_fasta_from_blast(self,string,outfile,seq_type='prot',evalue=10):
		'''
		Search for sequences with <string> in the Blast hit and write to
		fasta file.
		'''
		self._assert_kind(seq_type)
		search = self.con.execute('''
			SELECT t.name, {},u.title
			FROM TRINITY t JOIN Uniprot u ON t.name=u.name
			WHERE u.title LIKE ?
			AND u.evalue <= ?'''.format("t."+seq_type),('%'+string+'%',evalue))
		seq_gen = (SeqRecord(Seq(j),id=i,description=k) for i,j,k in search)
		SeqIO.write(seq_gen,outfile,'fasta')
