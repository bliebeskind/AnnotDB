#! /usr/bin/env python

from peewee import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from itertools import izip_longest
import pickle
from trinity_parser import parse_trinity
from uniprot_blast_parser import tsv_line_gen
from parse_hmmscan import byline_scan_generator
	
db = SqliteDatabase(None,threadlocals=True)

class SqliteModel(Model):
	'''Base model defining which database to use'''
	class Meta:
		database = db
	
class Trinity(SqliteModel):
	'''
	Primary model containing data for each transcript. Fields are:
	name: CharField(primary_key=True)
	seq: BlobField() -- pickled BioPython TrinityIO object
	prot: BlobField() -- pickled BioPython TrinityIO object
	dna_length: IntegerField()
	prot_length: IntegerField()
	is_canonical: BooleanField()
	'''
	name = CharField(primary_key=True)
	seq = CharField()
	orf = CharField()
	prot = CharField()
	dna_len = IntegerField()
	orf_len = IntegerField()
	prot_len = IntegerField()
	is_canonical = BooleanField()
	num_var = IntegerField()
		
class Uniprot(SqliteModel):
	'''
	Model of simple Uniprot annotation. Includes three fields (columns):
	name: ForeignKeyField
	uniprot_id: CharField
	title: 		CharField
	'''
	name = ForeignKeyField(Trinity,related_name='uniprot_annotation') # connect to Trinity
	uniprot_id = CharField(default=None,null=True)
	title = CharField(default=None,null=True)
	### Add evalue ###
	
class PFAM(SqliteModel):
	'''
	Model for PFAM annotation using hmmscan. Six fields:
	name = ForeignKeyField
	accession = Charfield
	description = Charfield
	evalue = FloatField
	expected_domains = FloatField
	target = Charfield
	'''
	name = ForeignKeyField(Trinity,related_name='pfam_annotations') # connect to Trinity
	accession = CharField(default=None,null=True)
	description = CharField(default=None,null=True)
	evalue = FloatField(default=None,null=True)
	expected_domains = FloatField(default=None,null=True)
	target = CharField(default=None,null=True)

class TrinityDB:
	
	'''
	Class for loading a database with output from a Trinity assembly and
	annotations from Uniprot using Blastp and PFAM using hmmscan.
	'''
	
	def __init__(self,db_name):
		db.init(db_name)
		db.connect()
		
	def close(self):
		db.close()
	
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
		Trinity.create_table()
		with db.transaction():
			Trinity.insert_many(parse_trinity(infile,min_length)).execute()
			
	def get_canonicals(self,kind='prot'):
		'''Returns an iterator of SeqRecord objects corresponding to the 
		longest sequences for each comp. Kinds can be "seq" for the whole
		sequence; "orf" for the longest open reading frame; and "prot" for
		translated ORFs'''
		assert kind.lower() in ["seq","orf","prot"],\
		"Kind must be 'seq', 'orf', or 'prot'"
		kind = kind.lower()
		if kind == 'prot':
			for i in Trinity.select().where(Trinity.is_canonical == True):
				seq = SeqRecord(Seq(i.prot),id=i.name,description='')
				yield seq
		elif kind == 'orf':
			for i in Trinity.select().where(Trinity.is_canonical == True):
				seq = SeqRecord(Seq(i.orf),id=i.name,description='')
				yield seq
		elif kind == 'seq':
			for i in Trinity.select().where(Trinity.is_canonical == True):
				sequence = SeqRecord(Seq(i.seq),id=i.name,description='')
				yield sequence
		else:
			raise

	def write_canonicals(self,outfile,kind='prot',format='fasta'):
		'''Write canonical (longest) transcripts to outfile. Kinds can be 
		"seq" for the whole sequence; "orf" for the longest open reading 
		frame; and "prot" for translated ORFs.'''
		SeqIO.write(self.get_canonicals(kind),outfile,format)
		

	def load_uniprot(self,infile):
		'''Load top hits from .xml blast of Uniprot database.'''
		Uniprot.create_table()
		with db.transaction():
			fields = ["name","uniprot_id","title"]
			for row in tsv_line_gen(infile):
				row = row.strip().split("\t")
				data_dict = dict(izip_longest(fields,row))
				transcript_name = Trinity.get(Trinity.name == str(data_dict['name']))
				Uniprot.create(name=transcript_name,
					uniprot_id=data_dict['uniprot_id'],
					title=data_dict['title'])
				
	def load_pfam(self,infile):
		'''Load PFAM table from hmmscan table (must use --tblout option).'''
		PFAM.create_table()
		with db.transaction():
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
