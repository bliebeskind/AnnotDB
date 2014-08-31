from peewee import *

db_proxy = Proxy()

class SqliteModel(Model):
	'''Base model defining which database to use'''
	class Meta:
		database = db_proxy
	
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
