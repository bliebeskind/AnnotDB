from peewee import *
from Models import *
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def uniprot_by_domain(domain):
	'''
	Searches for Uniprot annotations of protein match Pfam domain.
	Returns an iterators of csv rows for two columns:
	Sequence name and Uniprot description.
	
	Usage:
	>>> select_iter = uniprot_by_domain("7tm_1")
	>>> with open("oufile.csv",'a') as f:
	... 	for row in selec_iter:
	...			f.write(row)
	'''
	selection = (Trinity
				.select(Trinity,
					Uniprot.title.alias("title"))
				.join(PFAM)
				.switch(Trinity)
				.join(Uniprot)
				.where(PFAM.target == domain)
				.naive()
				.iterator())
	yield ",".join(["Sequence","Uniprot_name"]) + "\n"
	for row in selection:
		yield ",".join([row.name,row.title]) + "\n"
		
def sequence_by_domain(domain,kind='prot'):
	'''
	Get sequences based on domain annotation. Returns an iterator of
	Biopython SeqRecord objects.
	
	Usage:
	>>> SeqIO.write(sequence_by_domain("V1R"),"V1Rs.fas",'fasta')
	'''
	if kind == 'prot':
		expr = "row.prot"
	elif kind == 'orf':
		expr = "row.orf"
	elif kind == 'seq':
		expr = "row.seq"
	else:
		raise Exception("Kind %s not found" % kind)
	selection = (Trinity
				.select(Trinity)
				.join(PFAM)
				.where(PFAM.target == domain)
				.naive()
				.iterator())
	for row in selection:
		outseq = SeqRecord(Seq(eval(expr)),id=row.name,description='')
		yield outseq
