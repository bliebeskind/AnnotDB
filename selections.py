from peewee import *
from Models import *

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
