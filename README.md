AnnotDB
===========

This module allows you to upload transcriptomes and attendent annotations
into a SQLite database. If you find yourself making transcriptomes frequently
and would like to keep them organized, then this module might help. It also 
allows you to share your transcriptome with GUI-only users without having to 
set up a blast server. Simply load up your transcriptome and annotations, and 
send it to them. They can then use one of the easy and free SQLite browers to
search for genes of interest.

It has built in support for annotations from BLAST and HMMER's hmmscan program,
but assumes that these will be run against Uniprot and PFAM, respectively. 
There is also support for user-defined annotations, which must be put in a csv 
file. 

There are a few basic functions to help you search the database once you've
loaded it, and two SQL scripts which GUI users can use to search using, for
instance, [SQLiteBrowser] (http://sqlitebrowser.org/).

AnnotDB currently only has support for Trinity transcriptomes.
Requires BioPython.

## TrinityDB ##

TrinityDB will upload a  Trinity assembly and associated annotations into a
SQLite database. There are also functions for simple but helpful stuff
like writing out canonical (longest) transcripts and associated proteins
using BioPython's SeqIO methods.

TrinityDB divides the sequences and annotations into three tables. These are
listed here along with their associated fields.

* **Trinity** - Main table of sequences assembled by Trinity
  * **name**: Name of the sequence, e.g. "comp29167_c0_seq1"
  * **seq**: Full nucleotide sequence
  * **orf**: Longest open reading frame (only forward frame are checked)
  * **prot**: Protein associated with longest ORF
  * **dna\_len**: Length of full sequence
  * **orf\_len**: Length of longest ORF
  * **prot\_len**: Length of longest translation
  * **is_canonical**: Boolean, whether this is the longest splice variant
  * **num\_var**: Number of variants associated with this comp
  

* **Uniprot** - Uniprot annotation using Blastp. Input must be XML (outfmt 5)
  * **name**: Name of the sequence, Foreign key from Trinity
  * **uniprot_id**: Uniprot id of the top Blastp hit
  * **title**: Uniprot name of the top Blastp hit
  * **evalue**: E-value of the hit

  
* **PFAM** - PFAM annotation using hmmscan. Input must be table format (--tblout)
  * **name**: Name of the sequence, Foreign key from Trinity
  * **accession**: PFAM accession of target
  * **description**: Description of target
  * **evalue**: e-value of the target hit
  * **expected\_domains**: Expected number of this target in query
  * **target**: The PFAM domain hit
 
## Quick Start ##

```python
import TrinityDB
myTinyDb = TrinityDB.TrinityDB("myDb.db") # initialize database
myTinyDb.load_all("MyTrinity.fas","MyBlast.xml","MyHmmscan.txt") # load all annotations
#  Loading from Trinity fasta file
#  Loading Uniprot BLAST annotations
#  Loading PFAM hmmscan annotations
myTinyDB.table_info()
#  Table Columns Rows
#  Trinity	9	11
#  Uniprot	4	10
#  PFAM		6	8
```

## Load a user-defined table ##
This is intended to be flexible. Your table could easily be transcript counts
from various treatments or any other set of values you can put in a table.
The first column should always be called "name" and have the transcript names
from the original Trinity outfile.

```python
### Table should be something like the following
#
### name, 			 treatment1, 	treatment2
### comp297_c0_seq1	 2,				39
### comp491_c0_seq1	 6,				93
#
### The table must be in a flat text file, not Excel
#
### User must also provide a dictionary mapping
### the column names to their intended SQL types.
#
cols = {'name':'text', 'treatment1':'integer'...}
#
myTinyDB.load_user_table("myTable.csv",cols,"userTable",dialect='excel')
#
### You can now retrieve values from this or any other table based on values
### in any of those tables.
```

## Write out longest transcripts for each gene ##

```python
## Defaults: length cutoff:0 - sequence type:'prot' - format: fasta
myTinyDB.write_canonicals("My_canonical_prots.fas")

## Write ORFs to a nexus file
myTinyDB.write_canonicals("My_100bp_ORFs.nex",seq_type='orf',format='nexus')
```

## Write out Uniprot descriptions of all sequences matching a PFAM domain ##

```python
with open("myDescriptions.txt",'a') as f:
	search_generator = myTinyDB.uniprot_descriptions_from_domain("SomeDomain")
	for line in search_generator:
		f.write(line)
```

## Write out all sequences matching a PFAM domain ##

```python
## Defaults: sequence type: Protein - evalue: 10
myTinyDB.to_fasta_from_domain("SomeDomain","myFasta.fas")

## Write ORFs
myTinyDB.to_fasta_from_domain("SomeOtherDomain","myORFFasta.fas",seq_type='orf')

## Only write hit with e-value equal or better than .01
myTinyDB.to_fasta_from_domain("SomeDomain","myFasta.fas",evalue=.01)
```

## Write out sequences with a key word in the blast annotation

```python
## Defaults: sequence type: Protein - evalue: 10
## Will write all proteins with "channel" in the top blast hit
myTinyDB.to_fasta_from_blast("channel","MyBlastFasta.fas")
```
