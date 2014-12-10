AnnotDB
===========

Modules for working with protein annotation databases using sqlite3.
Requires BioPython.

## TrinityDB ##

TrinityDB will upload a  Trinity assembly and associated annotations into a
SQLite database. There are also functions for simple but helpful stuff
like writing out canonical (longest) transcripts and associated proteins
using BioPython's SeqIO methods.

There are three tables listed here along with their associated fields.

* **Trinity** - Main model of sequences assembled by Trinity
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
## Defaults: sequence type: Protein
myTinyDB.to_fasta_from_domain("SomeDomain","myFasta.fas")

## Write ORFs
myTinyDB.to_fasta_from_domain("SomeOtherDomain","myORFFasta.fas",seq_type='orf')
```
