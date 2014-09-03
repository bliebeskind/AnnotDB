AnnotDB
===========

Modules for working with protein annotation databases using sqlite and peewee.
Also requires BioPython

###TrinityDB

TrinitDB provides peewee models for a  Trinity assembly and associated 
annotations. There are functions for loading a database, and simple stuff
like writing out canonical (longest) transcripts and associated proteins
using BioPython's SeqIO methods.

The three models are listed here, along with their associated fields.

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
  

* **Uniprot** - Uniprot annotation using Blastp
  * **name**: Name of the sequence, Foreign key from Trinity
  * **uniprot_id**: Uniprot id of the top Blastp hit
  * **title**: Uniprot name of the top Blastp hit

  
* **PFAM** - PFAM annotation using hmmscan. Multiple rows per query.
  * **name**: Name of the sequence, Foreign key from Trinity
  * **accession**: PFAM accession of target
  * **description**: Description of target
  * **evalue**: e-value of the target hit
  * **expected\_domains**: Expected number of this target in query
  * **target**: The PFAM domain hit
