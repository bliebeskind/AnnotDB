#! /usr/bin/env python

### For packing values into AnnotDB.Seq

from Bio import SeqIO
from copy import copy
import pickle
	
def frames(rec, min_length=0):
	'''Given a Bio.Seq object, find forward reading frames and return Dictionary
	of structure: {orf_length: (start, end)}. Default minimum protein
	length if 100 AAs. Removes stop codons.'''
	reading_frames = {}
	seq_len = len(rec.seq)
	for frame in range(3):
		seq = rec.seq[frame:]
		seq = seq[:len(seq) - len(seq) %3]
		trans = str(seq.translate())
		trans_len = len(trans)
		aa_start, aa_end = 0, 0
		while aa_start < trans_len:
			aa_end = trans.find("*", aa_start)
			if aa_end == -1: # if "*" not found
				aa_end = trans_len
			if aa_end - aa_start >= min_length:
				start = frame+aa_start*3
				end = min(seq_len, frame+aa_end*3+3) - 3 #remove stop codon
				reading_frames[(end - start)] = (start, end)
			aa_start = aa_end + 1 #start again after "*"
	return reading_frames
	
def get_orf(rec, min_length=0):
	'''Calls frames() and returns the nucleotide sequence of the longest orf.'''
	D = frames(rec, min_length)
	top_orf = sorted(D.keys())[-1]
	sequence = rec.seq[D[top_orf][0]:D[top_orf][1]]
	return sequence
	
def get_orfs_from_file(infile, min_length=0):
	'''Open fasta infile and return iterator of SeqRecords with ORF sequences.'''
	records = SeqIO.parse(infile, 'fasta')
	for rec in records:
		try:
			orf = copy(rec)
			orf.seq = get_orf(rec, min_length)
		except IndexError: # frames() found nothing above min length
			continue
		yield rec, orf
		
def return_previous(seq_list,orf_list,splices=False):
	'''Handle dumping of orfs in find_canonicals'''
	if not splices: 	# if no splice variants
		assert len(seq_list) == 1
		yield seq_list[0], orf_list[0], True, 1 # rec, orf, is_longest, num_var
	else:
		longest = reduce(lambda x,y: x if len(x)>= len(y) else y, orf_list).id
		for i,j in enumerate(orf_list):
			yield seq_list[i], j, j.id == longest, len(orf_list)
		
def find_canonicals(infile,min_length=0):
	'''
	Parse trinity output file and return generator that yields four-part tuples:
	
	sequence (SeqRecord)
	ORF (SeqRecord)
	is_longest_transcript (Boolean)
	num_variants
	'''
	orf_D = {} 	# dictionary for identifying longest transcripts
	seq_D = {}	# dictionary holding full sequences
	comp_list = [] 	# list for checking whether comps are contiguous
	for rec, orf in get_orfs_from_file(infile, min_length):
		comp = orf.id.split('_')[0]
		if orf_D == {}: # is first
			orf_D[comp] = [orf]
			seq_D[comp] = [rec]
			continue
		elif comp in orf_D:	# is a variant of previous sequence
			orf_D[comp].append(orf)
			seq_D[comp].append(rec)
		elif comp not in orf_D and orf_D != {}: # new gene, yield previous
			assert len(orf_D) == 1 and len(seq_D) == 1 #should have just one key
			assert comp not in comp_list, "Comps not contiguous"
			seq_list = seq_D[seq_D.keys()[0]]
			orf_list = orf_D[orf_D.keys()[0]]
			if len(orf_list) == 1: 	# if no splice variants
				assert len(seq_list) == 1
				yield seq_list[0], orf_list[0], True, 1 # rec, orf, is_longest, num_var
			else: # yield all, with info on whether it's longest
				longest = reduce(lambda x,y: x if len(x)>= len(y) else y, orf_list).id
				for i,j in enumerate(orf_list):
					yield seq_list[i], j, j.id == longest, len(orf_list)
			seq_D, orf_D = {},{}		# flush
			orf_D[comp] = [orf]
			seq_D[comp] = [rec]
		else:
			raise
		comp_list.append(comp)
	seq_list = seq_D[seq_D.keys()[0]] # repeats line 84-92 above, but oh well.
	orf_list = orf_D[orf_D.keys()[0]]
	if len(orf_list) == 1: 	
		assert len(seq_list) == 1
		yield seq_list[0], orf_list[0], True, 1 # rec, orf, is_longest, num_var
	else: # yield all, with info on whether it's longest
		longest = reduce(lambda x,y: x if len(x)>= len(y) else y, orf_list).id
		for i,j in enumerate(orf_list):
			yield seq_list[i], j, j.id == longest, len(orf_list)
			
def parse_trinity(infile,min_length=0):
	'''
	Parse a Trinity output file and return a generator yielding dictionary of 
	info about each sequence. Dictionary fields are:
	
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
	
	**Note: only the forward frames are searched for protein sequences.
	'''
	fields = ["name","seq","orf","prot","is_canonical","dna_len",
		"orf_len","prot_len","num_var"]
	count = 0
	for rec,orf,canon,isoforms in find_canonicals(infile,min_length):
		assert rec.id == orf.id
		prot = copy(orf)
		prot.seq = str(orf.seq.translate())
		info = [rec.id, str(rec.seq), str(orf.seq), str(prot.seq), canon, len(rec.seq), 
			len(orf.seq), len(prot.seq), isoforms]
		yield dict(zip(fields,info))
		count +=1
		if count % 100 == 0:
			print str(count)
