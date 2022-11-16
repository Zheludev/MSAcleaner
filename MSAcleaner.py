## the purpose of this script is to aggressively clean an input (fasta-format) MSA given a set of reference sequences
## the premise is no doubt evolutionarily flawed:
	## the idea is that given a reference set, remove any sequences that produce weakly supported insertions in the references
	## the output is then a fasta file (NOT MSA) that needs to be re-aligned

## the usecase is when you have a reference(s) and you only care about homologues to that reference with minimal innovations
	## though deletions are ok

## overview:
	## 1) - load MSA
		## fix to one-line
		## check each seq is same length
	## 2) - load reference seqIDs
	## 3) - make a seqID alias dict to deal with combersome seqIDs
	## 4) - load MSA into dataframe with seqs as rows and positions as columns
	## 5) - for each col (posn) in the reference sequences, check to see if there's an insertion over all references
	## 6) - if reference insertion found, check abundance in all other seqs - if below threshold, flag any seqID that has insertion for removal
	## 7) - iterate 5) and 6) over whole sequence length
	## 8) - remove the flagged sequences
	## 9) - redefine seqIDs back to initial seqIDs
	## 10) - convert remaining MSA to regular FASTA and save

## written by INZ - 11/15/22
## Stanford Unversity
## provided with no acceptance of liability or promise of functionality
## version 0.1.0

## load libraries

import pandas as pd
import numpy as np
from tqdm import tqdm
import argparse

def main(input_msa_fn,input_ref_fn,min_fraction,output_seq_fn):
	
	## set variables
	
	#input_msa_fn = 'iter_1.msa'
	#input_ref_fn = 'reference.seqIDs'
	#min_fraction = 0.05
	#output_seq_fn = 'cleaned.fa'
	
	## open input files and cleanup
	
	def opener(filename):
		open_file = open(filename, mode='r')
		open_list = open_file.read().splitlines()
		for entry_ind, entry in enumerate(open_list):
			open_list[entry_ind] = entry
		##
		return open_list
	##
	
	def singleline(in_fasta_list):
		out_fasta_list = []
		new_seq_line = []
		for line in in_fasta_list:
			if line.startswith('>'):
				if new_seq_line:
					out_fasta_list.append(''.join(new_seq_line))
					new_seq_line = []
				out_fasta_list.append(line)
			else:
				new_seq_line.append(line.strip())
			##
		if new_seq_line:
			out_fasta_list.append(''.join(new_seq_line))
		##
		return out_fasta_list
	##
	
	input_msa_ls = opener(input_msa_fn)
	input_msa_ls = singleline(input_msa_ls)
	
	input_ref_ls = opener(input_ref_fn)
	
	## check that the input MSA is the right shape (all seqs same length)
		## compare to first seq
	
	ref_len = len(input_msa_ls[1])
	
	for msa_seq_ind in range(1,len(input_msa_ls),2):
		msa_seq_len = len(input_msa_ls[msa_seq_ind])
		if (msa_seq_len != ref_len):
			print('MSA sequences are not all same length, please re-compute the MSA')
			print('quitting.')
			quit()
	##
	
	## make alias dict
		## redo the reference seqIDs to give them different names
	
	alias_dict = {}
	
	for msa_seqID_ind in range(0,len(input_msa_ls),2):
		msa_alias = 'msa_' + str(int((msa_seqID_ind/2)))
		alias_dict[input_msa_ls[msa_seqID_ind]] = msa_alias
	##
	
	for ref_seqID_ind, ref_seqID in enumerate(input_ref_ls):
		ref_alias = 'ref_' + str(ref_seqID_ind)
		alias_dict[ref_seqID] = ref_alias
	##
	
	## rename the input MSA based on the alias dict
	
	for seqID_ind in range(0,len(input_msa_ls),2):
		input_msa_ls[seqID_ind] = alias_dict[input_msa_ls[seqID_ind]]
	##
	
	## load MSA into dataframe
		## via numpy, make a 2 column matrix
		## then convert to dataframe
		## then sort by seqID (index)
		## then split seq
	
	input_msa_np = np.array(input_msa_ls)
	input_msa_np = input_msa_np.reshape(-1, 2)
	
	input_msa_df = pd.DataFrame(input_msa_np , columns = ['seqID','seq'])
	input_msa_df = input_msa_df.set_index('seqID')
	
	input_msa_df.sort_index(ascending = False, inplace = True)
	
	input_msa_df = input_msa_df['seq'].str.split(pat = '', expand = True).add_prefix('pos_')
	
	## clean up empty columns
	
	nan_value = float("NaN")
	input_msa_df.replace("", nan_value, inplace=True)
	input_msa_df.dropna(how='all', axis=1, inplace=True)
	
	## identify which positions are unsupported by the reference sequences
		## where number of '-' matches number of reference sequences
			## cant work out how to summarize with a char count so will have to iterate
	
	ref_msa_df = input_msa_df.filter(like= 'ref_' , axis = 0)
	
	ref_count = len(ref_msa_df)
	unsupported_positions = []
	
	for position in ref_msa_df:
		insertion_count = ref_msa_df[position].str.count('-').sum()
		if (ref_count == insertion_count):
			unsupported_positions.append(position)
	##
	
	## iterate over the unsupported positions and check if they're inusufficiently supported
	
	que_msa_df = input_msa_df.filter(like= 'msa_' , axis = 0)
	
	que_count = len(que_msa_df)
	
	seqID_omit_set = set()

	threshold = min_fraction * que_count
		
	for position in tqdm(unsupported_positions):
		insertion_count = (que_msa_df[position] != '-').sum() ## sum of sequences that contribute to the insertion (not '-')
		if (insertion_count <= threshold):
			to_omit_ls = que_msa_df.index[(que_msa_df[position] != '-')].tolist() ## list of sequences that contribute to insertion
			seqID_omit_set.update(to_omit_ls)
		elif not threshold: ## if threshold is set to zero, omit all offending sequences
			to_omit_ls = que_msa_df.index[(que_msa_df[position] != '-')].tolist() ## list of sequences that contribute to insertion
			seqID_omit_set.update(to_omit_ls)
	##
	
	## find the set of sequences to keep
	
	seqID_start_set = set()
	seqID_start_set.update(list(alias_dict.values()))
	
	seqID_keep_set = seqID_start_set.difference(seqID_omit_set)
	
	## extract sequences, re-name them back from the alias, and strip any insertions
	
	def getkey(value):
		key = list(alias_dict.keys())[list(alias_dict.values()).index(value)]
		return key
	##
	
	output_fasta_ls = []
	
	for seqID in seqID_keep_set:
		for msa_seqID_ind in range(0,len(input_msa_ls),2):
			if (seqID == input_msa_ls[msa_seqID_ind]):
				old_seqID = getkey(seqID)
				seq = input_msa_ls[msa_seqID_ind+1].replace('-','')
				output_fasta_ls.append(old_seqID)
				output_fasta_ls.append(seq)
	##
	
	## print output stats
	
	print('started with ' + str(len(seqID_start_set)) + ' aligned sequences')
	print('identified ' + str(len(seqID_omit_set)) + ' sequences to omit based on insertions at <' + str(min_fraction) + ' fxn positions')
	print('saving ' + str(len(seqID_keep_set)) + ' un-aligned sequences to ' + output_seq_fn)
	
	## save output
	
	def saver(input_name, input_list):
		name_obj = open(input_name, "w")
		for element in input_list:
			if isinstance(element, list):
				element = "\t".join(element)
			name_obj.write(element + "\n")
		##
		name_obj.close()
	##
	
	saver(output_seq_fn, output_fasta_ls)
##

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', type=str, help='input FASTA format MSA to be cleaned e.g.: "iter_1.msa"')
	parser.add_argument('-ref', type=str, help='input newline delim list of reference sequence IDs (with >) to compare to e.g.: "reference.seqIDs"')
	parser.add_argument('-fxn', type=float, default = 0.05, help='min necessary fraction of insertions needed to keep sequence given insertion at a given position (default = 0.05)')
	parser.add_argument('-o', type=str, default = 'cleaned.fa', help='output filename for cleaned FASTA (default = cleaned.fa)')
	args = parser.parse_args()
	main(args.i, args.ref, args.fxn, args.o)
##















