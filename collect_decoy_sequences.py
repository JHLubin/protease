import argparse
import joey_utils as ut
import pandas as pd
from os.path import basename

def parse_args():
	info = """
		Generates a table of all decoys in a specified directory, comparing  
		their sequences to that of a reference. 
		"""
	parser = argparse.ArgumentParser(description=info)
	parser.add_argument("-s", "--start_struct", required=True,
		help="Pick starting PDB")
	parser.add_argument("-d", "--decoys_dir", required=True,
		help="Specify the directory containing the decoys to be analyzed")
	parser.add_argument("-od", "--out_dir", 
		help="Name an output directory for the csv(s)")
	parser.add_argument("-name", "--name", type=str, default='decoy_sequences',
		help="How would you like to name your output(s)? \
		(Default will be a csv called decoy_sequences.)")
	parser.add_argument("-cons", "--constraints", type=str, 
		help="If EnzDes constraints are to be used, specify a constraints file")
	parser.add_argument("-hbn", "--use_hb_net", action="store_true", 
		help="Option to include HBnet score term in design.")
	parser.add_argument('-ic', '--ignore_chains', type=str, nargs='*', \
        help='Specify chain(s) to exclude from comparison.')
	parser.add_argument("-part", "--parallel_partition", type=int, nargs=2, 
        default=[1,1], help="To parallelize the analysis, enter two integers. \
        The first is the number of partitions, the second is which partition \
        member to run on this processor, from 1 to the number of partitions.")
	args = parser.parse_args()
	return args


def make_sequences_table_headers(atom_table, ignore_chains=None):
	"""
	Makes headers list for sequences_table
	"""
	headers = ['name', 'total_energy', 'total_cst', 'subs_count', 
		'substitutions']
	for index, row in atom_table.iterrows():
		chain = row['chain_id']
		if ignore_chains:
			if chain in ignore_chains:
				continue
		pdb_number = row['res_number']
		native = ut.convert_aa_name(row['res_name'])
		headers.append('{}_{}{}'.format(chain, native, str(pdb_number)))

	return headers


def main(args):
	# Get reference pose and residue and scores tables
	ref_pdb = ut.find_maybe_compressed_file(args.start_struct, gunzip=True)
	ref_tabulated_atoms = ut.tabulate_pdb_atom_lines(ref_pdb)
	ref_ca_atoms = ref_tabulated_atoms[ref_tabulated_atoms['atom_name']=='CA']
	
	# Identify list of target PDBs to be analyzed by this cpu
	all_pdbs = ut.collect_pdbs_list(args.decoys_dir)
	pdbs_to_check = ut.partition_list(all_pdbs, *args.parallel_partition)

	# Create collection dataframe and headers list
	sequences_table = pd.DataFrame([])
	headers = make_sequences_table_headers(ref_ca_atoms, 
		ignore_chains=args.ignore_chains)

	# Iterate through all PDBs to assemble table of sequence changes
	for p in pdbs_to_check:
		print('checking pdb {} of {}: {}'.format(
			pdbs_to_check.index(p) + 1, len(pdbs_to_check), p))

		# Get pdb residue and scores tables
		pdb = ut.find_maybe_compressed_file(p, gunzip=True)
		tabulated_atoms = ut.tabulate_pdb_atom_lines(pdb)
		tabulated_scores = ut.tabulate_pdb_energy_lines(pdb)
		score_totals = tabulated_scores[tabulated_scores['pdb_number']=='pose']

		# Collect dict 
		decoy_props = {}
		decoy_props['subs_list'] = []
		decoy_props['subs_count'] = 0
		decoy_props['name'] = basename(pdb)	
	
		# Any sequence offsets (for loop swaps)
			# Add later

		# total_energy
		decoy_props['total_energy'] = float(score_totals['total'])

		# constraints_energy
		cst_cols = score_totals.columns[
			score_totals.columns.str.contains('constraint')]
		decoy_props['total_cst'] = score_totals[cst_cols].to_numpy().sum()

		# interface_energy
			# Add later
		# ddG
			# Add later
		# HBNet
			# Add later

		# Iterate through reference residues and collect sequence changes 
		# Currently requires poses to be from the same input
		for index, row in ref_ca_atoms.iterrows():
			chain = row['chain_id']
			if args.ignore_chains:
				if chain in args.ignore_chains:
					continue
			#site = row['pose_number']
			pdb_number = row['res_number']
			native = ut.convert_aa_name(row['res_name'])
			candidate = ut.find_aa_in_pdb_atom_table(
				tabulated_atoms, chain, pdb_number)

			nres = '{}_{}{}'.format(chain, native, str(pdb_number))
			decoy_props[nres] = candidate
			if candidate != native:
				decoy_props['subs_list'].append(nres + candidate)
				decoy_props['subs_count'] += 1

		decoy_props['substitutions'] = ';'.join(decoy_props['subs_list'])

		# Add dict to main DataFrame
		sequences_table = sequences_table.append(decoy_props, ignore_index=True)
	
	# Set output name		
	suf = None
	if args.parallel_partition != [1, 1]:
		suf = '{}_of_{}'.format(*args.parallel_partition[::-1])
	out_name = ut.output_file_name(args.name, path=args.out_dir, 
		extension='csv', suffix=suf)

	# Trim/reorder columns and output CSV
	sequences_table = sequences_table[headers]
	sequences_table.to_csv(out_name, index=False)


if __name__ == '__main__':
	args = parse_args()
	main(args)