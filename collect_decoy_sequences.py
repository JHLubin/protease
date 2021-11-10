import argparse
import joey_utils as ut
from os.path import basename, join
import pandas as pd
import pyrosetta as pr

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
	parser.add_argument("-name", "--name", type=str,
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


def make_sequences_table_headers(tabulated_residues, ignore_chains=None):
	"""
	Makes headers list for sequences_table
	"""
	headers = ['name', 'total_energy', 'subs_count', 'substitutions']
	for index, row in tabulated_residues.iterrows():
		chain = row['pdb_chain']
		if ignore_chains:
			if chain in ignore_chains:
				continue
		pdb_number = row['pdb_number']
		native = row['residue']
		headers.append('{}_{}{}'.format(chain, native, str(pdb_number)))

	return headers


def main(args):
	# Get reference pose and residue table
	ref_pose = ut.load_pose(args.start_struct, enzdes_cst=args.constraints)
	ref_tabulated_residues = ut.tabulate_pose_residues(ref_pose)
	
	# Define score function
	sf = ut.get_sf(hbnet=args.use_hb_net)

	# Identify list of target PDBs to be analyzed by this cpu
	all_pdbs = ut.collect_pdbs_list(args.decoys_dir)
	pdbs_to_check = ut.partition_list(all_pdbs, *args.parallel_partition)

	# Create collection dataframe
	sequences_table = pd.DataFrame([])

	# Iterate through all PDBs to assemble table of sequence changes
	for pdb in pdbs_to_check:
		# Prepare Pose
		pose = ut.load_pose(pdb, enzdes_cst=args.constraints)

		# Collect dict 
		decoy_props = {}
		decoy_props['subs_list'] = []
		decoy_props['subs_count'] = 0
		decoy_props['name'] = basename(pdb)	
	
		# Any sequence offsets (for loop swaps)
			# Add later

		# total_energy
		decoy_props['total_energy'] = ut.total_energy(pose, sf)

		# constraints_energy
			# Add later
		# interface_energy
			# Add later
		# ddG
			# Add later
		# HBNet
			# Add later

		# All sequence changes from reference 
		# Currently requires poses to be from the same input
		for index, row in ref_tabulated_residues.iterrows():
			chain = row['pdb_chain']
			if args.ignore_chains:
				if chain in args.ignore_chains:
					continue
			site = row['pose_number']
			pdb_number = row['pdb_number']
			native = row['residue']
			candidate = pose.residue(site).name1()

			nres = '{}_{}{}'.format(chain, native, str(pdb_number))
			decoy_props[nres] = candidate
			if candidate != native:
				decoy_props['subs_list'].append(nres + candidate)
				decoy_props['subs_count'] += 1

		decoy_props['substitutions'] = ';'.join(decoy_props['subs_list'])

		# Add dict to main DataFrame
		sequences_table = sequences_table.append(decoy_props, ignore_index=True)
	
	# Set output name
	od = ut.out_directory(args.out_dir)
	out_name = 'decoy_sequences'
	if args.name:
		out_name = args.name
		
	if args.parallel_partition != [1, 1]:
		out_name += '_{}_of_{}'.format(*args.parallel_partition[::-1])
	out_name += '.csv'
	out_name = join(od, out_name)

	# Set headers and output CSV
	headers = make_sequences_table_headers(ref_tabulated_residues, 
		ignore_chains=args.ignore_chains)
	sequences_table = sequences_table[headers]
	sequences_table.to_csv(out_name, index=False)


if __name__ == '__main__':
	args = parse_args()
	pr.init('-run:preserve_header')
	main(args)