import argparse
import joey_utils as ut
import pandas as pd

def parse_args():
	info = """
		Generates a table of representative design decoys, with their sequences
		and interfacial energy changes.
		"""
	parser = argparse.ArgumentParser(description=info)
	parser.add_argument("-w", "--wild_canonical", type=str, required=True,
		help="Specify PDB file with wild-type protein with canonical substrate")
	parser.add_argument("-t", "--wild_target", type=str, required=True,
		help="Specify PDB file with wild-type protein with target substrate")
	parser.add_argument("-c", "--representatives_csv", type=str, required=True,
		help="Specify CSV file listing representative models.")
	parser.add_argument("-d", "--decoys_dir", type=str, required=True,
		help="Specify the directory containing representative design decoys")
	parser.add_argument("-od", "--out_dir",  type=str, 
		help="Name an output directory for the analysis csv. If none \
		specified, uses the same directory as the ")
	parser.add_argument("-pre", "--name_prefix", type=str,
		help="Add a prefix name to the output csv.")
	parser.add_argument("-suf", "--name_suffix", type=str,
		help="Add a suffix name to the output csv.")
	parser.add_argument("-apo", "--make_apo_model", action="store_true", 
		help="Option to generate unbound relaxed structures for per-residue \
		total energy comparison. Significantly increases runtime.")
	parser.add_argument("-cons", "--constraints", type=str, 
		help="If EnzDes constraints are to be used, specify a constraints file")
	parser.add_argument("-hbn", "--use_hb_net", action="store_true", 
		help="Option to include HBnet score term in design.")
	parser.add_argument('-mc', '--main_chain', nargs='*', type=str, default='A',
        help='Specify the main chain(s) of the interacting complex. \
        (Default: A)')
	parser.add_argument('-sc', '--substrate_chain', nargs='*', type=str, 
        default='B', help='Specify the substrate chain(s) of the interacting \
        complex. (Default: B)')
	parser.add_argument("-part", "--parallel_partition", type=int, nargs=2, 
        default=[1,1], help="To parallelize the analysis, enter two integers. \
        The first is the number of partitions, the second is which partition \
        member to run on this processor, from 1 to the number of partitions.")
	args = parser.parse_args()
	return args


def tabulate_per_residue_interactions(pose, scorefunction, cluster, 
	substitutions, main_chain='A', substrate_chain='B'):
	"""
	Given a pose, makes a table that includes the per-residue interfacial 
	energies of each residue. The energy encompases the sum of pairwise 
	energies between a given residue and all residues across the interface.
	"""
	# Copy and score pose
	score_pose = ut.pose_copy(pose)
	scorefunction(score_pose)

	# Make table for energies from pose
	pose_table = ut.tabulate_pose_residues(score_pose)

	# Make a dataframe with interface energies
	## Initialize
	int_e_df = pd.DataFrame([])
	
	## Collect overall column values
	int_e_df['cluster'] = [cluster]
	int_e_df['substitutions'] = [substitutions]
	if substitutions:
		int_e_df['subs_count'] = [len(substitutions.split(';'))]
	else:
		int_e_df['subs_count'] = [0]
	int_e_df['total'] = [ut.total_energy(score_pose, scorefunction)]
	
	## Collect per-residue interfacial energies
	for index, row in pose_table.iterrows():
		# Identify residue number
		target_res = row['pose_number']

		# Make selector for single residue
		target_res_selection = ut.index_selector(target_res)

		# Make selector for interface chain
		if row['pdb_chain'] in main_chain:
			interface_chain = ut.chain_selector(','.join(main_chain))
		elif row['pdb_chain'] in substrate_chain:
			interface_chain = ut.chain_selector(','.join(substrate_chain))
		else:
			print('Main chain:', main_chain)
			print('Substrate chain:', substrate_chain)
			print("Residue's chain", row['pdb_chain'])
			raise ValueError('PDB chain of residue is not among listed chains')

		# Get interaction energy and add it to the table
		int_e = ut.interaction_energy(score_pose, scorefunction, 
			target_res_selection, interface_chain)
		int_e_df[target_res] = [int_e]

	return int_e_df


def main(args):	
	# Read in representatives csv and take partition to be analyzed by this cpu
	full_reps_df = pd.read_csv(args.representatives_csv)
	rows_to_analyze = ut.partition_list(range(len(full_reps_df)), 
		*args.parallel_partition)
	reps_df = full_reps_df.iloc[rows_to_analyze]

	# Define score function
	sf = ut.get_sf(hbnet=args.use_hb_net)

	# Read in reference poses
	wc_pose = ut.load_pose(args.wild_canonical, enzdes_cst=args.constraints)
	wt_pose = ut.load_pose(args.wild_target, enzdes_cst=args.constraints)

	# Tabulate interface energies of reference poses
	wc_table = tabulate_per_residue_interactions(wc_pose, sf, 'wild_canonical', 
		'', args.main_chain, args.substrate_chain)
	wt_table = tabulate_per_residue_interactions(wt_pose, sf, 'wild_target', 
		'', args.main_chain, args.substrate_chain)

	# Initialize main per-residue energy table
	per_res_table = pd.DataFrame([])

	# Add reference pose tables to main table
	per_res_table = per_res_table.append(wc_table, ignore_index=True)
	per_res_table = per_res_table.append(wt_table, ignore_index=True)

	# Iterate through all rows in reps_df and add to main table
	for index, row in reps_df.iterrows():
		# Identify representative PDB
		rep_pdb = row['name']

		# Prepare pose
		pose = ut.load_pose(rep_pdb, path=args.decoys_dir, 
			enzdes_cst=args.constraints)

		# Tabulate interface energies of reference poses
		sf(pose)
		p_table = tabulate_per_residue_interactions(pose, sf, rep_pdb, 
			row['substitutions'], args.main_chain, args.substrate_chain)

		# Add to the main table
		per_res_table = per_res_table.append(p_table, ignore_index=True)
		
	# Set output name
	out_name = ut.output_file_name('interface_energies', path=args.out_dir, 
		extension='csv', prefix=args.name_prefix, suffix=args.name_suffix)
	if args.parallel_partition != [1, 1]:
		out_ext = '_{}_of_{}'.format(*args.parallel_partition[::-1])
		out_name = ut.output_file_name(out_name, suffix=out_ext)

	# Output CSV
	per_res_table.to_csv(out_name, index=False)


if __name__ == '__main__':
	args = parse_args()
	ut.pyrosetta_init(preserve_header=True)
	main(args)