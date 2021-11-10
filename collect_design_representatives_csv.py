import argparse
import joey_utils as ut
from os.path import join
import pandas as pd
from shutil import copyfile

def parse_args():
	info = """
		Takes one or multiple CSV files output by collect_decoy_sequences and 
		collapses them to a single representative of each unique sequence, 
		taking the lowest energy model of the set as the representative. The 
		output is a reduced table, which can be further reduced to exclude 
		residue columns in which no decoys have substitutions. Also, as an 
		option, the representative decoys can be copied into another directory.
		"""
	parser = argparse.ArgumentParser(description=info)
	parser.add_argument("-c", "--csvs", required=True, type=str, nargs='*', 
		help="Identify the input csv(s) from collect_decoy_sequences to merge \
		and identify unique design sequences.")
	parser.add_argument("-od", "--out_dir", type=str, 
		help="Name an output directory for the consolidated csv and decoys")
	parser.add_argument("-e", "--copy_extract", action='store_true',
		help="Copy the representative models to a folder specified by \
		the decoys_dir argument.")
	parser.add_argument("-d", "--decoys_dir", type=str, 
		help="Specify the directory containing the decoys to be copied. \
		(Unnecessary if copy_extract argument is not given, required if it is.")
	parser.add_argument("-name", "--name", type=str,
		help="Add a prefix name to the output csv. (By default, it is \
		design_models.)")
	parser.add_argument("-trim", "--trim_unchanged_cols", action='store_true',
		help="Exclude columns with no sequence changes from final output.")
	args = parser.parse_args()
	if args.copy_extract and not args.decoys_dir:
	    parser.error('If using copy_extract, decoys_dir is required.')
	return args


def make_models_table(csvs, trim=False):
	"""
	Reads csv(s) from collect_decoy_sequences and produces a reduced version
	with only one representative for each unique sequence (the one with the 
	lowest total score) and optionally trimming out all sites with no 
	substitutions. Adds a column with the count of decoys that had that unique
	sequence. 
	"""
	# Aggregate input csvs into models table
	sequences_table = pd.DataFrame([])
	for csv in csvs:
		df = pd.read_csv(csv, engine='python')
		sequences_table = sequences_table.append(df)

	# Collapse models table to single representative of each sequence
	rep_set = sequences_table.loc[
		sequences_table.groupby('substitutions')['total_energy'].idxmin()]

	# Add count columns to models table
	rep_set['count'] = rep_set.apply(lambda row: 
		len(sequences_table[sequences_table['substitutions'] == 
		row['substitutions']]), axis='columns')

	# Reorder columns so count and substitutions come before residues
	cmain = ['name', 'total_energy', 'subs_count', 'substitutions', 'count']
	csite = [i for i in list(rep_set) if i not in cmain]
	if trim: # Eliminate unchanged columns from models table
		cols = cmain + [i for i in csite if len(rep_set[i].unique()) > 1]
	else:
		cols = cmain + csite	
	rep_set = rep_set[cols]

	return rep_set


def main(args):
	# Make models table
	models_df = make_models_table(args.csvs, trim=args.trim_unchanged_cols)

	# Output models CSV file
	out_name = 'designs_models.csv'
	if args.name:
		out_name = '_'.join([args.name, out_name])
	out_name = join(ut.out_directory(args.out_dir), out_name)
	models_df.to_csv(out_name, index=False)

	# Copy out representative models
	if args.copy_extract:
		for decoy in models_df['name'].unique():
			d_in = join(args.decoys_dir, decoy)
			d_out = join(args.out_dir, decoy)
			copyfile(d_in, d_out)


if __name__ == '__main__':
	args = parse_args()
	main(args)