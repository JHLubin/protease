#!/usr/bin/python
""" 
PyRosetta4, Python 3.5
Joseph Lubin, 2019

Pipeline for rapidly modeling protease-substrate combinations, and using 
FastRelax and FastDesign to explore potentially better interacting variants

It is assumed in the program that the input PDB structure will have a set of
enzdes constraint comments at the beginning of the document, and that the 
protease is chain A and the substrate is chain B.

Sample commands:
python design_protease.py -s protease_design/start_proteases/HCV.pdb 
-od test -name despep -seq DVDAR -site 198 -ps "198-202" -cons 
protease_design/ly104.cst -cr 72 96 154 -dprot 0 -dpep 1 -n 100 -mm 138 I 
-mm 170 Q -mm 171 S -mm 173 I -mm 175 K -mm 183 R

python design_protease.py  -name htra1_protease_asyn92 
-s fibrils_collaboration/htra1_protease_p6-2p.pdb -seq SIAAATGF 
-od fibrils_collaboration/relax_protease_on_92 -site 212 -cr 61 91 169 
-cons fibrils_collaboration/htra1_protease.cst -n 10 -dprot 0 -dpep 0
"""
from __future__ import print_function # For compatability with Python 2.7
import argparse
from os import makedirs
from os.path import basename, isdir, isfile, join
from pyrosetta import *
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import \
	IncludeCurrent, ExtraRotamers, OperateOnResidueSubset, \
	PreventRepackingRLT, RestrictAbsentCanonicalAASRLT, RestrictToRepackingRLT
from pyrosetta.rosetta.core.scoring import ScoreType
from pyrosetta.rosetta.core.select.residue_selector import \
	AndResidueSelector, ChainSelector, InterGroupInterfaceByVectorSelector,\
	NeighborhoodResidueSelector, NotResidueSelector, OrResidueSelector, \
	ResidueIndexSelector
from pyrosetta.rosetta.core.simple_metrics.metrics import \
	SelectedResiduesMetric
from pyrosetta.rosetta.protocols.constraint_generator import \
	AddConstraints, CoordinateConstraintGenerator, \
	HydrogenBondConstraintGenerator
from pyrosetta.rosetta.protocols.enzdes import ADD_NEW, AddOrRemoveMatchCsts
from pyrosetta.rosetta.protocols.minimization_packing import \
	MinMover, PackRotamersMover
from pyrosetta.rosetta.protocols.protein_interface_design import \
	FavorNativeResidue
from pyrosetta.rosetta.protocols.relax import FastRelax
from random import randint
from sys import exit

def parse_args():
	info = "Design a protease around a peptide sequence"
	parser = argparse.ArgumentParser(description=info)
	parser.add_argument("-s", "--start_struct", required=True,
		default='start_proteases/HCV.pdb', help="Pick starting PDB")
	parser.add_argument("-od", "--out_dir", required=True,
		help="Name an output directory for decoys")
	parser.add_argument("-name", "--name", type=str,
		help="How would you like to name your outputs? \
		(Default will use the name of the input PDB file.)")
	parser.add_argument("-seq", "--sequence", required=True, type=str,
		help="What substrate sequence do you want to thread")
	parser.add_argument("-site", "--subst_site", required=True, type=int,
		help="Where in the pose should the substrate sequence begin \
		threading? (A 5-letter substitution ending with 202 should start at \
		198).")
	parser.add_argument("-cr", "--cat_res", type=int, nargs='+', 
		default=None, help="The catalytic residues of the protease, \
		excluded from design. (By default, no residues are so designated.)")
	parser.add_argument("-ps", "--pep_subset", type=str, default=None, 
		help='Select the subset of the peptide around which to design, as a \
		string of "first_res-last_res". (Ex: "198-202") Otherwise, design \
		will be performed around the full peptide. These numbers should be in \
		pose numbers, not PDB numbers, if the two differ.')
	parser.add_argument("-cons", "--constraints", type=str, 
		default='ly104.cst', help="Pick constraints file")
	parser.add_argument("-dprot", "--design_protease", type=str2bool, default=1, 
		help="Allow design on the protease near the peptide? 0 \
		for False, 1 for True. (Default: True)")
	parser.add_argument("-dpep", "--design_peptide", type=str2bool, default=0,
		help="Allow design on the peptide? 0 for False, 1 for \
		True. (Default: False)")
	parser.add_argument("-hbn", "--use_hb_net", action="store_true", 
		help="Option to include HBnet score term in design.")
	parser.add_argument("-n", "--number_decoys", type=int, default=10, 
		help="How many decoys should be made? (Default is 10.)")
	parser.add_argument("-mm", "--mutations", nargs=2, action='append',
		help="Manually input mutations in the format [site] [one-letter res]. \
		Accepts multiple uses. (Ex: -mm 138 I -mm 183 R) Note, if you intend \
		to change the catalytic residues, you must edit the PDB's enzdes \
		comments as well, or applying constraints won't work properly. Uses \
		pose numbering, which may differ from PDB.")
	parser.add_argument("-cp", "--constrain_peptide", action="store_true",
		help="Option to add coordinate constraints to the substrate peptide \
		backbone atoms. False by default.")
	parser.add_argument("-rtc", "--res_penalty", type=float, default=0.5,
		help="Set the constraint penalty for favoring native residue type in \
		design. (Default is 0.5.)")
	parser.add_argument("-hbc", "--hbond_constraints", nargs=2, type=str, 
		help="Apply H-bond constraints to preserve beta sheets. Requires two \
		string inputs in the form used by ResidueIndexSelector for lists of \
		residues to bond. (Use pose numbers, not PDB numbers.) '215,217' and \
		'185,187' for Htra1 with the extended peptide, '199-203' and '170-174' \
		for HCV.")
	parser.add_argument("-init", "--extra_init_options", type=str, 
		action='append', help='Extra init options for Rosetta (ex: \
		"-extra_res_fa LG1.params")')
	parser.add_argument("-nrc", "--no_relax_comparison", action="store_true",
		help="Prevent generation of a relaxed decoy along with the designed \
		model. By default, one will be produced, but this option will prevent \
		that from happening. Doing so may save time if many trajectories are \
		being run, and a smaller set of relaxed models is generated for \
		comparison. This can be done by executing the same command, but with \
		-dprot 0 and -dpep 0")
	parser.add_argument("-test", "--test_mode", action="store_true", 
		help="For debugging: test protocol, exiting before generating decoys.")
	parser.add_argument("-v", "--verbose", action="store_true", 
		help="For debugging: Don't mute Rosetta output.")
	args = parser.parse_args()

	if args.verbose:
		print(args)

	return args


def str2bool(v):
	""" Converts a number of potential string inputs to boolean """
	if isinstance(v, bool):
		return v
	if v.lower() in ('yes', 'true', 't', 'y', '1'):
		return True
	elif v.lower() in ('no', 'false', 'f', 'n', '0'):
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected.')


def init_opts(extra_opts, cst_file='ly104.cst', verbose=False):
	""" Produces a list of init options for PyRosetta, including cst file """
	ros_opts = '-ex1 -ex2  -use_input_sc -flip_HNQ'
	ros_opts += ' -cst_fa_weight 1.0 -run:preserve_header'
	ros_opts += ' -write_pdb_link_records false'

	if not verbose:
		ros_opts += ' -mute all'

	if extra_opts:
		for i in extra_opts:
			ros_opts += ' {}'.format(i)
	print(ros_opts)
	return ros_opts


def readfile(file_name):
	""" Opens a file in read-mode and returns a list of the text lines """
	with open(file_name, 'r') as r:
		lines = r.readlines()

	return lines

######### Threading ##########################################################

def make_residue_changes(pose, sf, subst_seq, subst_start, cat_res, manual_muts):
	"""
	Applies substrate sequence changes and manual mutations to a given pose.
	This is done through repacking, so unlike SimpleThreadingMover, the side 
	chains don't begin clashing. This means that the residue selectors will be 
	more accurate, and that design can begin without an initial relax step.

	pose is a Rosetta pose
	sf is a Rosetta scorefunction
	subst_seq is a string (doesn't need to be uppercase)
	subst_start is an integer corresponding to the first of a contiguous block 
		of  residues to re-sequence
	manual_muts is a list of two-member lists, of the following form: 
		[site, single-letter residue name]
	"""
	# Create dict of {res: AA} for changes to make
	res_changes = {}

	# Add substrate threading to list of res changes
	print("\nInserting substrate sequence:\n{}".format(subst_seq))
	subst_range = range(subst_start, subst_start + len(subst_seq))
	for n, i in enumerate(subst_range):
		res_changes[i] = subst_seq[n].upper()

	# Add manual mutations list
	if manual_muts:
		print("\nApplying point substitutions:")
		for m in manual_muts:
			res_changes[int(m[0])] = m[1].upper()
			print(m[0], m[1].upper())
			if cat_res:
				if int(m[0]) in cat_res:
					print('Warning: {} is a catalytic residue'.format(m[0]))

	# Make TaskFactory to input changes
	mobile_residues = OrResidueSelector() # Keep list of mobile residues
	tf = TaskFactory()

	# Force packing to target residue for each desired change
	for r, aa in res_changes.items():
		res_selection = ResidueIndexSelector(str(r))
		restriction = RestrictAbsentCanonicalAASRLT()
		restriction.aas_to_keep(aa.upper())
		tf.push_back(OperateOnResidueSubset(restriction,res_selection))
		mobile_residues.add_residue_selector(res_selection)

	# Repack nearby residues to accommodate substitutions
	shell = NeighborhoodResidueSelector()
	shell.set_focus_selector(mobile_residues)
	shell.set_include_focus_in_subset(False)
	shell.set_distance(8)
	
	# Exclude catalytic residues from repackable shell if not mutated
	if cat_res:
		skip_res = [i for i in cat_res if i not in res_changes.keys()]
		catalytic = ResidueIndexSelector(','.join([str(i) for i in skip_res]))
		not_catalytic = NotResidueSelector(catalytic)
		shell = selector_intersection(shell, not_catalytic)
	
	restrict = RestrictToRepackingRLT()
	tf.push_back(OperateOnResidueSubset(restrict, shell))
	
	# Prevent repacking of all other residues
	unchanging = NotResidueSelector(OrResidueSelector(mobile_residues, shell))
	prevent = PreventRepackingRLT()
	tf.push_back(OperateOnResidueSubset(prevent, unchanging))

	# Apply changes with PackRotamersMover
	pt = tf.create_task_and_apply_taskoperations(pose)
	prm = PackRotamersMover(sf, pt)
	mutated_pose = Pose(pose)
	prm.apply(mutated_pose)

	# Set up fixed-backbone movemap for minimization
	movemap = MoveMap()
	movemap.set_bb(False)
	movemap.set_chi(True)
	movemap.set_jump(False)

	# Minimize
	min_mover = MinMover()
	min_mover.movemap(movemap)
	min_mover.score_function(sf)

	# Apply the MinMover to the modified Pose
	min_mover.apply(mutated_pose)

	return mutated_pose


def random_aa(length):
	""" 
	Returns a string of random 1-letter amino acid names from the cannonical 
	20, to a specified length.
	"""

	aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
				'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
	
	aa_string = ""

	for aa in range(length):
		rand_index = randint(0,19)
		aa_string += aa_list[rand_index]

	return aa_string

######### Residue selection ##################################################

def selector_intersection(*selectors):
	""" Returns the intersection of any set of selectors """
	intersect_selection = AndResidueSelector()
	for s in selectors:
		intersect_selection.add_residue_selector(s)

	return intersect_selection


def selector_union(*selectors):
	""" Returns the intersection of any set of selectors """
	union_selection = OrResidueSelector()
	for s in selectors:
		union_selection.add_residue_selector(s)

	return union_selection


def mutable_residues_selector(protease_selection, peptide_selection,
	catalytic_selection, design_peptide=False):
	"""
	Selects the residues in a shell around the peptide using the 
	InterGroupInterfaceByVectorSelector residue selector
	"""
	# Making protease shell selector (includes peptide)
	first_shell_select = InterGroupInterfaceByVectorSelector()
	first_shell_select.group1_selector(protease_selection) 
	first_shell_select.group2_selector(peptide_selection) 

	# Excluding the catalytic residues, peptide (if not designed)
	not_cats_sel = NotResidueSelector(catalytic_selection)
	if design_peptide:
		mutable_selection = selector_intersection(
			not_cats_sel, first_shell_select)
	else:
		mutable_selection = selector_intersection(
			not_cats_sel, first_shell_select, protease_selection)

	return mutable_selection


def packable_residues_selector(
	peptide_selection, mutable_selection, catalytic_selection):
	"""
	Selects the shell of neighbor residues to repack. Packable set should not
	include the mutable set, since the action is RestrictToRepacking.
	"""
	# Making negative selections for mutable and catalytic
	not_mutable = NotResidueSelector(mutable_selection)

	# Selecting residues near mutable shell
	near_mutable = InterGroupInterfaceByVectorSelector()
	near_mutable.group1_selector(not_mutable)
	near_mutable.group2_selector(mutable_selection)

	# Selecting residues near the peptide, with wider range for BB mobility
	near_pep = InterGroupInterfaceByVectorSelector()
	near_pep.group1_selector(not_mutable)
	near_pep.group2_selector(peptide_selection)

	# Combining selections for peptide and near peptide and near mutable
	inclusive_packable = selector_union(
		catalytic_selection, near_mutable, near_pep, peptide_selection)
	
	# Setting up exclusion of mutable residues
	exclusive_packable = selector_intersection(inclusive_packable, not_mutable)

	return exclusive_packable


def select_residues( 
	cat_res, peptide_subset, design_protease=True, design_peptide=False):
	""" 
	Makes residue selectors for protease sections. Requires manual input for 
	which residues are catalytic and whether only part of the peptide should be 
	selected. Options for whether the peptide is designable (false by default) 
	and whether the protease is designable (true by default). Assumes that the 
	protease is chain A and the peptide is chain B.
	"""
	residue_selectors = {}

	# Protease residues. Protease assumed to be chain A
	protease = ChainSelector("A")
	residue_selectors['protease'] = protease
	
	# Peptide residues. Peptide assumed to be chain B, unless range specified
	if peptide_subset:
		peptide = ResidueIndexSelector(peptide_subset)
	else:
		peptide = ChainSelector("B")
	residue_selectors['peptide'] = peptide

	# Catalytic residues. ResidueIndexSelector needs a string, not a list.
	if cat_res:
		cats_as_str = ','.join([str(i) for i in cat_res]) 
		catalytic = ResidueIndexSelector(cats_as_str)
	else:
		# If no catalytic residues are given, return a null selector
		catalytic = selector_intersection(protease, peptide) # Empty set
	residue_selectors['catalytic'] = catalytic
			
	# Designable residues. May include protease and peptide, just one, or none
	if design_protease:
		mutable = mutable_residues_selector(protease, peptide,
			catalytic, design_peptide)
	elif design_peptide:
		mutable = peptide
	else: # Neither protease not peptide designable
		mutable = selector_intersection(protease, peptide) # Empty set
	residue_selectors['mutable'] = mutable

	# Packable residues. Centered around the peptide and designable set
	packable = packable_residues_selector(peptide, mutable, catalytic)
	residue_selectors['packable'] = packable
	
	# Immobile residues. Catalytic residues and everything that isn't mutable 
	# or packable
	immobile = NotResidueSelector(selector_union(mutable, packable))
	residue_selectors['immobile'] = immobile

	return residue_selectors


def selector_to_list(pose, selector):
	""" Converts a selector output vector to a list of selected residues """
	# Set up SelectedResiduesMetric
	srm = SelectedResiduesMetric()
	srm.set_residue_selector(selector)
	srm.set_output_in_rosetta_num(True)
	
	# Collect selection, and convert to a list
	sel_res_str = srm.calculate(pose)
	sel_res_str_list = sel_res_str.split(',')
	if sel_res_str_list == ['']: # Avoid errors when list is empty
		sel_res_str_list = []
	selection_list = [int(i) for i in sel_res_str_list]

	return selection_list 

######### Setup ##############################################################

def apply_constraints(pose, cst_file='protease_design/ly104.cst'):
	""" Applies the constraints form the input CST file to a pose """
	cstm = AddOrRemoveMatchCsts()
	cstm.set_cst_action(ADD_NEW)
	cstm.cstfile(cst_file)
	cstm.apply(pose)
	return pose


def coord_constrain_peptide(pose, selection=ChainSelector('B')):
	""" Applies backbone coordinate constraints to a selection of a pose """
	cg = CoordinateConstraintGenerator()
	if selection:
		cg.set_residue_selector(selection)
	ac = AddConstraints()
	ac.add_generator(cg)
	ac.apply(pose)
	return pose


def res_type_constrain(pose, penalty=0.5):
	""" Apply residue type constraints to a pose """
	FavorNativeResidue(pose, penalty)
	return 


def apply_hbond_constraints(pose, first_strand, second_strand):
	"""
	Adds H-bond constraints to a pose, given lists of residues in two beta 
	strands. Lists should be given as strings of residues, such as '215,217' 
	and '185,187' for Htra1 with the extended peptide, or '199-203' and 
	'170-174' for HCV. 
	"""
	# Create selectors for residues in each set
	hb_res_1 = ResidueIndexSelector(first_strand)
	hb_res_2 = ResidueIndexSelector(second_strand)

	# Set up constraints generation
	hbcg = HydrogenBondConstraintGenerator()
	hbcg.set_residue_selector1(hb_res_1)
	hbcg.set_residue_selector2(hb_res_2)

	# Apply, return pose
	hbcg.apply(pose)
	return pose


def make_move_map(pose, selectors):
	""" 
	Makes a movemap for a protease-peptide system, with all non-peptide 
	residue backbones fixed, and side chains mobile for all packable and 
	designable residues. 

	Takes a dict of selectors generated by select_residues.
	"""
	mm = MoveMap()

	# Mobile backbone for peptide
	for i in selector_to_list(pose, selectors['peptide']):
		mm.set_bb(i, True)
	
	# Mobile side chains for all mutable and packable residues
	for i in selector_to_list(pose, selectors['mutable']):
		mm.set_chi(i, True)
	for i in selector_to_list(pose, selectors['packable']):
		mm.set_chi(i, True)

	return mm


def make_task_factory(pose, residue_selectors):
	""" 
	Makes a TaskFactory with operations that leave the mutable residues 
	designable, restricts the nearby residues to repacking, and prevents 
	repacking of other residues.
	"""
	mutable_set = residue_selectors['mutable']
	repack_set = residue_selectors['packable']
	immobile_set = residue_selectors['immobile']

	prevent = PreventRepackingRLT() # No repack, no design
	repack = RestrictToRepackingRLT() # No design

	tf = TaskFactory()

	# Expand rotamer library
	tf.push_back(IncludeCurrent())
	tf.push_back(ExtraRotamers(0, 1, 1))
	tf.push_back(ExtraRotamers(0, 2, 1))

	# Prevent repacking of immobile and restrict to repacking of not designable
	tf.push_back(OperateOnResidueSubset(prevent, immobile_set))
	tf.push_back(OperateOnResidueSubset(repack, repack_set))

	# Prevent cysteines in design
	restriction = RestrictAbsentCanonicalAASRLT()
	restriction.aas_to_keep('ADEFGHIKLMNPQRSTVWY')
	for res in selector_to_list(pose, mutable_set):
		if pose.residue(res).name1() != 'C':
			res_sel = ResidueIndexSelector(str(res))
			tf.push_back(OperateOnResidueSubset(restriction, res_sel))

	return tf


def get_score_function(ref15=True, constraints=True, hbnet=False):
	""" Returns either default or weighted REF2015 with or without hbnet """
	# If including default REF2015, start from there, otherwise start from null
	if ref15:
		sf = get_fa_scorefxn()
	else:
		sf = ScoreFunction()

	# Picking between constraints and not
	if constraints:
		sf.set_weight(ScoreType.atom_pair_constraint, 1)
		sf.set_weight(ScoreType.coordinate_constraint, 1)
		sf.set_weight(ScoreType.angle_constraint, 1)
		sf.set_weight(ScoreType.dihedral_constraint, 1)
		sf.set_weight(ScoreType.res_type_constraint, 1)

	# Optionally adding in hbnet
	if hbnet:
		sf.set_weight(ScoreType.hbnet, 1)

	return sf

######### Design Protocols ###################################################

def fastrelax(pose, score_function, movemap, taskfactory=None):
	""" 
	Runs the FastRelax protocol on a pose, using given score function and 
	movemap, and optionally a task factory. By default, FastRelax will not do 
	design. However, given a task factory that enables design, it functions 
	like FastDesign.
	"""
	relax = FastRelax()
	relax.set_scorefxn(score_function)
	relax.set_movemap(movemap)
	if taskfactory:
		relax.set_task_factory(taskfactory)

	pp = Pose(pose)
	relax.apply(pp)

	return pp


def jd_design(name, decoy_count, pose, score_function, movemap, task_factory, 
	save_wt=True):
	""" Runs job distributor with relax and design protocols """
	print('\n')

	current_decoy = 1
	jd = PyJobDistributor(name, decoy_count, score_function)
	while not jd.job_complete:
		pp = Pose(pose)

		# Relaxing to provide a WT comparison to design
		if save_wt:
			print('Relaxing...')
			relax_name = jd.current_name.replace('designed', 'relaxed')
			relaxed_pose = fastrelax(pp, score_function, movemap)
			relaxed_pose.dump_pdb(relax_name)

		# Doing design and outputting decoy
		print('Designing...')
		pp = fastrelax(pp, score_function, movemap, taskfactory=task_factory)

		print('Decoy {}/{} complete\n'.format(current_decoy, decoy_count))
		current_decoy += 1
		jd.output_decoy(pp)

	return

######### Main ###############################################################

def write_opts(name, args):
	""" Write an output file listing the options """
	# Making output name
	opts_out_name = name + '_options.txt'

	# Make list of arguments
	#arg_list = []
	#template = '{:30}{}'
	#for arg in args:
	#	arg_list.append(template.format(arg, args[arg]))

	# Writing file
	with open(opts_out_name, 'w') as w:
		#w.writelines(arg_list)
		w.write(str(args))


def test_and_exit(args, opts, residue_selectors, pose, name):
	""" Prints info then exits """
	print('\n\nArgs:')
	print(args)
	print('\nRosetta options:')
	print(opts)
	print('\nSelectors:')
	for k, v in residue_selectors.items():
		print('\t',k)
		print('\t',selector_to_list(pose,v))
	print('\nSequence')
	print(pose.sequence())
	print('\nName:')
	print(name)
	print('\n')
	pose.dump_pdb(name.replace('designed', 'test.pdb'))

	exit()


def main(args):
	# Initializing PyRosetta
	ros_opts = init_opts(args.extra_init_options, cst_file=args.constraints,
		verbose=args.verbose)
	init(options=ros_opts)

	# Destination folder for PDB files
	dir_name = args.out_dir
	if not args.test_mode:
		if not isdir(dir_name):
			print('\nMaking directory: {}'.format(dir_name))
			makedirs(dir_name)

	# Getting name for outputs
	if args.name:
		out_name = args.name
	else:
		out_name = basename(args.start_struct)
		# strip out .pdb or .pdb.gz extension
		out_name = out_name.replace('.pdb', '').replace('.gz', '')
	dec_name = join(dir_name, out_name)

	# Writing inputs list file
	if not args.test_mode:
		write_opts(join(dir_name, out_name), args)

	# Getting score function
	sf = get_score_function(constraints=True, hbnet=args.use_hb_net)	

	# Preparing pose, with substrate threading, manual mutations, constraints
	pose = pose_from_pdb(args.start_struct)

	if args.constrain_peptide:
		pose = coord_constrain_peptide(pose)
	if args.hbond_constraints:
		pose = apply_hbond_constraints(pose, *args.hbond_constraints)
	pose = make_residue_changes(pose, sf, args.sequence, 
		args.subst_site, args.cat_res, args.mutations)
	res_type_constrain(pose, penalty=args.res_penalty)
	pose = apply_constraints(pose, args.constraints)

	# Making residue selectors
	residue_selectors = select_residues(args.cat_res, args.pep_subset, 
		design_protease=args.design_protease, 
		design_peptide=args.design_peptide)

	# Creating movemap, and taskfactory for design
	mm = make_move_map(pose, residue_selectors)
	tf = make_task_factory(pose, residue_selectors)

	# Running relax and design protocol
	if args.design_protease or args.design_peptide:
		dec_name += '_designed'
	else: 
		dec_name += '_relaxed'

	save_wt = False
	if args.design_protease or args.design_peptide:
		if not args.no_relax_comparison:
			save_wt = True

	if args.test_mode:
		test_and_exit(args, ros_opts, residue_selectors, pose, dec_name)

	#with open(dec_name + '_inputs.txt', 'w') as w:
	#	for arg in vars(args):
	#		w.write('{}\t{}\n'.format([arg, getattr(args, arg)]))

	jd_design(dec_name, args.number_decoys, pose, sf, mm, tf, save_wt=save_wt)
	#relaxed_pose = fastrelax(pose, sf, mm)
	#io.poses_to_silent(relaxed_pose, 'test/trial')


if __name__ == '__main__':
	args = parse_args()
	main(args)