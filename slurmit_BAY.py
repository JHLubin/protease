#!/usr/bin/python

import os, sys, argparse, subprocess, glob

#This is an exception class to catch Argument errors.
class ArgError(Exception):
	def __init__(self, msg):
		self.msg = msg
		print(self.msg)
		
#This function returns a list of tuples.  Each item in the list corresponds to one argument taken from the command line.
#The members of the tuple describe how the variable should be handled in various places.
#To change the behaviour of an argument or to add new variables, add or modify the tuples here.
def gen_varkey():
	vars = []
	
	#(0 arg.varname, 1 SLURM variable name, 2 variable type, 3 default, 4 required, 5 class, 6 help)
	#Classes: -regular: A SLURM variable that takes an argument.
	#		  -boolean: A SLURM variable that doesn't take an argument.  Decides whether to include or exclude the line in the sbatch script.
	#		  -array: A SLURM variable related to arrays.
	#		  -other: Something that isn't a SLURM variable.  This includes non-SBATCH commands and script execution commands.
	vars.append(('job', 'job-name', str, None, True, 'regular', 'The SLURM name for your job (--job-name).  Will be used to generate log/err filenames if not given, in which case it will genereate the files in the current directory.'))
	vars.append(('partition', 'partition', str, 'p_sdk94_1,main', False, 'regular', 'The SLURM partition to run on (--partition). Default: p_sdk94_1,p_sdk94_0,main.'))
	vars.append(('requeue', 'requeue', bool, True, False, 'boolean', 'Run SLURM with the requeue option turned on (--requeue).'))
	vars.append(('openmode', 'open-mode', str, None, False, 'regular', 'Should the log and err files be overwritten (truncate) or appended (append) if the file already exists.  The default is truncate, unless requeue is set, in which case the default is append.'))
	vars.append(('tasks', 'ntasks', int, 140, False, 'regular', 'The number of SLURM tasks to request (--ntasks). Default 140.'))
	vars.append(('cpus', 'cpus-per-task', int, 1, False, 'regular', 'The number of CPUs to request for each task (--cpus-per-task). Default 1.'))
	vars.append(('usearray', None, bool, False, False, 'other', "Is this job being submitted as an array?  If so, either the --array flag must be given with the number of jobs, or --arraygen must be used to generate the array members.  The former case is for running a number of identical jobs (or jobs that differ only by a sequential number in an input filename), while the latter is for running an array with a list of arbitrary filenames."))
	vars.append(('array', 'array', str, None, False, 'array', "Enter the number of array elements to include.  For example, an array with 10 subjobs numbered 1-10 should be given --array 1-10.  If you include this without --arraygen (i.e. array of numbered elements), you must put 'percenta' in --outfiles/--log/--err and 'dollarjob' in --command.  If --arraygen is included, --array can be autogenerated, or alternatively a subset of indices to run can be specified (useful to repeat specific array indices).  In the latter case, the indices start at 0."))
	vars.append(('arraygen', None, str, None, False, 'other', "Specify what files to loop over in the array.  For example, set to '*.pdb' or 'mydir/*.pdb'.  Surround with single quotes to avoid shell expansion.  Use 'dollarfile' in your --command parameter to reference the current file from arraygen."))
	vars.append(('arrayformat', None, str, '${arr[$job]}', False, 'other', "Specify bash code to format the array strings.  For example, to remove the pdb extension, set to ${arr[$job]%%.pdb} .  If left out, the default is to use the array elements unmodified.  You must use arr and job as the array and array element variable names."))
	vars.append(('mem', 'mem', str, "16000", False, 'regular', 'The memory to reserve (--mem). Default 16000'))
	vars.append(('outfiles', None, str, "", False, 'other', 'The name of the log and error files (--output and --err).  This will use the same name for both files, with the extensions .log and .err.'))
	vars.append(('log', 'output', str, "", False, 'regular', 'The name of the log file (--output).  Do not use this together with the outfiles option.'))
	vars.append(('err', 'error', str, "", False, 'regular', 'The name of the error file (--err).  Do not use this together with the outfiles option.'))
	vars.append(('time', 'time', str, "3-00:00:00", False, 'regular', 'The maximum walltime allowed for the job (--time).'))
	vars.append(('begin', 'begin', str, "now", False, 'regular', 'The time to start the script.  By default, starts immediately. May overload the job scheduler.'))
	vars.append(('execute', None, bool, True, False, 'other', 'Do you want to execute the script, or just make the script file?'))
	vars.append(('cleanup', None, bool, True, False, 'other', 'Do you want to delete the script after it is submitted?'))
	vars.append(('command', None, str, None, True, 'other', "The job's command: in other words what you would type into a bash shell to run it normally.  Surround with single quotes or use backslahes to escape symbols to make sure it is parsed correctly.  The variable dollarfile will be populated from --arraygen, when the latter is used."))
	#--export=ALL
	
	return vars

def gen_argparser(argv, key):
	import distutils.util
	#Parse command line arguments (argv) using the information stored in the list of tuples.
	parser = argparse.ArgumentParser(description='This script generates a SLURM script for a desired command.  You must provide the job name and the command.  SLURM options can be changed to alter your typical #SBATCH variables.  The options below list the SLURM options being controlled in parentheses.')

	for argnum in range(0, len(key)):
		if key[argnum][2] == bool:
			parser.add_argument('--'+key[argnum][0], type=distutils.util.strtobool, default=key[argnum][3], required=key[argnum][4], help=key[argnum][6])
		else:
			parser.add_argument('--'+key[argnum][0], type=key[argnum][2], default=key[argnum][3], required=key[argnum][4], help=key[argnum][6])
	
	#--export=ALL
	args = parser.parse_args()
	return(args)
	
#This function should examine the input and determine if there are any problems.
#Hopefully return a helpful error message and quit if a problem is detected.
def validate(args):
	#If the user specifies --outfiles (identical prefixes for the log and the error files), they cannot specifiy either --output or --error.
	try:
		if args.outfiles != '' and (args.log != '' or args.err != ''):
			raise ArgError('You have provided both a value for both --outfiles and --log and/or --err.  If you provide --outfiles, you cannot provide --log or --err.')
	except ArgError:
		exit()
		
	#--mem is input as a string, as it can be listed as an integer in MB (eg. --mem 10000) or a value with a suffix (eg. --mem 10G).
	#Check to what is input is valid.
	try:
		if args.mem[-1] in ['K', 'M', 'G', 'T']:
			int(args.mem[:-1])
		else:
			int(args.mem)
	except ValueError:
		exit("--mem must be provided as an integer or an integer followed by the suffix K, M, G, or T for kilobyte, megabyte, gigabyte, or terrabyte.")
		
	#--openmode must be either None, append, or truncate.
	try:
		if args.openmode != None and args.openmode != 'append' and args.openmode != 'truncate':
			raise ArgError('The only valid options for --openmode are "append" or "truncate".')
	except ArgError:
		exit()
		
	#If using arrays, make sure %a is in output and error filenames, and '$job' is in the command script.
	try:
		if args.usearray:
			#Make sure --array and/or --arraygen are given.
			if args.arraygen == None and args.array == None:
				raise ArgError("You are creating an array submission without providing either --arraygen or --array.")
				
			#Check for %a in the logs.
			if (args.outfiles != '') and ('%a' not in args.outfiles):
				raise ArgError('You have created an array submission without putting %a in the --outfiles parameter.')
			if (args.log != '') and ('%a' not in args.log):
				raise ArgError('You have created an array submission without putting %a in the --log parameter.')
			if (args.err != '') and ('%a' not in args.err):
				raise ArgError('You have created an array submission without putting %a in the --err parameter.')
			
			#If we are using numeric arrays only, make sure $job is in --command.		
			#if args.arraygen == None and "$job" not in args.command:
			#	raise ArgError("You have created a numeric array submission (ie. without --arraygen) without putting '$job' in the --command parameter.")
				
			#If we are using file-based arrays, make sure $file is in --command.
			if args.arraygen != None and "$file" not in args.command:
				raise ArgError("You have created an array submission with --arraygen without putting '$file' in the --command parameter.")
				
			#If we are using special filename arrays, make sure arrayformat contains arr and [$job].
			if "arr" not in args.arrayformat:
				raise ArgError("You have given --arrayformat without referencing $arr.  The $arr variable is used to access the array elements, and must be included in --arrayformat.")
			if "[$job]" not in args.arrayformat:
				raise ArgError("You have given --arrayformat without referencing $job.  You must include $job for substitution of the array index.")
	except ArgError:
		exit()

#Make final variable values depending on the user input.
def arg_logic(args, key):
	#If outfiles, outfiles, log, and err are all empty, generate log and err based on jobname.
	if (args.outfiles == '' and args.log == '' and args.err == ''):
		args.log = args.job + '.%N.%j.out'
		args.err = args.job + '.%N.%j.err'
		
	#If outfiles is defined, set log and err based on outfiles.
	if (args.outfiles != ''):
		args.log = args.outfiles + '.%N.%j.out'
		args.err = args.outfiles + '.%N.%j.err'
		
	#If args.usearray, check if args.array is empty, and autogenerate it if so.
	if (args.usearray):
		if args.array == None:
			arraylen = len(glob.glob(args.arraygen))
			args.array = "0-" + str(arraylen-1)
			
	#If args.openmode is not defined, set it according to whether requeue is set.
	#If it is defined, use whatever the user said.
	if args.openmode == None:
		if args.requeue == True:
			args.openmode = 'append'
		else:
			args.openmode = 'truncate'

	return(args)
		
#This function creates and returns the SLURM script as a list, where each list item is a line in the file.
#It uses the command line arguments (args) and the list of tuples to figure out what to do with it all.
def make_script(args, key):
	#Add the #! and --export=ALL.
	script = ['#!/bin/bash\n', '#SBATCH --export=ALL\n']
	
	#Loop over all arguments in the key.
	for argnum in range(len(key)):
		#Do not add a line unless we have processed this variable and switched the flag on addline.
		addline = False
		#If this is a "regular" variable that takes an argument, process it this way.
		if key[argnum][5] == 'regular':
			line = "#SBATCH --" + key[argnum][1] + " " + str(getattr(args, key[argnum][0])) + "\n"
			addline = True
		#If this is an "array" variable, process it only if usearray is given.
		if key[argnum][5] == 'array' and args.usearray:
			line = "#SBATCH --" + key[argnum][1] + " " + str(getattr(args, key[argnum][0])) + "\n"
			addline = True
		#If this is a "boolean" variable that shows up or not in the SLURM script, but doesn't take an argment, process this way.
		if key[argnum][5] == 'boolean':
			if getattr(args, key[argnum][0]) == True:
				line = "#SBATCH --" + key[argnum][1] + "\n"
				addline = True
		
		#If we processed the line for this variable, add it to the script list of lines.
		if addline == True:
			script.append(line)
		
	#If we are using arrays, add the definition of job.
	if args.array != None:
		script.append("job=$SLURM_ARRAY_TASK_ID\n")
			
	#Add the command line to the bottom of the script.  If using a "special" array, run it with srun.
	if args.usearray and args.arraygen != None:
		script.append("arr=(" + args.arraygen + ")\n")
		script.append("file=" + args.arrayformat + "\n")
		script.append("srun " + args.command + "\n")
	else:
		script.append(args.command + "\n")
	
	return(script)
			
#Write all the data in script to the file scriptname.
def write_script(script, scriptname):
	outfile = open(scriptname, 'w')
	for line in script:
		outfile.write(line)
	outfile.close()
			
def main(argv):
	#Setup a list of tuples.  Each tuple should contain:
	#(0 arg.varname, 1 SLURM variable name, 2 variable type, 3 default, 4 required, 5 class)
	#class should be 'regular' for SLURM parameters that receive an argument, 'boolean' for those that either appear or do not,
	#'array' for array-related arguments, 'other' for other types of arguments.
	varkey = gen_varkey()

	#Generate the argument list using argparser
	args = gen_argparser(argv, varkey)
	
	#Validate arguments
	validate(args)
	
	#Generate "logical defaults" based on the input arguments.
	procargs = arg_logic(args, varkey)
	
	#Process the arguments and generate the slurm script.
	slurm_script = make_script(procargs, varkey)
	
	#Print the script to screen.
	for line in slurm_script:
		print(line)
		
	#Make a filename from the jobname, and then write the script.
	scriptname = 'tmp_' + args.job + '.sh'
	write_script(slurm_script, scriptname)
		
	#If desired, execute the script.
	if args.execute == True:
		#Run the script, waiting for it to finish.
		subprocess.call(['sbatch', scriptname])
		#subprocess.call(['cat', scriptname])
		
	#If desired, cleanup the script.
	if args.cleanup == True:
		#Delete the script.
		os.remove(scriptname)
		
if __name__ == "__main__":
	main(sys.argv[1:])