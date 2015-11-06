#!/usr/bin/env python
'''Instantiate the ENCODE ChIP-seq workflow'''

import pdb
import os.path, sys, subprocess, logging, re
import dxpy

EPILOG = '''Notes:

Examples:
	Build blank workflow from fastq to peaks (no IDR)
	%(prog)s

	Build a blank workflow that includes both naive peak calling and IDR.
	%(prog)s --idr

	Build and run a workflow, specifying fastq's for two replicates and matched controls, including naive peaks and IDR.
	%(prog)s --rep1 r1.fastq.gz --rep2 r2.fastq.gz --ctl1 c1.fastq.gz --ctl2 c2.fastq.gz --idr --yes

	Build and run a workflow, specifying fastq's for two replicates and matched controls, reporting only IDR-processed peaks.
	%(prog)s --rep1 r1.fastq.gz --rep2 r2.fastq.gz --ctl1 c1.fastq.gz --ctl2 c2.fastq.gz --idronly --yes

	Build and run a workflow, skipping mapping and starting from tagAligns from paired-end data, reporting both naive and IDR-processed peaks.
	%(prog)s --rep1 f1.tagAlign.gz --rep2 r2.tagAlign.gz --ctl1 c1.tagAlign.gz --ctl2 c2.tagAlign.gz --rep1_ended PE --rep2_ended PE --idr --yes

'''

WF_NAME = 'tf_chip_seq'
WF_DESCRIPTION = 'ENCODE TF ChIP-Seq Pipeline'

MAPPING_APPLET_NAME = 'encode_bwa'
FILTER_QC_APPLET_NAME = 'filter_qc'
XCOR_APPLET_NAME = 'xcor'
XCOR_ONLY_APPLET_NAME = 'xcor_only'
SPP_APPLET_NAME = 'spp'
POOL_APPLET_NAME = 'pool'
PSEUDOREPLICATOR_APPLET_NAME = 'pseudoreplicator'
ENCODE_SPP_APPLET_NAME = 'encode_spp'
ENCODE_MACS2_APPLET_NAME = 'encode_macs2'
IDR_APPLET_NAME='idr'
IDR2_APPLET_NAME='idr2'
ENCODE_IDR_APPLET_NAME='encode_idr'

APPLETS = {}

def get_args():
	import argparse
	parser = argparse.ArgumentParser(
		description=__doc__, epilog=EPILOG,
		formatter_class=argparse.RawDescriptionHelpFormatter)

	parser.add_argument('--debug',   help="Print debug messages and hold jobs for ssh", 				default=False, action='store_true')
	parser.add_argument('--reference', help="Reference tar to map to")
	parser.add_argument('--chrom_sizes', help="chrom.sizes file for bed to bigbed")
	parser.add_argument('--genomesize', help="Genome size string for MACS2, e.g. mm or hs")
	parser.add_argument('--narrowpeak_as', help=".as file for bed to bigbed", default='ENCODE Reference Files:narrowPeak.as')
	parser.add_argument('--gappedpeak_as', help=".as file for bed to bigbed", default='ENCODE Reference Files:gappedPeak.as')
	parser.add_argument('--broadpeak_as', help=".as file for bed to bigbed", default='ENCODE Reference Files:broadPeak.as')
	parser.add_argument('--rep1',    help="Replicate 1 fastq or tagAlign", 				default=None, nargs='*')
	parser.add_argument('--rep2',    help="Replicate 2 fastq or tagAlign", 				default=None, nargs='*')
	parser.add_argument('--ctl1',    help="Control for replicate 1 fastq or tagAlign", 	default=None, nargs='*')
	parser.add_argument('--ctl2',    help="Control for replicate 2 fastq or tagAlign", 	default=None, nargs='*')
	parser.add_argument('--outp',    help="Output project name or ID", 			default=dxpy.WORKSPACE_ID)
	parser.add_argument('--outf',    help="Output folder name or ID", 			default="/analysis_run")
	parser.add_argument('--name',    help="Title of new workflow", 				default="TF ChIP-Seq")
	parser.add_argument('--applets', help="Name of project containing applets", default="E3 ChIP-seq")
	parser.add_argument('--nomap',   help='Given tagAligns, skip to peak calling', default=False, action='store_true')
	parser.add_argument('--rep1pe', help='Specify if rep1 is PE (required only if --nomap)', type=bool, default=None)
	parser.add_argument('--rep2pe', help='Specify if rep2 is PE (required only if --nomap)', type=bool, default=None)
	parser.add_argument('--blacklist', help="Blacklist to filter IDR peaks")
	parser.add_argument('--idr', 	 help='Report peaks with and without IDR analysis',					default=False, action='store_true')
	#parser.add_argument('--idronly',  help='Only report IDR peaks', default=None, action='store_true')
	parser.add_argument('--idrversion', help='Version of IDR to use (1 or 2)', default=2)
	parser.add_argument('--yes', 	 help='Run the workflow',					default=False, action='store_true')

	args = parser.parse_args()

	global DEBUG
	DEBUG = args.debug
	if DEBUG:
		logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
	else: #use the defaulf logging level
		logging.basicConfig(format='%(levelname)s:%(message)s')

	logging.debug("rep1 is: %s" %(args.rep1))

	if args.nomap and (args.rep1pe == None or args.rep2pe == None):
		logging.error("With --nomap, endedness of replicates must be specified with --rep1pe and --rep2pe")
		raise ValueError

	return args

def blank_workflow(args):
	return

def map_and_filter(infile, args):
	if not infile:
		return {None}
	
	stages = {None}

	return stages

def call_peaks(expvsctl, args):
	if not expvsctl:
		return {None}

	stages = {None}

	return stages

def resolve_project(identifier, privs='r'):
	project = dxpy.find_one_project(name=identifier, level='VIEW', name_mode='exact', return_handler=True, zero_ok=True)
	if project == None:
		try:
			project = dxpy.get_handler(identifier)
		except:
			logging.error('Could not find a unique project with name or id %s' %(identifier))
			raise ValueError(identifier)
	logging.debug('Project %s access level is %s' %(project.name, project.describe()['level']))
	if privs == 'w' and project.describe()['level'] == 'VIEW':
		logging.error('Output project %s is read-only' %(identifier))
		raise ValueError(identifier)
	return project

def resolve_folder(project, identifier):
	if not identifier.startswith('/'):
		identifier = '/' + identifier
	try:
		project_id = project.list_folder(identifier)
	except:
		try:
			project_id = project.new_folder(identifier, parents=True)
		except:
			logging.error("Cannot create folder %s in project %s" %(identifier, project.name))
			raise ValueError('%s:%s' %(project.name, identifier))
		else:
			logging.info("New folder %s created in project %s" %(identifier, project.name))
	return identifier

def resolve_file(identifier):
	logging.debug("resolve_file: %s" %(identifier))

	if not identifier:
		return None

	m = re.match(r'''^([\w\-\ \.]+):([\w\-\ /\.]+)''', identifier)
	if m:
		project_identifier = m.group(1)
		file_identifier = m.group(2)
	else:
		logging.debug("Defaulting to the current project")
		project_identifier = dxpy.WORKSPACE_ID
		file_identifier = identifier	

	project = resolve_project(project_identifier)
	logging.debug("Got project %s" %(project.name))
	logging.debug("Now looking for file %s" %(file_identifier))

	m = re.match(r'''(^[\w\-\ /\.]+)/([\w\-\ \.]+)''', file_identifier)
	if m:
		folder_name = m.group(1)
		if not folder_name.startswith('/'):
			folder_name = '/' + folder_name
		file_name = m.group(2)
	else:
		# folder_name = '/'
		folder_name = None
		file_name = file_identifier

	logging.debug("Looking for file %s in folder %s" %(file_name, folder_name))

	try:
		if folder_name:
			file_handler = dxpy.find_one_data_object(name=file_name, folder=folder_name, project=project.get_id(),
				recurse=False, more_ok=False, zero_ok=False, return_handler=True)
		else:
			file_handler = dxpy.find_one_data_object(name=file_name, project=project.get_id(), folder='/',
				recurse=True, more_ok=False, zero_ok=False, return_handler=True)
	except dxpy.DXSearchError:
		logging.debug('%s not found in project %s folder %s.  Trying as file ID' %(file_name, project.get_id(), folder_name))
		try:
			file_handler = dxpy.DXFile(dxid=identifier, mode='r')
		except:
			logging.debug('%s not found as a dxid' %(identifier))
			try:
				file_handler = resolve_accession(identifier)
			except:
				logging.debug('%s not found as an accession' %(identifier))
				logging.warning('Could not find file %s.' %(identifier))
				return None
	except:
		raise
	logging.info("Resolved file identifier %s to %s" %(identifier, file_handler.get_id()))
	return file_handler

def resolve_accession(accession):
	logging.debug("Looking for accession %s" %(accession))
	
	if not re.match(r'''^ENCFF\d{3}[A-Z]{3}''', accession):
		logging.debug("%s is not a valid accession format" %(accession))
		raise ValueError(accession)
	
	DNANEXUS_ENCODE_SNAPSHOT = 'ENCODE-SDSC-snapshot-20140505'
	logging.debug('Testing')

	try:
		snapshot_project
	except:
		logging.debug('Looking for snapshot project %s' %(DNANEXUS_ENCODE_SNAPSHOT))
		try:
			project_handler = resolve_project(DNANEXUS_ENCODE_SNAPSHOT)
			global snapshot_project
			snapshot_project = project_handler
		except:
			logging.error("Cannot find snapshot project %s" %(DNANEXUS_ENCODE_SNAPSHOT))
			raise ValueError(DNANEXUS_ENCODE_SNAPSHOT)
		logging.debug('Found snapshot project %s' %(snapshot_project.name))

	try:
		accession_search = accession + '*'
		logging.debug('Looking recursively for %s in %s' %(accession_search, snapshot_project.name))
		file_handler = dxpy.find_one_data_object(
			name=accession_search, name_mode='glob', more_ok=False, classname='file', recurse=True, return_handler=True,
			folder='/', project=snapshot_project.get_id())
		logging.debug('Got file handler for %s' %(file_handler.name))
		return file_handler
	except:
		logging.error("Cannot find accession %s in project %s" %(accession, snapshot_project.name))
		raise ValueError(accession)

def find_applet_by_name(applet_name, applets_project_id):
	'''Looks up an applet by name in the project that holds tools.  From Joe Dale's code.'''
	cached = '*'
	if (applet_name, applets_project_id) not in APPLETS:
		found = dxpy.find_one_data_object(classname="applet", name=applet_name,
										  project=applets_project_id,
										  zero_ok=False, more_ok=False, return_handler=True)
		APPLETS[(applet_name, applets_project_id)] = found
		cached = ''

	logging.info(cached + "Resolved applet %s to %s" %(applet_name, APPLETS[(applet_name, applets_project_id)].get_id()))
	return APPLETS[(applet_name, applets_project_id)]

def main():
	args = get_args()

	output_project = resolve_project(args.outp, 'w')
	logging.info('Found output project %s' %(output_project.name))
	output_folder = resolve_folder(output_project, args.outf)
	logging.info('Using output folder %s' %(output_folder))
	applet_project = resolve_project(args.applets, 'r')
	logging.info('Found applet project %s' %(applet_project.name))

	workflow = dxpy.new_dxworkflow(
		name=WF_NAME,
		title=args.name,
		description=WF_DESCRIPTION,
		project=output_project.get_id(),
		folder=output_folder)

	blank_workflow = not (args.rep1 or args.rep2 or args.ctl1 or args.ctl2)

	if not args.nomap:
		#a "superstage" is just a dict with a name, name(s) of input files, and then names and id's of stages that process that input
		#each superstage here could be implemented as a stage in a more abstract workflow.  That stage would then call the various applets that are separate
		#stages here.
		mapping_superstages = [
			{'name': 'Rep1', 'input_args': args.rep1},
			{'name': 'Rep2', 'input_args': args.rep2},
			{'name': 'Ctl1', 'input_args': args.ctl1},
			{'name': 'Ctl2', 'input_args': args.ctl2}
			# {'name': 'Pooled Reps', 'input_args': (args.rep1 and args.rep2)},
			# {'name': 'Pooled Controls', 'input_args': (args.ctl1 and args.ctl2)} ##idea is to create a "stub" stage and then populate it's input with the output of the pool stage, defined below
		]

		mapping_applet = find_applet_by_name(MAPPING_APPLET_NAME, applet_project.get_id())
		mapping_output_folder = resolve_folder(output_project, output_folder + '/' + mapping_applet.name)
		reference_tar = resolve_file(args.reference)
		filter_qc_applet = find_applet_by_name(FILTER_QC_APPLET_NAME, applet_project.get_id())
		filter_qc_output_folder = mapping_output_folder
		xcor_applet = find_applet_by_name(XCOR_APPLET_NAME, applet_project.get_id())
		xcor_output_folder = mapping_output_folder

		for mapping_superstage in mapping_superstages:
			superstage_name = mapping_superstage.get('name')

			if mapping_superstage.get('input_args') or blank_workflow:
				if blank_workflow:
					mapping_stage_input = None
				else:
					mapping_stage_input = {'reference_tar' : dxpy.dxlink(reference_tar.get_id())}
					for arg_index,input_arg in enumerate(mapping_superstage['input_args']): #read pairs assumed be in order read1,read2
						reads = dxpy.dxlink(resolve_file(input_arg).get_id())
						mapping_stage_input.update({'reads%d' %(arg_index+1): reads})

				mapped_stage_id = workflow.add_stage(
					mapping_applet,
					name='Map %s' %(superstage_name),
					folder=mapping_output_folder,
					stage_input=mapping_stage_input
				)
				mapping_superstage.update({'map_stage_id': mapped_stage_id})

				filter_qc_stage_id = workflow.add_stage(
					filter_qc_applet,
					name='Filter_QC %s' %(superstage_name),
					folder=filter_qc_output_folder,
					stage_input={
						'input_bam': dxpy.dxlink({'stage': mapped_stage_id, 'outputField': 'mapped_reads'}),
						'paired_end': dxpy.dxlink({'stage': mapped_stage_id, 'outputField': 'paired_end'})
					}
				)
				mapping_superstage.update({'filter_qc_stage_id': filter_qc_stage_id})

				xcor_stage_id = workflow.add_stage(
					xcor_applet,
					name='Xcor %s' %(superstage_name),
					folder=xcor_output_folder,
					stage_input={
						'input_bam': dxpy.dxlink({'stage': filter_qc_stage_id, 'outputField': 'filtered_bam'}),
						'paired_end': dxpy.dxlink({'stage': filter_qc_stage_id, 'outputField': 'paired_end'})
					}
				)
				mapping_superstage.update({'xcor_stage_id': xcor_stage_id})

		exp_rep1_ta = dxpy.dxlink(
					{'stage': next(ss.get('xcor_stage_id') for ss in mapping_superstages if ss['name'] == 'Rep1'),
					 'outputField': 'tagAlign_file'})
		exp_rep1_cc = dxpy.dxlink(
					{'stage': next(ss.get('xcor_stage_id') for ss in mapping_superstages if ss['name'] == 'Rep1'),
					 'outputField': 'CC_scores_file'})
		exp_rep2_ta = dxpy.dxlink(
					{'stage': next(ss.get('xcor_stage_id') for ss in mapping_superstages if ss['name'] == 'Rep2'),
					 'outputField': 'tagAlign_file'})
		exp_rep2_cc = dxpy.dxlink(
					{'stage': next(ss.get('xcor_stage_id') for ss in mapping_superstages if ss['name'] == 'Rep2'),
					 'outputField': 'CC_scores_file'})
		ctl_rep1_ta = dxpy.dxlink(
					{'stage' : next(ss.get('xcor_stage_id') for ss in mapping_superstages if ss['name'] == 'Ctl1'),
					 'outputField': 'tagAlign_file'})
		ctl_rep2_ta = dxpy.dxlink(
					{'stage' : next(ss.get('xcor_stage_id') for ss in mapping_superstages if ss['name'] == 'Ctl2'),
					 'outputField': 'tagAlign_file'})
		rep1_paired_end = dxpy.dxlink(
						{'stage': next(ss.get('xcor_stage_id') for ss in mapping_superstages if ss['name'] == 'Rep1'),
						 'outputField': 'paired_end'})
		rep2_paired_end = dxpy.dxlink(
						{'stage': next(ss.get('xcor_stage_id') for ss in mapping_superstages if ss['name'] == 'Rep2'),
						 'outputField': 'paired_end'})
	else: #skipped the mapping, so just bring in the inputs from arguments
		exp_rep1_ta = dxpy.dxlink(resolve_file(args.rep1[0]).get_id())
		exp_rep2_ta = dxpy.dxlink(resolve_file(args.rep2[0]).get_id())
		ctl_rep1_ta = dxpy.dxlink(resolve_file(args.ctl1[0]).get_id())
		ctl_rep2_ta = dxpy.dxlink(resolve_file(args.ctl2[0]).get_id())
		rep1_paired_end = args.rep1pe
		rep2_paired_end = args.rep2pe

		#here we need to calculate the cc scores files, because we're only being supplied tagAligns
		#if we had mapped everything above we'd already have a handle to the cc file
		xcor_only_applet = find_applet_by_name(XCOR_ONLY_APPLET_NAME, applet_project.get_id())
		xcor_output_folder = resolve_folder(output_project, output_folder + '/' + xcor_only_applet.name)
		xcor_only_stages = []

		exp_rep1_cc_stage_id = workflow.add_stage(
			xcor_only_applet,
			name="Rep1 cross-correlation",
			folder=xcor_output_folder,
			stage_input={
				'input_tagAlign': exp_rep1_ta,
				'paired_end': rep1_paired_end
			}
		)
		xcor_only_stages.append({'xcor_only_rep1_id': exp_rep1_cc_stage_id})
		exp_rep1_cc = dxpy.dxlink(
					{'stage': exp_rep1_cc_stage_id,
					 'outputField': 'CC_scores_file'})

		exp_rep2_cc_stage_id = workflow.add_stage(
			xcor_only_applet,
			name="Rep2 cross-correlation",
			folder=xcor_output_folder,
			stage_input={
				'input_tagAlign': exp_rep2_ta,
				'paired_end': rep2_paired_end
			}
		)
		xcor_only_stages.append({'xcor_only_rep2_id': exp_rep2_cc_stage_id})
		exp_rep2_cc = dxpy.dxlink(
					{'stage': exp_rep2_cc_stage_id,
					 'outputField': 'CC_scores_file'})

	# if not args.idronly:
	# 	spp_applet = find_applet_by_name(SPP_APPLET_NAME, applet_project.get_id())
	# 	peaks_output_folder = resolve_folder(output_project, output_folder + '/' + spp_applet.name)
	# 	spp_stages = []
	# 	if (args.rep1 and args.ctl1) or blank_workflow:
	# 		rep1_spp_stage_id = workflow.add_stage(
	# 			spp_applet,
	# 			name='Peaks Rep1',
	# 			folder=peaks_output_folder,
	# 			stage_input={
	# 				'experiment': exp_rep1_ta,
	# 				'control': ctl_rep1_ta,
	# 				'xcor_scores_input': exp_rep1_cc,
	# 				'bigbed': True,
	# 				'chrom_sizes': dxpy.dxlink(resolve_file(args.chrom_sizes)),
	# 				'as_file': dxpy.dxlink(resolve_file(args.as_file))
	# 			}

	# 		)
	# 		spp_stages.append({'name': 'Peaks Rep1', 'stage_id': rep1_spp_stage_id})
	# 	if (args.rep2 and args.ctl2) or blank_workflow:
	# 		rep2_spp_stage_id = workflow.add_stage(
	# 			spp_applet,
	# 			name='Peaks Rep2',
	# 			folder=peaks_output_folder,
	# 			stage_input={
	# 				'experiment': exp_rep2_ta,
	# 				'control': ctl_rep2_ta,
	# 				'xcor_scores_input': exp_rep2_cc,
	# 				'bigbed': True,
	# 				'chrom_sizes': dxpy.dxlink(resolve_file(args.chrom_sizes)),
	# 				'as_file': dxpy.dxlink(resolve_file(args.as_file))
	# 			}
	# 		)
	# 		spp_stages.append({'name': 'Peaks Rep2', 'stage_id': rep2_spp_stage_id})


	encode_spp_applet = find_applet_by_name(ENCODE_SPP_APPLET_NAME, applet_project.get_id())
	encode_spp_stages = []
	idr_peaks_output_folder = resolve_folder(output_project, output_folder + '/' + encode_spp_applet.name)
	PEAKS_STAGE_NAME = 'SPP Peaks'
	if (args.rep1 and args.ctl1 and args.rep2 and args.ctl2) or blank_workflow:
		encode_spp_stage_id = workflow.add_stage(
			encode_spp_applet,
			name=PEAKS_STAGE_NAME,
			folder=idr_peaks_output_folder,
			stage_input={
				'rep1_ta' : exp_rep1_ta,
				'rep2_ta' : exp_rep2_ta,
				'ctl1_ta': ctl_rep1_ta,
				'ctl2_ta' : ctl_rep2_ta,
				'rep1_xcor' : exp_rep1_cc,
				'rep2_xcor' : exp_rep2_cc,
				'rep1_paired_end': rep1_paired_end,
				'rep2_paired_end': rep2_paired_end,
				'chrom_sizes': dxpy.dxlink(resolve_file(args.chrom_sizes)),
				'as_file': dxpy.dxlink(resolve_file(args.narrowpeak_as)),
				'idr_peaks': args.idr
				}
			)

	encode_spp_stages.append({'name': PEAKS_STAGE_NAME, 'stage_id': encode_spp_stage_id})

	encode_macs2_applet = find_applet_by_name(ENCODE_MACS2_APPLET_NAME, applet_project.get_id())
	encode_macs2_stages = []
	peaks_output_folder = resolve_folder(output_project, output_folder + '/' + encode_macs2_applet.name)
	if (args.rep1 and args.ctl1 and args.rep2 and args.ctl2) or blank_workflow:
		encode_macs2_stage_id = workflow.add_stage(
			encode_macs2_applet,
			name='ENCODE Peaks',
			folder=peaks_output_folder,
			stage_input={
				'rep1_ta' : exp_rep1_ta,
				'rep2_ta' : exp_rep2_ta,
				'ctl1_ta': ctl_rep1_ta,
				'ctl2_ta' : ctl_rep2_ta,
				'rep1_xcor' : exp_rep1_cc,
				'rep2_xcor' : exp_rep2_cc,
				'rep1_paired_end': rep1_paired_end,
				'rep2_paired_end': rep2_paired_end,
				'chrom_sizes': dxpy.dxlink(resolve_file(args.chrom_sizes)),
				'narrowpeak_as': dxpy.dxlink(resolve_file(args.narrowpeak_as)),
				'gappedpeak_as': dxpy.dxlink(resolve_file(args.gappedpeak_as)),
				'broadpeak_as':  dxpy.dxlink(resolve_file(args.broadpeak_as)),
				'genomesize': args.genomesize
			}
		)
		encode_macs2_stages.append({'name': 'ENCODE Peaks', 'stage_id': encode_macs2_stage_id})

	if args.idr:
		if args.idrversion == "1":
			idr_applet = find_applet_by_name(IDR_APPLET_NAME, applet_project.get_id())
		elif args.idrversion == "2":
			idr_applet = find_applet_by_name(IDR2_APPLET_NAME, applet_project.get_id())
		else:
			logging.error("Invalid IDR version: %s" %(args.idrversion))
			idr_applet = None
		encode_idr_applet = find_applet_by_name(ENCODE_IDR_APPLET_NAME, applet_project.get_id())
		idr_stages = []
		idr_output_folder = resolve_folder(output_project, output_folder + '/' + idr_applet.name)
		if (args.rep1 and args.ctl1 and args.rep2 and args.ctl2) or blank_workflow:
			idr_stage_id = workflow.add_stage(
				idr_applet,
				name='IDR True Replicates',
				folder=idr_output_folder,
				stage_input={
					'rep1_peaks' : dxpy.dxlink(
						{'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == PEAKS_STAGE_NAME),
						 'outputField': 'rep1_peaks'}),
					'rep2_peaks' : dxpy.dxlink(
						{'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == PEAKS_STAGE_NAME),
						 'outputField': 'rep2_peaks'}),
					'pooled_peaks': dxpy.dxlink(
						{'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == PEAKS_STAGE_NAME),
						 'outputField': 'pooled_peaks'})
				}
			)
			idr_stages.append({'name': 'IDR True Replicates', 'stage_id': idr_stage_id})

			idr_stage_id = workflow.add_stage(
				idr_applet,
				name='IDR Rep 1 Self-pseudoreplicates',
				folder=idr_output_folder,
				stage_input={
					'rep1_peaks' : dxpy.dxlink(
						{'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == PEAKS_STAGE_NAME),
						 'outputField': 'rep1pr1_peaks'}),
					'rep2_peaks' : dxpy.dxlink(
						{'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == PEAKS_STAGE_NAME),
						 'outputField': 'rep1pr2_peaks'}),
					'pooled_peaks': dxpy.dxlink(
						{'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == PEAKS_STAGE_NAME),
						 'outputField': 'rep1_peaks'})
				}
			)
			idr_stages.append({'name': 'IDR Rep 1 Self-pseudoreplicates', 'stage_id': idr_stage_id})

			idr_stage_id = workflow.add_stage(
				idr_applet,
				name='IDR Rep 2 Self-pseudoreplicates',
				folder=idr_output_folder,
				stage_input={
					'rep1_peaks' : dxpy.dxlink(
						{'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == PEAKS_STAGE_NAME),
						 'outputField': 'rep2pr1_peaks'}),
					'rep2_peaks' : dxpy.dxlink(
						{'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == PEAKS_STAGE_NAME),
						 'outputField': 'rep2pr2_peaks'}),
					'pooled_peaks': dxpy.dxlink(
						{'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == PEAKS_STAGE_NAME),
						 'outputField': 'rep2_peaks'})
				}
			)
			idr_stages.append({'name': 'IDR Rep 2 Self-pseudoreplicates', 'stage_id': idr_stage_id})

			idr_stage_id = workflow.add_stage(
				idr_applet,
				name='IDR Pooled Pseudoeplicates',
				folder=idr_output_folder,
				stage_input={
					'rep1_peaks' : dxpy.dxlink(
						{'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == PEAKS_STAGE_NAME),
						 'outputField': 'pooledpr1_peaks'}),
					'rep2_peaks' : dxpy.dxlink(
						{'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == PEAKS_STAGE_NAME),
						 'outputField': 'pooledpr2_peaks'}),
					'pooled_peaks': dxpy.dxlink(
						{'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == PEAKS_STAGE_NAME),
						 'outputField': 'pooled_peaks'})
				}
			)
			idr_stages.append({'name': 'IDR Pooled Pseudoreplicates', 'stage_id': idr_stage_id})

			blacklist = resolve_file(args.blacklist)
			idr_stage_id = workflow.add_stage(
				encode_idr_applet,
				name='Final IDR peak calls',
				folder=idr_output_folder,
				stage_input={
					'reps_peaks' : dxpy.dxlink(
						{'stage': next(ss.get('stage_id') for ss in idr_stages if ss['name'] == 'IDR True Replicates'),
						 'outputField': 'IDR_peaks'}),
					'r1pr_peaks' : dxpy.dxlink(
						{'stage': next(ss.get('stage_id') for ss in idr_stages if ss['name'] == 'IDR Rep 1 Self-pseudoreplicates'),
						 'outputField': 'IDR_peaks'}),
					'r2pr_peaks' : dxpy.dxlink(
						{'stage': next(ss.get('stage_id') for ss in idr_stages if ss['name'] == 'IDR Rep 2 Self-pseudoreplicates'),
						 'outputField': 'IDR_peaks'}),
					'pooledpr_peaks': dxpy.dxlink(
						{'stage': next(ss.get('stage_id') for ss in idr_stages if ss['name'] == 'IDR Pooled Pseudoreplicates'),
						 'outputField': 'IDR_peaks'}),
					'blacklist': dxpy.dxlink(blacklist.get_id()),
					'chrom_sizes': dxpy.dxlink(resolve_file(args.chrom_sizes)),
					'as_file': dxpy.dxlink(resolve_file(args.narrowpeak_as))
				}
			)
			idr_stages.append({'name': 'Final IDR peak calls', 'stage_id': idr_stage_id})

	if not (args.nomap):
		logging.debug("Mapping stages: %s" %(mapping_superstages))
	else:
		logging.debug("xcor only stages: %s" %(xcor_only_stages))
	# if not args.idronly:
	# 	logging.debug("Peak stages: %s" %(spp_stages))
	logging.debug("Peak stages: %s" %(encode_spp_stages))
	if args.idr:
		logging.debug("IDR stages: %s" %(idr_stages))

	if args.yes:
		if args.debug:
			job_id = workflow.run({}, priority='high', debug={'debugOn': ['AppInternalError', 'AppError']}, delay_workspace_destruction=True, allow_ssh=['255.255.255.255'])
		else:
			job_id = workflow.run({}, priority='high')
		logging.info("Running as job %s" %(job_id))

if __name__ == '__main__':
	main()
