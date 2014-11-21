#!/usr/bin/env python
'''Instantiate the ENCODE ChIP-seq workflow'''

import pdb
import os.path, sys, subprocess, logging, re
import dxpy

EPILOG = '''Notes:

Examples:

	%(prog)s
'''

WF_TITLE = 'tf_chip_seq'
WF_DESCRIPTION = 'ENCODE TF ChIP-Seq Pipeline'

MAPPING_APPLET_NAME = 'encode_bwa'
FILTER_QC_APPLET_NAME = 'filter_qc'
XCOR_APPLET_NAME = 'xcor'
SPP_APPLET_NAME = 'spp'
POOL_APPLET_NAME = 'pool'
PSEUDOREPLICATOR_APPLET_NAME = 'pseudoreplicator'
ENCODE_SPP_APPLET_NAME = 'encode_spp'
IDR_APPLET_NAME='idr'

APPLETS = {}

def get_args():
	import argparse
	parser = argparse.ArgumentParser(
		description=__doc__, epilog=EPILOG,
		formatter_class=argparse.RawDescriptionHelpFormatter)

	parser.add_argument('--debug',   help="Print debug messages", 				default=False, action='store_true')
	parser.add_argument('--reference', help="Reference tar to map to", 			default='ENCODE Reference Files:male.hg19.fa.gz')
	parser.add_argument('--rep1',    help="Replicate 1 fastq.gz", 				default=None, nargs='*')
	parser.add_argument('--rep2',    help="Replicate 2 fastq.gz", 				default=None, nargs='*')
	parser.add_argument('--ctl1',    help="Control for Replicate 1 fastq.gz", 	default=None, nargs='*')
	parser.add_argument('--ctl2',    help="Control for Replicate 2 fastq.gz", 	default=None, nargs='*')
	parser.add_argument('--outp',    help="Output project name or ID", 			default=dxpy.WORKSPACE_ID)
	parser.add_argument('--outf',    help="Output folder name or ID", 			default="/analysis_run")
	parser.add_argument('--name',    help="Name of new workflow", 				default="TF ChIP-Seq")
	parser.add_argument('--applets', help="Name of project containing applets", default="E3 ChIP-seq")
	parser.add_argument('--instance_type', help="Instance type for applets",	default=None)
	parser.add_argument('--idr', 	 help='Run IDR analysis',					default=False, action='store_true')

	args = parser.parse_args()

	global DEBUG
	DEBUG = args.debug
	if DEBUG:
		logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
	else: #use the defaulf logging level
		logging.basicConfig(format='%(levelname)s:%(message)s')

	logging.debug("rep1 is: %s" %(args.rep1))

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
		folder_name = '/'
		file_name = file_identifier

	logging.debug("Looking for file %s in folder %s" %(file_name, folder_name))

	try:
		file_handler = dxpy.find_one_data_object(name=file_name, folder=folder_name, project=project.get_id(),
			more_ok=False, zero_ok=False, return_handler=True)
	except:
		logging.debug('%s not found in project %s folder %s' %(file_name, project.get_id(), folder_name))
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
		title=WF_TITLE,
		name=args.name,
		description=WF_DESCRIPTION,
		project=output_project.get_id(),
		folder=output_folder)

	blank_workflow = not (args.rep1 or args.rep2 or args.ctl1 or args.ctl2)

	#this whole strategy is fragile and unsatisfying
	#subsequent code assumes reps come before contols
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
				stage_input=mapping_stage_input,
				instance_type=args.instance_type
			)
			mapping_superstage.update({'map_stage_id': mapped_stage_id})

			filter_qc_stage_id = workflow.add_stage(
				filter_qc_applet,
				name='Filter_QC %s' %(superstage_name),
				folder=filter_qc_output_folder,
				stage_input={
					'input_bam': dxpy.dxlink({'stage': mapped_stage_id, 'outputField': 'mapped_reads'}),
					'paired_end': dxpy.dxlink({'stage': mapped_stage_id, 'outputField': 'paired_end'})
				},
				instance_type=args.instance_type
			)
			mapping_superstage.update({'filter_qc_stage_id': filter_qc_stage_id})

			xcor_stage_id = workflow.add_stage(
				xcor_applet,
				name='Xcor %s' %(superstage_name),
				folder=xcor_output_folder,
				stage_input={
					'input_bam': dxpy.dxlink({'stage': filter_qc_stage_id, 'outputField': 'filtered_bam'}),
					'paired_end': dxpy.dxlink({'stage': filter_qc_stage_id, 'outputField': 'paired_end'})
				},
				instance_type=args.instance_type
			)
			mapping_superstage.update({'xcor_stage_id': xcor_stage_id})

	spp_applet = find_applet_by_name(SPP_APPLET_NAME, applet_project.get_id())
	spp_stages = []
	peaks_output_folder = resolve_folder(output_project, output_folder + '/' + spp_applet.name)
	if (args.rep1 and args.ctl1) or blank_workflow:
		rep1_spp_stage_id = workflow.add_stage(
			spp_applet,
			name='Peaks Rep1',
			folder=peaks_output_folder,
			stage_input={
				'experiment' : dxpy.dxlink(
					{'stage': next(ss.get('xcor_stage_id') for ss in mapping_superstages if ss['name'] == 'Rep1'),
					 'outputField': 'tagAlign_file'}),
				'control': dxpy.dxlink(
					{'stage' : next(ss.get('xcor_stage_id') for ss in mapping_superstages if ss['name'] == 'Ctl1'),
					 'outputField': 'tagAlign_file'}),
				'xcor_scores_input': dxpy.dxlink(
					{'stage': next(ss.get('xcor_stage_id') for ss in mapping_superstages if ss['name'] == 'Rep1'),
					 'outputField': 'CC_scores_file'})
			},
			instance_type=args.instance_type
		)
		spp_stages.append({'name': 'Peaks Rep1', 'stage_id': rep1_spp_stage_id})
	if (args.rep2 and args.ctl2) or blank_workflow:
		rep2_spp_stage_id = workflow.add_stage(
			spp_applet,
			name='Peaks Rep2',
			folder=peaks_output_folder,
			stage_input={
				'experiment' : dxpy.dxlink(
					{'stage': next(ss.get('xcor_stage_id') for ss in mapping_superstages if ss['name'] == 'Rep2'),
					 'outputField': 'tagAlign_file'}),
				'control': dxpy.dxlink(
					{'stage' : next(ss.get('xcor_stage_id') for ss in mapping_superstages if ss['name'] == 'Ctl2'),
					 'outputField': 'tagAlign_file'}),
				'xcor_scores_input': dxpy.dxlink(
					{'stage': next(ss.get('xcor_stage_id') for ss in mapping_superstages if ss['name'] == 'Rep2'),
					 'outputField': 'CC_scores_file'})
			},
			instance_type=args.instance_type
		)
		spp_stages.append({'name': 'Peaks Rep2', 'stage_id': rep2_spp_stage_id})

	encode_spp_applet = find_applet_by_name(ENCODE_SPP_APPLET_NAME, applet_project.get_id())
	encode_spp_stages = []
	if args.idr:
		idr_peaks_output_folder = resolve_folder(output_project, output_folder + '/' + encode_spp_applet.name)
		if (args.rep1 and args.ctl1 and args.rep2 and args.ctl2) or blank_workflow:
			encode_spp_stage_id = workflow.add_stage(
				encode_spp_applet,
				name='Peaks for IDR',
				folder=idr_peaks_output_folder,
				stage_input={
					'rep1_ta' : dxpy.dxlink(
						{'stage': next(ss.get('xcor_stage_id') for ss in mapping_superstages if ss['name'] == 'Rep1'),
						 'outputField': 'tagAlign_file'}),
					'rep2_ta' : dxpy.dxlink(
						{'stage': next(ss.get('xcor_stage_id') for ss in mapping_superstages if ss['name'] == 'Rep2'),
						 'outputField': 'tagAlign_file'}),
					'ctl1_ta': dxpy.dxlink(
						{'stage': next(ss.get('xcor_stage_id') for ss in mapping_superstages if ss['name'] == 'Ctl1'),
						 'outputField': 'tagAlign_file'}),
					'ctl2_ta' : dxpy.dxlink(
						{'stage': next(ss.get('xcor_stage_id') for ss in mapping_superstages if ss['name'] == 'Ctl2'),
						 'outputField': 'tagAlign_file'}),
					'rep1_xcor' : dxpy.dxlink(
						{'stage': next(ss.get('xcor_stage_id') for ss in mapping_superstages if ss['name'] == 'Rep1'),
						 'outputField': '"CC_scores_file"'}),
					'rep2_xcor' : dxpy.dxlink(
						{'stage': next(ss.get('xcor_stage_id') for ss in mapping_superstages if ss['name'] == 'Rep2'),
						 'outputField': '"CC_scores_file"'}),
					'paired_end': dxpy.dxlink(
						{'stage': next(ss.get('xcor_stage_id') for ss in mapping_superstages if ss['name'] == 'Rep1'),
						 'outputField': 'paired_end'}) #here we're assuming if rep1 is PE it's a PE experiment - need better error checking
				},
				instance_type=args.instance_type
			)
			encode_spp_stages.append({'name': 'Peaks for IDR', 'stage_id': encode_spp_stage_id})

	idr_applet = find_applet_by_name(IDR_APPLET_NAME, applet_project.get_id())
	idr_stages = []
	if args.idr:
		idr_output_folder = resolve_folder(output_project, output_folder + '/' + idr_applet.name)
		if (args.rep1 and args.ctl1 and args.rep2 and args.ctl2) or blank_workflow:
			idr_stage_id = workflow.add_stage(
				idr_applet,
				name='IDR True Replicates',
				folder=idr_output_folder,
				stage_input={
					'rep1_peaks' : dxpy.dxlink(
						{'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == 'Peaks for IDR'),
						 'outputField': 'rep1_peaks'}),
					'rep2_peaks' : dxpy.dxlink(
						{'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == 'Peaks for IDR'),
						 'outputField': 'rep2_peaks'}),
					'pooled_peaks': dxpy.dxlink(
						{'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == 'Peaks for IDR'),
						 'outputField': 'pooled_peaks'})
				},
				instance_type=args.instance_type
			)
			idr_stages.append({'name': 'IDR True Replicates', 'stage_id': idr_stage_id})

			idr_stage_id = workflow.add_stage(
				idr_applet,
				name='IDR Rep 1 Self-pseudoreplicates',
				folder=idr_output_folder,
				stage_input={
					'rep1_peaks' : dxpy.dxlink(
						{'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == 'Peaks for IDR'),
						 'outputField': 'rep1pr1_peaks'}),
					'rep2_peaks' : dxpy.dxlink(
						{'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == 'Peaks for IDR'),
						 'outputField': 'rep1pr2_peaks'}),
					'pooled_peaks': dxpy.dxlink(
						{'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == 'Peaks for IDR'),
						 'outputField': 'rep1_peaks'})
				},
				instance_type=args.instance_type
			)
			idr_stages.append({'name': 'IDR Rep 1 Self-pseudoreplicates', 'stage_id': idr_stage_id})

			idr_stage_id = workflow.add_stage(
				idr_applet,
				name='IDR Rep 2 Self-pseudoreplicates',
				folder=idr_output_folder,
				stage_input={
					'rep1_peaks' : dxpy.dxlink(
						{'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == 'Peaks for IDR'),
						 'outputField': 'rep2pr1_peaks'}),
					'rep2_peaks' : dxpy.dxlink(
						{'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == 'Peaks for IDR'),
						 'outputField': 'rep2pr2_peaks'}),
					'pooled_peaks': dxpy.dxlink(
						{'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == 'Peaks for IDR'),
						 'outputField': 'rep2_peaks'})
				},
				instance_type=args.instance_type
			)
			idr_stages.append({'name': 'IDR Rep 2 Self-pseudoreplicates', 'stage_id': idr_stage_id})

			idr_stage_id = workflow.add_stage(
				idr_applet,
				name='IDR Pooled Pseudoeplicates',
				folder=idr_output_folder,
				stage_input={
					'rep1_peaks' : dxpy.dxlink(
						{'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == 'Peaks for IDR'),
						 'outputField': 'pooledpr1_peaks'}),
					'rep2_peaks' : dxpy.dxlink(
						{'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == 'Peaks for IDR'),
						 'outputField': 'pooledpr2_peaks'}),
					'pooled_peaks': dxpy.dxlink(
						{'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == 'Peaks for IDR'),
						 'outputField': 'pooled_peaks'})
				},
				instance_type=args.instance_type
			)
			idr_stages.append({'name': 'IDR Pooled Pseudoeplicates', 'stage_id': idr_stage_id})


	logging.debug("Mapping stages: %s" %(mapping_superstages))
	logging.debug("Peak stages: %s" %(spp_stages))
	logging.debug("Peaks for IDR stages: %s" %(encode_spp_stages))
	logging.debug("IDR stages: %s" %(idr_stages))

if __name__ == '__main__':
	main()
