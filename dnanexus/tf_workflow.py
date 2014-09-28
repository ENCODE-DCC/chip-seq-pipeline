#!/usr/bin/env python
'''Instantiate the ENCODE ChIP-seq workflow'''

import pdb
import os.path, sys, subprocess, logging, re
import dxpy

EPILOG = '''Notes:

Examples:

	%(prog)s
'''

def get_args():
	import argparse
	parser = argparse.ArgumentParser(
		description=__doc__, epilog=EPILOG,
		formatter_class=argparse.RawDescriptionHelpFormatter)

	parser.add_argument('--debug',   help="Print debug messages", 				default=False, action='store_true')
	parser.add_argument('--rep1',    help="Replicate 1 fastq.gz", 				default=None)
	parser.add_argument('--rep2',    help="Replicate 2 fastq.gz", 				default=None)
	parser.add_argument('--ctl1',    help="Control for Replicate 1 fastq.gz", 	default=None)
	parser.add_argument('--ctl2',    help="Control for Replicate 2 fastq.gz", 	default=None)
	parser.add_argument('--outp',    help="Output project name or ID", 			default=dxpy.WORKSPACE_ID)
	parser.add_argument('--outf',    help="Output folder name or ID", 			default="/")
	parser.add_argument('--name',    help="Name of new workflow", 				default="TF ChIP-Seq")
	parser.add_argument('--applets', help="Name of project containing applets", default="E3 ChIP-seq")

	args = parser.parse_args()

	global DEBUG
	DEBUG = args.debug
	if DEBUG:
		logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
	else: #use the defaulf logging level
		logging.basicConfig(format='%(levelname)s:%(message)s')

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

	m = re.match(r'''^([A-Za-z0-9_\-\ ]+):(\S+)''', identifier)
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

	m = re.match(r'''(^\S+)/(\S+)''', file_identifier)
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
				logging.warning('Could not find file %s. Skipping' %(identifier))
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

def main():
	args = get_args()

	output_project = resolve_project(args.outp, 'w')
	logging.info('Found output project %s' %(output_project.name))
	output_folder = resolve_folder(output_project, args.outf)
	logging.info('Using output folder %s' %(output_folder))
	applet_project = resolve_project(args.applets, 'r')
	logging.info('Found applet project %s' %(applet_project.name))

	rep1_fastq = resolve_file(args.rep1)
	rep2_fastq = resolve_file(args.rep2)
	ctl1_fastq = resolve_file(args.ctl1)
	ctl2_fastq = resolve_file(args.ctl2)

	if not (args.rep1 or args.rep2 or args.ctl1 or args.ctl2):
		wf = blank_workflow(args)
		return
	else:
		wf = None

	'''
	stages = {}

	for infile in [args.rep1, args.rep2, args.ctl1, args.ctl2]:
		stages.update(map_and_filter(infile,args))

	for expvsctl in [(args.rep1, args.ctl1), (args.rep2, args.ctl2)]:
		stages.update(call_peaks(expvsctl, args))
	'''

if __name__ == '__main__':
	main()
