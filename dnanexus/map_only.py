#!/usr/bin/env python

import pdb
import os.path, sys, subprocess, logging, re, json, urlparse, requests
import dxpy

EPILOG = '''Notes:

Examples:

	%(prog)s
'''

DEFAULT_APPLET_PROJECT = 'E3 ChIP-seq'
MAPPING_APPLET_NAME = 'encode_bwa'
POOL_APPLET_NAME = 'pool'
XX_REFERENCE = 'E3 ChIP-seq:/GRCh38/GRCh38_minimal_X.tar.gz'
XY_REFERENCE = 'E3 ChIP-seq:/GRCh38/GRCh38_minimal_XY.tar.gz'
KEYFILE = os.path.expanduser("~/keypairs.json")
DEFAULT_SERVER = 'https://www.encodeproject.org'

APPLETS = {}

def get_args():
	import argparse
	parser = argparse.ArgumentParser(
		description=__doc__, epilog=EPILOG,
		formatter_class=argparse.RawDescriptionHelpFormatter)

	parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
	parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
	parser.add_argument('--debug',   help="Print debug messages", 				default=False, action='store_true')
	parser.add_argument('--outp',    help="Output project name or ID", 			default=dxpy.WORKSPACE_ID)
	parser.add_argument('--outf',    help="Output folder name or ID", 			default="/")
	parser.add_argument('--applets', help="Name of project containing applets", default=DEFAULT_APPLET_PROJECT)
	parser.add_argument('--instance_type', help="Instance type for applets",	default=None)
	parser.add_argument('--key', help="The keypair identifier from the keyfile.  Default is --key=default",
		default='default')

	args = parser.parse_args()

	global DEBUG
	DEBUG = args.debug
	if DEBUG:
		logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
	else: #use the defaulf logging level
		logging.basicConfig(format='%(levelname)s:%(message)s')

	return args

def processkey(key):

    if key:
        keysf = open(KEYFILE,'r')
        keys_json_string = keysf.read()
        keysf.close()
        keys = json.loads(keys_json_string)
        key_dict = keys[key]
    else:
        key_dict = {}
    AUTHID = key_dict.get('key')
    AUTHPW = key_dict.get('secret')
    if key:
        SERVER = key_dict.get('server')
    else:
        SERVER = DEFAULT_SERVER

    if not SERVER.endswith("/"):
        SERVER += "/"

    return (AUTHID,AUTHPW,SERVER)

def encoded_get(url, keypair=None):
    HEADERS = {'content-type': 'application/json'}
    url = urlparse.urljoin(url,'?format=json&frame=embedded')
    if keypair:
        response = requests.get(url, auth=keypair, headers=HEADERS)
    else:
        response = requests.get(url, headers=HEADERS)
    return response.json()

def blank_workflow(args):
	return

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

def filenames_in(files=None):
	if not len(files):
		return []
	else:
		return [f.get('submitted_file_name') for f in files]

def files_to_map(exp_obj):
	if not exp_obj:
		return []
	else:
		files = []
		for file_obj in exp_obj.get('files'):
			if file_obj.get('output_type') == 'reads' and \
			   file_obj.get('file_format') == 'fastq' and \
			   file_obj.get('submitted_file_name') not in filenames_in(files):
			   files.extend([file_obj])
			elif file_obj.get('submitted_file_name') in filenames_in(files):
				logging.warning('%s:%s Duplicate filename, ignoring.' %(exp_obj.get('accession'),file_obj.get('accession')))
		return files

def replicates_to_map(files):
	if not files:
		return []
	else:
		return [f.get('replicate') for f in files]

def map_paired(experiment, biorep_n, files=(None,None)):
	pass

def map_unpaired(experiment, biorep_n, files=None):
	pass

def main():
	args = get_args()
	authid, authpw, server = processkey(args.key)
	keypair = (authid,authpw)

	for exp_id in args.infile:
		outstrings = []
		encode_url = urlparse.urljoin(server,exp_id.rstrip())
		experiment = encoded_get(encode_url, keypair)
		outstrings.extend([exp_id.rstrip()])
		outstrings.extend([experiment.get('accession')])
		files = files_to_map(experiment)
		replicates = replicates_to_map(files)
		if files:
			for biorep_n in set([rep.get('biological_replicate_number') for rep in replicates]):
				biorep_files = [f for f in files if f.get('replicate').get('biological_replicate_number') == biorep_n]
				paired_files = []
				unpaired_files = []
				while biorep_files:
					file_object = biorep_files.pop()
					if file_object.get('paired_end') == None: # group all the unpaired reads for this biorep together
						unpaired_files.extend([file_object])
					elif file_object.get('paired_end') in ['1','2']:
						if file_object.get('paired_with'):
							mate = next((f for f in biorep_files if f.get('@id') == file_object.get('paired_with')), None)
							biorep_files.remove(mate)
						else: #have to find the file that is paired with this one
							mate = next((f for f in biorep_files if f.get('paired_with') == file_object.get('@id')), None)
							biorep_files.remove(mate)
						paired_files.extend([(file_object,mate)])
				if paired_files:
					logging.debug('paired_files: %s' %([(a.get('accession'), b.get('accession')) for (a,b) in paired_files]))
				else:
					logging.debug('paired_files: %s' %(None))
				if unpaired_files:
					logging.debug('unpaired_files: %s' %([f.get('accession') for f in unpaired_files]))
				else:
					logging.debug('unpaired_files: %s' %(None))
				if biorep_files:
					logging.warning('%s: leftover file(s) %s' %(experiment.get('accession'), biorep_files))
				#map_paired(experiment, biorep_n, paired_files)
				#map_unpaired(experiment, biorep_n, unpaired_files)
		else:
			logging.warning('%s: No files to map' %experiment.get('accession'))


		outstrings.extend([str(len(files))])
		outstrings.extend([str([f.get('accession') for f in files])])
		print '\t'.join(outstrings)
if __name__ == '__main__':
	main()
