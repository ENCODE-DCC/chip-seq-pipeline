#!/usr/bin/env python

import os, sys, subprocess, logging, dxpy

EPILOG = '''Notes:

Examples:

	%(prog)s
'''

DEFAULT_APPLET_PROJECT = 'E3 ChIP-seq'
KEYFILE = os.path.expanduser("~/keypairs.json")

def get_args():
	import argparse
	parser = argparse.ArgumentParser(
		description=__doc__, epilog=EPILOG,
		formatter_class=argparse.RawDescriptionHelpFormatter)

	parser.add_argument('infile', help="Experiment accessions", nargs='?', type=argparse.FileType('r'), default=sys.stdin)
	parser.add_argument('--debug',   help="Print debug messages", 				default=False, action='store_true')
	parser.add_argument('--project',    help="Project name or ID", 			default=dxpy.WORKSPACE_ID)
	parser.add_argument('--outf',    help="Output folder name or ID", 			default="/")
	parser.add_argument('--inf', nargs='*',    help="Folder(s) name or ID with tagAligns", 			default="/")
	parser.add_argument('--yes',   help="Run the workflows created", 			default=False, action='store_true')
	parser.add_argument('--tag',   help="String to add to the workflow name")
	parser.add_argument('--key', help="The keypair identifier from the keyfile.  Default is --key=default", default='default')
	parser.add_argument('--gsize', help="Genome size string for MACS2, e.g. mm or hs", required=True)
	parser.add_argument('--csizes', help="chrom.sizes file for bedtobigbed, e.g. ENCODE Reference Files:/mm10/male.mm10.chrom.sizes", required=True)
	parser.add_argument('--idr', help="Run IDR", default=False, action='store_true')
	parser.add_argument('--idrversion', help="IDR version (relevant only if --idr is specified", default="2")

	args = parser.parse_args()

	if args.debug:
		logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
	else: #use the defaulf logging level
		logging.basicConfig(format='%(levelname)s:%(message)s')

	return args

def processkey(key):

	import json

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
	import urlparse, requests
	HEADERS = {'content-type': 'application/json'}
	url = urlparse.urljoin(url,'?format=json&frame=embedded&datastore=database')
	if keypair:
		response = requests.get(url, auth=keypair, headers=HEADERS)
	else:
		response = requests.get(url, headers=HEADERS)
	return response.json()

def get_control_id(experiment):
	# url = server + '/experiments/%s/' %(exp_id)
	# experiment = encoded_get(url, keypair)
	possible_controls = experiment.get('possible_controls')
	if not possible_controls or len(possible_controls) != 1:
		logging.error("Tried to find one possible control, found %s" %(possible_controls))
		return None
	return possible_controls[0].get('accession')

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

def get_tas(exp_id, default_project, ta_folders):
	possible_files = []
	for base_folder in ta_folders:
		if ':' in base_folder:
			project_name, path = base_folder.split(':')
			project = resolve_project(project_name)
			project = project.get_id()
			project_name += ":"
		else:
			project = default_project
			project_name = ""
			path = base_folder
		if not path.startswith('/'):
			path = '/' + path
		print project, project_name, path
		for dxfile in dxpy.find_data_objects(classname='file', state='closed', folder=path, describe=True, recurse=True, project=project):
			desc = dxfile.get('describe')
			if exp_id in desc.get('folder') and '/bams' in desc.get('folder') and desc.get('name').endswith(('tagAlign', 'tagAlign.gz')):
				possible_files.append(desc)
	print "%s %i possible files" %(exp_id, len(possible_files))
	rep1_files = [f for f in possible_files if 'rep1' in f.get('folder')]
	rep2_files = [f for f in possible_files if 'rep2' in f.get('folder')]
	if len(rep1_files) != 1:
		print "Tried to find one rep1 ta, found %d" %(len(rep1_files))
		rep1 = None
	else:
		rep1 = rep1_files[0].get('project') + ':' + rep1_files[0].get('folder') + '/' + rep1_files[0].get('name')
	if len(rep2_files) != 1:
		print "Tried to find one rep2 ta, found %d" %(len(rep2_files))
		rep2 = None
	else:
		rep2 = rep2_files[0].get('project') + ':' + rep2_files[0].get('folder') + '/' + rep2_files[0].get('name')
	
	return rep1, rep2

def get_exp_tas(exp_id, default_project, ta_folders):
	return get_tas(exp_id, default_project, ta_folders)

def get_ctl_tas(exp_id, default_project, ta_folders):
	rep1, rep2 = get_tas(exp_id, default_project, ta_folders)
	if rep1 and rep2:
		return rep1,rep2
	elif rep1:
		return rep1,rep1
	elif rep2:
		return rep2,rep2
	else: #both are None - just pass it on and let the calling code deal with it
		return rep1,rep2

def main():
	args = get_args()
	authid, authpw, server = processkey(args.key)
	keypair = (authid,authpw)

	for exp_id in args.infile:
		exp_id = exp_id.rstrip()
		print "Experiment %s" %(exp_id)
		url = server + '/experiments/%s/' %(exp_id)
		experiment = encoded_get(url, keypair)
		print "%s %s %s" %(experiment['accession'], experiment.get('target')['investigated_as'], experiment.get('description'))
		ctl_id = get_control_id(experiment)
		if ctl_id:
			print "Control %s" %(ctl_id)
		else:
			print "Found no control ... skipping %s" %(exp_id)
			continue
		rep1_ta, rep2_ta = get_exp_tas(exp_id, args.project, args.inf)
		ctl1_ta, ctl2_ta = get_ctl_tas(ctl_id, args.project, args.inf)
		if not (rep1_ta and rep2_ta and ctl1_ta and ctl2_ta):
			print "Skipping %s" %(exp_id)
			continue
		workflow_name = '%s Peaks' %(exp_id)
		if args.tag:
			workflow_name += ' %s' %(args.tag)
		outf = args.outf
		if not outf.startswith('/'):
			outf = '/'+outf
		outf += '/%s/peaks/' %(exp_id)
		try:
			investigated_as = experiment.get('target')['investigated_as']
			print investigated_as
		except:
			print "Failed to determine target type ... skipping %s" %(exp_id)
			continue
		if any('histone' in target_type for target_type in investigated_as):
			print "Found to be histone"
			workflow_spinner = '~/chip-seq-pipeline/dnanexus/histone_workflow.py'
		else:
			print "Assumed to be tf"
			workflow_spinner = '~/chip-seq-pipeline/dnanexus/tf_workflow.py'
		run_command = \
			'%s --debug --name "%s" --outf "%s" --nomap --yes ' %(workflow_spinner, workflow_name, outf) + \
			'--rep1pe false --rep2pe false ' + \
			'--rep1 %s --rep2 %s ' %(rep1_ta, rep2_ta) + \
			'--ctl1 %s --ctl2 %s ' %(ctl1_ta, ctl2_ta) + \
			'--genomesize %s --chrom_sizes "%s"' %(args.gsize, args.csizes)
		if args.idr:
			run_command += ' --idr --idrversion %s' %(args.idrversion)

		print run_command
		try:
			subprocess.check_call(run_command, shell=True)
		except subprocess.CalledProcessError as e:
			logging.error("%s exited with non-zero code %d" %(workflow_spinner, e.returncode))

if __name__ == '__main__':
	main()
