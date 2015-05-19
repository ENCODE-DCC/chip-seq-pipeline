#!/usr/bin/env python

import os, sys, subprocess, logging, dxpy, json, re, socket, getpass, urlparse, datetime, requests, time
import common
import dateutil.parser

logger = logging.getLogger(__name__)

EPILOG = '''Notes:

Examples:

	%(prog)s
'''

DEFAULT_APPLET_PROJECT = 'E3 ChIP-seq'

def get_args():
	import argparse
	parser = argparse.ArgumentParser(
		description=__doc__, epilog=EPILOG,
		formatter_class=argparse.RawDescriptionHelpFormatter)

	parser.add_argument('analysis_ids',	nargs='*', default=None)
	parser.add_argument('--infile',		help='File containing analysis IDs to accession', type=argparse.FileType('r'), default=sys.stdin)
	parser.add_argument('--outfile', 	help='tsv table of files with metadata', type=argparse.FileType('w'), default=sys.stdout)
	parser.add_argument('--assembly', 	help='Genome assembly like hg19 or mm10', required=True)
	parser.add_argument('--debug',		help="Print debug messages", default=False, action='store_true')
	parser.add_argument('--project',	help="Project name or ID", default=dxpy.WORKSPACE_ID)
	parser.add_argument('--key',		help="The keypair identifier from the keyfile.", default='www')
	parser.add_argument('--keyfile',	help="The keyfile.", default=os.path.expanduser("~/keypairs.json"))
	parser.add_argument('--tag',		help="A short string to add to the composite track longLabel")
	parser.add_argument('--dryrun',		help="Don't POST or upload", default=False, action='store_true')
	parser.add_argument('--force',		help="Try to POST and upload even if the file appears to be a duplicate", default=False, action='store_true')

	args = parser.parse_args()

	if args.debug:
		logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
	else: #use the defaulf logging level
		logging.basicConfig(format='%(levelname)s:%(message)s')

	return args

def after(date1, date2):
	return(dateutil.parser.parse(date1) > dateutil.parser.parse(date2))

def get_rep_bams(experiment, keypair, server):

	original_files = [common.encoded_get(urlparse.urljoin(server,'%s' %(uri)), keypair) for uri in experiment.get('original_files')]

	#resolve the biorep_n for each fastq
	for fastq in [f for f in original_files if f.get('file_format') == 'fastq']:
		replicate = common.encoded_get(urlparse.urljoin(server,'%s' %(fastq.get('replicate'))), keypair)
		fastq.update({'biorep_n' : replicate.get('biological_replicate_number')})
	#resolve the biorep_n's from derived_from for each bam
	for bam in [f for f in original_files if f.get('file_format') == 'bam']:
		biorep_ns = set()
		for derived_from_uri in bam.get('derived_from'):
			derived_from_accession = os.path.basename(derived_from_uri.strip('/')) #this assumes frame=object
			biorep_ns.add(next(f.get('biorep_n') for f in original_files if f.get('accession') == derived_from_accession))
		if len(biorep_ns) != 1:
			print >> sys.stderr, "%s %s expected 1 biorep_n, found %d, skipping." %(experiment_accession, bam.get('accession'))
			return
		else:
			biorep_n = biorep_ns.pop()
			bam.update({'biorep_n': biorep_n})
	#remove any bams that are older than another bam (resultsing in only the most recent surviving)
	for bam in [f for f in original_files if f.get('file_format') == 'bam' and f.get('biorep_n') == biorep_n and after(bam.get('date_created'), f.get('date_created'))]:
		original_files.remove(bam)

	rep1_bam = next(f for f in original_files if f.get('file_format') == 'bam' and f.get('biorep_n') == 1)
	rep2_bam = next(f for f in original_files if f.get('file_format') == 'bam' and f.get('biorep_n') == 2)

	return rep1_bam, rep2_bam

def accession_file(f, keypair, server, dryrun, force):
	#check for duplication
	#download
	#calculate md5 and add to f.md5sum
	#post file and get accession, upload credentials
	#upload to S3
	#remove the local file (to save space)
	#return the ENCODEd file object
	already_accessioned = False
	dx = f.pop('dx')
	for tag in dx.tags:
		m = re.search(r'(ENCFF\d{3}\D{3})|(TSTFF\D{6})', tag)
		if m:
			logger.info('%s appears to contain ENCODE accession number in tag %s ... skipping' %(dx.get_id(),m.group(0)))
			already_accessioned = True
			break
	if already_accessioned and not force:
		return
	url = urlparse.urljoin(server, 'search/?type=file&submitted_file_name=%s&format=json&frame=object' %(f.get('submitted_file_name')))
	r = requests.get(url,auth=keypair)
	try:
		r.raise_for_status()
		if r.json()['@graph']:
			for duplicate_item in r.json()['@graph']:
				if duplicate_item.get('status')  == 'deleted':
					logger.info("A potential duplicate file was found but its status=deleted ... proceeding")
					duplicate_found = False
				else:
					logger.info("Found potential duplicate: %s" %(duplicate_item.get('accession')))
					if submitted_file_size ==  duplicate_item.get('file_size'):
						logger.info("%s %s: File sizes match, assuming duplicate." %(str(submitted_file_size), duplicate_item.get('file_size')))
						duplicate_found = True
						break
					else:
						logger.info("%s %s: File sizes differ, assuming new file." %(str(submitted_file_size), duplicate_item.get('file_size')))
						duplicate_found = False
		else:
			duplicate_found = False
	except:
		logger.warning('Duplicate accession check failed: %s %s' % (r.status_code, r.reason))
		logger.debug(r.text)
		duplicate_found = False

	if duplicate_found:
		if force:
			logger.info("Duplicate detected, but force=true, so continuing")
		else:
			logger.info("Duplicate detected, skipping")
			return

	local_fname = dx.name
	logger.info("Downloading %s" %(local_fname))
	dxpy.download_dxfile(dx.get_id(),local_fname)
	f.update({'md5sum': common.md5(local_fname)})
	f['notes'] = json.dumps(f.get('notes'))

	url = urlparse.urljoin(server,'files/')
	if dryrun:
		logger.info("Dry run.  Would POST %s" %(f))
		new_file_object = {}
	else:
		r = requests.post(url, auth=keypair, headers={'content-type': 'application/json'}, data=json.dumps(f))
		try:
			r.raise_for_status()
			new_file_object = r.json()['@graph'][0]
			logger.info("New accession: %s" %(new_file_object.get('accession')))
		except:
			logger.warning('POST file object failed: %s %s' % (r.status_code, r.reason))
			logger.debug(r.text)
			new_file_object = {}
			if r.status_code == 409:
				try: #cautiously add a tag with the existing accession number
					if calculated_md5 in r.json().get('detail'):
						url = urlparse.urljoin(server,'/search/?type=file&md5sum=%s' %(calculated_md5))
						r = requests.get(url,auth=keypair)
						r.raise_for_status()
						accessioned_file = r.json()['@graph'][0]
						existing_accession = accessioned_file['accession']
						dx.add_tags([existing_accession])
						logger.info('Already accessioned.  Added %s to dxfile tags' %(existing_accession))
				except:
					logger.info('Conflict does not appear to be md5 ... continuing')
		if new_file_object:
			creds = new_file_object['upload_credentials']
			env = os.environ.copy()
			env.update({
				'AWS_ACCESS_KEY_ID': creds['access_key'],
				'AWS_SECRET_ACCESS_KEY': creds['secret_key'],
				'AWS_SECURITY_TOKEN': creds['session_token'],
			})

			logger.info("Uploading file.")
			start = time.time()
			try:
				subprocess.check_call(['aws', 's3', 'cp', local_fname, creds['upload_url'], '--quiet'], env=env)
			except subprocess.CalledProcessError as e:
				# The aws command returns a non-zero exit code on error.
				logger.error("Upload failed with exit code %d" % e.returncode)
				upload_returncode = e.returncode
			else:
				upload_returncode = 0
				end = time.time()
				duration = end - start
				logger.info("Uploaded in %.2f seconds" % duration)
				dx.add_tags([new_file_object.get('accession')])
		else:
			upload_returncode = -1

	try:
		os.remove(local_fname)
	except:
		pass

	return common.encoded_get(urlparse.urljoin(server,'/files/%s' %(new_file_object.get('accession')), keypair))

def accession_analysis(analysis_id, keypair, server, assembly, dryrun, force):
	analysis_id = analysis_id.strip()
	analysis = dxpy.describe(analysis_id)
	project = analysis.get('project')

	m = re.match('^(ENCSR[0-9]{3}[A-Z]{3}) Peaks',analysis['executableName'])
	if m:
		experiment_accession = m.group(1)
	else:
		logger.info("No accession in %s, skipping." %(analysis['executableName']))
		return

	experiment = common.encoded_get(urlparse.urljoin(server,'/experiments/%s' %(experiment_accession)), keypair)
	bams = get_rep_bams(experiment, keypair, server)
	rep1_bam = bams[0]['accession']
	rep2_bam = bams[1]['accession']

	common_metadata = {
		'assembly': assembly,
		'lab': 'encode-processing-pipeline',
		'award': 'U41HG006992',
		}

	narrowpeak_metadata = common.merge_dicts(
		{'file_format': 'bed_narrowPeak', 'file_format_specifications': ['ENCODE:narrowPeak.as'], 'output_type': 'peaks'}, common_metadata)
	replicated_narrowpeak_metadata = common.merge_dicts(
		{'file_format': 'bed_narrowPeak', 'file_format_specifications': ['ENCODE:narrowPeak.as'], 'output_type': 'replicated peaks'}, common_metadata)

	gappedpeak_metadata = common.merge_dicts(
		{'file_format': 'bed_gappedPeak', 'file_format_specifications': ['ENCODE:gappedPeak.as'], 'output_type': 'peaks'}, common_metadata)
	replicated_gappedpeak_metadata = common.merge_dicts(
		{'file_format': 'bed_gappedPeak', 'file_format_specifications': ['ENCODE:gappedPeak.as'], 'output_type': 'replicated peaks'}, common_metadata)

	narrowpeak_bb_metadata = common.merge_dicts(
		{'file_format': 'narrowPeak', 'file_format_specifications': ['ENCODE:narrowPeak.as'], 'output_type': 'peaks'}, common_metadata)
	replicated_narrowpeak_bb_metadata = common.merge_dicts(
		{'file_format': 'narrowPeak', 'file_format_specifications': ['ENCODE:narrowPeak.as'], 'output_type': 'replicated peaks'}, common_metadata)

	gappedpeak_bb_metadata = common.merge_dicts(
		{'file_format': 'gappedPeak', 'file_format_specifications': ['ENCODE:gappedPeak.as'], 'output_type': 'peaks'}, common_metadata)
	replicated_gappedpeak_bb_metadata = common.merge_dicts(
		{'file_format': 'gappedPeak', 'file_format_specifications': ['ENCODE:gappedPeak.as'], 'output_type': 'replicated peaks'}, common_metadata)

	fc_signal_metadata = 	common.merge_dicts(
		{'file_format': 'bigWig', 'output_type': 'fold change over control'}, common_metadata)
	pvalue_signal_metadata = common.merge_dicts(
		{'file_format': 'bigWig', 'output_type': 'signal p-value'}, common_metadata)

	stage_outputs = {
		"ENCODE Peaks" : {
			'files': [
				common.merge_dicts({'name': 'rep1_narrowpeaks', 		'derived_from': [rep1_bam]},			narrowpeak_metadata),
				common.merge_dicts({'name': 'rep2_narrowpeaks', 		'derived_from': [rep2_bam]},			narrowpeak_metadata),
				common.merge_dicts({'name': 'pooled_narrowpeaks',		'derived_from': [rep1_bam, rep2_bam]},	narrowpeak_metadata),
				common.merge_dicts({'name': 'rep1_narrowpeaks_bb', 		'derived_from': [rep1_bam]},			narrowpeak_bb_metadata),
				common.merge_dicts({'name': 'rep2_narrowpeaks_bb', 		'derived_from': [rep2_bam]},			narrowpeak_bb_metadata),
				common.merge_dicts({'name': 'pooled_narrowpeaks_bb',	'derived_from': [rep1_bam, rep2_bam]},	narrowpeak_bb_metadata),
				common.merge_dicts({'name': 'rep1_gappedpeaks', 		'derived_from': [rep1_bam]},			gappedpeak_metadata),
				common.merge_dicts({'name': 'rep2_gappedpeaks', 		'derived_from': [rep2_bam]},			gappedpeak_metadata),
				common.merge_dicts({'name': 'pooled_gappedpeaks', 		'derived_from': [rep1_bam, rep2_bam]},	gappedpeak_metadata),
				common.merge_dicts({'name': 'rep1_gappedpeaks_bb', 		'derived_from': [rep1_bam]},			gappedpeak_bb_metadata),
				common.merge_dicts({'name': 'rep2_gappedpeaks_bb', 		'derived_from': [rep2_bam]},			gappedpeak_bb_metadata),
				common.merge_dicts({'name': 'pooled_gappedpeaks_bb', 	'derived_from': [rep1_bam, rep2_bam]},	gappedpeak_bb_metadata),
				common.merge_dicts({'name': 'rep1_pvalue_signal',		'derived_from': [rep1_bam]},			pvalue_signal_metadata),
				common.merge_dicts({'name': 'rep2_pvalue_signal',		'derived_from': [rep2_bam]},			pvalue_signal_metadata),
				common.merge_dicts({'name': 'pooled_pvalue_signal',		'derived_from': [rep1_bam, rep2_bam]},	pvalue_signal_metadata),
				common.merge_dicts({'name': 'rep1_fc_signal',			'derived_from': [rep1_bam]},			fc_signal_metadata),
				common.merge_dicts({'name': 'rep2_fc_signal',			'derived_from': [rep2_bam]},			fc_signal_metadata),
				common.merge_dicts({'name': 'pooled_fc_signal',			'derived_from': [rep1_bam, rep2_bam]},	fc_signal_metadata)],
			'qc': []},
		"Overlap narrowpeaks": {
			'files': [
				common.merge_dicts({'name': 'overlapping_peaks',		'derived_from': [rep1_bam, rep2_bam]},	replicated_narrowpeak_metadata),
				common.merge_dicts({'name': 'overlapping_peaks_bb',		'derived_from': [rep1_bam, rep2_bam]},	replicated_narrowpeak_bb_metadata)],
			'qc': ['npeaks_in', 'npeaks_out', 'npeaks_rejected']},
		"Overlap gappedpeaks": {
			'files': [
				common.merge_dicts({'name': 'overlapping_peaks',		'derived_from': [rep1_bam, rep2_bam]},	replicated_gappedpeak_metadata),
				common.merge_dicts({'name': 'overlapping_peaks_bb',		'derived_from': [rep1_bam, rep2_bam]},	replicated_gappedpeak_bb_metadata)],
			'qc': ['npeaks_in', 'npeaks_out', 'npeaks_rejected']}
		}

	experiment = common.encoded_get(urlparse.urljoin(server,'/experiments/%s' %(experiment_accession)), keypair)
	rep1_bam, rep2_bam = get_rep_bams(experiment, keypair, server)

	files = []
	for (stage_name, outputs) in stage_outputs.iteritems():
		stage_metadata = next(s['execution'] for s in analysis.get('stages') if s['execution']['name'] == stage_name)
		for static_metadata in outputs['files']:
			output_name = static_metadata['name']
			dx = dxpy.DXFile(stage_metadata['output'][output_name], project=project)
			file_metadata = {
				'dx': dx,
				'notes': {
					'dx-id': dx.get_id(),
					'dx-createdBy': {
						'job': stage_metadata['id'],
						'executable': stage_metadata['executable'], #todo get applet ID
						'user': stage_metadata['launchedBy']},
					'qc': dict(zip(outputs['qc'],[stage_metadata['output'][metric] for metric in outputs['qc']]))}, #'aliases': ['ENCODE:%s-%s' %(experiment.get('accession'), static_metadata.pop('name'))],
				'dataset': experiment.get('accession'),
				'file_size': dx.describe().get('size'),
				'submitted_file_name': dx.get_proj_id() + ':' + '/'.join([dx.folder,dx.name])}
			static_metadata.pop('name')
			file_metadata.update(static_metadata)
			files.append(file_metadata)

	for f in files:
		f.update({'new_file_object' : accession_file(f, keypair, server, dryrun, force)})

	return files

def main():

	args = get_args()
	if args.debug:
		logger.setLevel(logging.DEBUG)
	else:
		logger.setLevel(logging.INFO)

	authid, authpw, server = common.processkey(args.key, args.keyfile)
	keypair = (authid,authpw)

	if args.analysis_ids:
		ids = args.analysis_ids
	else:
		ids = args.infile

	for (i, analysis_id) in enumerate(ids):
		logger.info('%s' %(analysis_id))
		accessioned_files = accession_analysis(analysis_id, keypair, server, args.assembly, args.dryrun, args.force)

if __name__ == '__main__':
	main()
