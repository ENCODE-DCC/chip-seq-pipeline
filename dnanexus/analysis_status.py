#!/usr/bin/env python

import logging, os, sys, urlparse, requests, json
import common

EPILOG = '''Notes:

Examples:

	%(prog)s
'''

def get_args():
	import argparse
	parser = argparse.ArgumentParser(
		description=__doc__, epilog=EPILOG,
		formatter_class=argparse.RawDescriptionHelpFormatter)

	parser.add_argument('experiments', help='Experiment accessions', nargs='*')
	parser.add_argument('--infile', help="Infile with experiment accessions", nargs='?', type=argparse.FileType('r'), default=None)
	parser.add_argument('--status', help='The status to patch', default='replaced')
	parser.add_argument('--debug',   help="Print debug messages", 				default=False, action='store_true')
	parser.add_argument('--key', help="The keypair identifier from the keyfile.  Default is --key=default", default='default')
	parser.add_argument('--keyfile', help="The JSON-formatted file of keypairs.  Default is ~/keypairs.json", default=os.path.expanduser("~/keypairs.json"))
	parser.add_argument('--dryrun', help="Show the patch payload, but don't change anything", default=False, action='store_true')
	parser.add_argument('--force', help="Patch even released files", default=False, action='store_true')

	args = parser.parse_args()

	if args.debug:
		logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
	else: #use the defaulf logging level
		logging.basicConfig(format='%(levelname)s:%(message)s')

	return args


def main():
	args = get_args()
	authid, authpw, server = common.processkey(args.key, args.keyfile)
	keypair = (authid,authpw)

	if args.infile and args.experiments:
		experiments = args.experiments
		experiments.extend([e.strip() for e in args.infile if e.strip()])
	elif args.infile:
		experiments = args.infile
	else:
		experiments = args.experiments

	for exp_id in experiments:
		uri = '/experiments/%s' %(exp_id)
		experiment = common.encoded_get(urlparse.urljoin(server,'%s' %(uri)), keypair)
		if experiment.get('status') == 'error':
			print experiment
			print "Error fetching %s ... skipping" %(exp_id)
			continue

		print experiment.get('accession')
		for uri in experiment['original_files']:
			url = urlparse.urljoin(server,'%s' %(uri))
			file_obj = common.encoded_get(url, keypair)
			print "%s, %s, %s, %s, %s, %s" %(file_obj.get('accession'),file_obj.get('file_type'),file_obj.get('file_format'),file_obj.get('file_format_type'),file_obj.get('output_type'),file_obj.get('status'))
			if file_obj.get('file_format') in ['bed', 'bigBed', 'bigWig']:
				if file_obj.get('status') != 'released' or args.force:
					patch_payload = {'status': args.status}
					if args.dryrun:
						print "--dryrun:  would have patched %s" %(json.dumps(patch_payload))
					else:
						r = requests.patch(url, auth=keypair, data=json.dumps(patch_payload), headers={'content-type': 'application/json', 'accept': 'application/json'})
						try:
							r.raise_for_status()
						except:
							print(r.text)
							print('Patch failed: %s %s ... skipping' % (r.status_code, r.reason))
							continue
						else:
							print "Patched %s" %(json.dumps(patch_payload))

if __name__ == '__main__':
	main()
