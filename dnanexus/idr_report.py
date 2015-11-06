#!/usr/bin/env python

import os.path, sys, subprocess, logging, re, json, urlparse, requests, csv, time, pprint
import common
import dxpy

logger = logging.getLogger(__name__)

EPILOG = '''Notes:

Examples:

	%(prog)s
'''

class InputError(Exception):
	pass

def get_args():
	import argparse
	parser = argparse.ArgumentParser(
		description=__doc__, epilog=EPILOG,
		formatter_class=argparse.RawDescriptionHelpFormatter)

	parser.add_argument('analysis_ids',	help='List of analysis IDs to report on', nargs='*', default=None)
	parser.add_argument('--infile',		help='File containing analysis IDs', type=argparse.FileType('r'), default=sys.stdin)
	parser.add_argument('--outfile', 	help='csv output', type=argparse.FileType('wb'), default=sys.stdout)
	parser.add_argument('--assembly', 	help='Genome assembly to report on (like hg19 or mm10)', required=True)
	parser.add_argument('--debug',		help="Print debug messages", default=False, action='store_true')
	parser.add_argument('--key',		help="The keypair identifier from the keyfile.", default='www')
	parser.add_argument('--keyfile',	help="The keyfile.", default=os.path.expanduser("~/keypairs.json"))
	parser.add_argument('--created_after', help="String to search for analyses (instead of looking in --infile or arguments) in the DNAnexus form like -5d", default=None)

	args = parser.parse_args()

	if args.debug:
		logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
	else: #use the defaulf logging level
		logging.basicConfig(format='%(levelname)s:%(message)s')

	return args

def get_experiment_accession(analysis):
	m_executableName = re.search('(ENCSR[0-9]{3}[A-Z]{3})',analysis['executableName'])
	m_name = re.search('(ENCSR[0-9]{3}[A-Z]{3})',analysis['name'])
	if not (m_executableName or m_name):
		logger.error("No experiment accession in name %s or executableName %s." %(analysis['name'], analysis['executableName']))
		return
	elif (m_executableName and m_name):
		executableName_accession = m_executableName.group(1)
		name_accession = m_name.group(1)
		if executableName_accession == name_accession:
			return executableName_accession
		else:
			logger.error('Different experiment accessions in name %s and executableName %s.' %(analysis['name'], analysis['executableName']))
			return None
	else:
		m = (m_executableName or m_name)
		experiment_accession = m.group(1)
		logger.debug("get_experiment_accession returning %s" %(experiment_accession))
		return experiment_accession


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
	elif args.created_after:
		analyses = dxpy.find_analyses(name="ENCSR*",name_mode='glob',state='done',include_subjobs=True,return_handler=True,created_after="%s" %(args.created_after))
		ids = [analysis.get_id() for analysis in analyses if analysis.describe()['executableName'] == 'tf_chip_seq']
	elif args.infile:
		ids = args.infile
	else:
		#never reached because inile defaults to stdin
		raise InputError("Must supply analysis id's in arguments, --infile or supply search string in --created_after")

	fieldnames = [	'date','analysis','experiment','target','biosample_term_name','biosample_type','lab','rfa','assembly',
					'Nt','Np','N1','N2','rescue_ratio','self_consistency_ratio','reproducibility_test',
					'state','total price','notes']
	writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, delimiter='\t', quotechar='"')
	writer.writeheader()

	for (i, analysis_id) in enumerate(ids):
		if analysis_id.startswith('#'):
			continue
		analysis_id = analysis_id.rstrip()
		logger.debug('%s' %(analysis_id))
		analysis = dxpy.DXAnalysis(analysis_id)
		desc = analysis.describe()
		project = desc.get('project')

		m = re.match('^(ENCSR[0-9]{3}[A-Z]{3}) Peaks',desc['name'])
		if m:
			experiment_accession = m.group(1)
		else:
			logger.error("No accession in %s, skipping." %(desc['name']))
			continue

		experiment = common.encoded_get(urlparse.urljoin(server,'/experiments/%s' %(experiment_accession)), keypair)
		logger.debug('ENCODEd experiment %s' %(experiment['accession']))
		try:
			idr_stage = next(s['execution'] for s in desc['stages'] if s['execution']['name'] == "Final IDR peak calls")
		except:
			logging.error('Failed to find IDR stage in %s' %(analysis_id))
		else:
			Np = idr_stage['output'].get('Np')
			N1 = idr_stage['output'].get('N1')
			N2 = idr_stage['output'].get('N2')
			Nt = idr_stage['output'].get('Nt')
			rescue_ratio = idr_stage['output'].get('rescue_ratio')
			self_consistency_ratio = idr_stage['output'].get('self_consistency_ratio')
			reproducibility_test = idr_stage['output'].get('reproducibility_test')

		done_time = next(transition['setAt'] for transition in desc['stateTransitions'] if transition['newState'] == "done")

		row = {
			'date': time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(done_time/1000)),
			'analysis':		analysis.get_id(),
			'experiment': 	experiment.get('accession'),
			'target':		experiment['target'].split('/')[2],
			'biosample_term_name':	experiment.get('biosample_term_name'),
			'biosample_type':	experiment.get('biosample_type'),
			'lab':			experiment['lab'].split('/')[2],
			'rfa':			common.encoded_get(server+experiment.get('award'),keypair).get('rfa'),
			'assembly':		args.assembly, #TODO ... derive this from the analysis
			'Np':			Np,
			'N1':			N1,
			'N2':			N2,
			'Nt':			Nt,
			'rescue_ratio':	rescue_ratio,
			'self_consistency_ratio': self_consistency_ratio,
			'reproducibility_test': reproducibility_test,
			'state': 		desc.get('state'),
			'total price': 	desc.get('totalPrice')
		}
		notes = []

		# if int(np_stage.get('output').get('npeaks_in')) - int(np_stage.get('output').get('npeaks_out')) != int(np_stage.get('output').get('npeaks_rejected')):
		# 	notes.append("in-out!=rej delta=%i" %(int(np_stage.get('output').get('npeaks_in')) - int(np_stage.get('output').get('npeaks_out'))))
		# else:
		# 	notes.append("in-out=rej OK")

		# bb_check_notes = []
		# for stage in [np_stage, gp_stage]:
		# 	bb_dxf = dxpy.DXFile(stage['output']['overlapping_peaks_bb'])
		# 	if int(bb_dxf.describe()['size']) < 200000:
		# 		bb_check_notes.append("%s bb size=%i" %(stage['name'], int(bb_dxf.describe()['size'])))
		# if not bb_check_notes:
		# 	notes.append("bb check OK")
		# else:
		# 	notes.append(bb_check_notes)




		if notes:
			row.update({'notes': '%s' %(notes)})
		else:
			row.update({'notes': '%s' %('OK')})
		#log = subprocess.check_output('dx watch %s' %(analysis.))
		writer.writerow(row)

if __name__ == '__main__':
	main()
