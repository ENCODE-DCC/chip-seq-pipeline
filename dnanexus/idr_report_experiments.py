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

    parser.add_argument('experiments',  help='List of experiment accessions to report on', nargs='*', default=None)
    parser.add_argument('--infile',     help='File containing experiment accessions', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('--all',        help='Report on all possible IDR experiments', default=False, action='store_true')
    parser.add_argument('--outfile',    help='csv output', type=argparse.FileType('wb'), default=sys.stdout)
    parser.add_argument('--debug',      help="Print debug messages", default=False, action='store_true')
    parser.add_argument('--key',        help="The keypair identifier from the keyfile.", default='www')
    parser.add_argument('--keyfile',    help="The keyfile.", default=os.path.expanduser("~/keypairs.json"))
    parser.add_argument('--created_after', help="String to search for analyses (instead of looking in --infile or arguments) in the DNAnexus form like -5d", default=None)
    parser.add_argument('--state',      help="One or more analysis states to report on (only with --created_after)", nargs='*', default=["done"])
    parser.add_argument('--lab',        help="One or more labs to limit the reporting to", nargs='*', default=[])
    parser.add_argument('--assembly',   help="Genome assembly to report on (e.g. hg19 or GRCg38", required=True)

    args = parser.parse_args()

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
        logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
        logger.setLevel(logging.DEBUG)
    else:  # use the defaulf logging level
        logging.basicConfig(format='%(levelname)s:%(message)s')
        logger.setLevel(logging.INFO)

    authid, authpw, server = common.processkey(args.key, args.keyfile)
    keypair = (authid,authpw)

    if args.experiments:
        ids = args.experiments
    # elif args.created_after:
    #   analyses = []
    #   for state in args.state:
    #       analyses.extend(dxpy.find_analyses(name="ENCSR*",name_mode='glob',state=state,include_subjobs=True,return_handler=True,created_after="%s" %(args.created_after)))
    #   ids = [analysis.get_id() for analysis in analyses if analysis.describe()['executableName'] == 'tf_chip_seq' or analysis.describe()['executableName'].startswith('ENCSR783QUL Peaks')]
    elif args.all:
        exp_query = \
            "/search/?type=Experiment" + \
            "&assay_title=ChIP-seq" + \
            "&award.project=ENCODE" + \
            "&status=released&status=submitted&status=in+progress&status=started&status=release+ready"
        all_experiments = common.encoded_get(server+exp_query, keypair)['@graph']
        ids = [exp.get('accession') for exp in all_experiments]
    elif args.infile:
        ids = args.infile
    else:
        #never reached because inile defaults to stdin
        raise InputError("Must supply experiment id's in arguments or --infile")

    fieldnames = [  'date','analysis','analysis id','experiment','target','biosample_term_name','biosample_type','lab','rfa','assembly',
                    'Nt','Np','N1','N2','rescue_ratio','self_consistency_ratio','reproducibility_test',
                    'state','release','total price','notes']
    writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, delimiter='\t', quotechar='"')
    writer.writeheader()

    idr_query = \
        "/search/?type=File" + \
        "&assembly=%s" % (args.assembly) + \
        "&file_format=bed" + \
        "&output_type=optimal+idr+thresholded+peaks" + \
        "&output_type=conservative+idr+thresholded+peaks" + \
        "&lab.title=ENCODE+Processing+Pipeline" + \
        "&lab.title=J.+Michael+Cherry,+Stanford" + \
        "&status=in+progress&status=released&status=uploading&status=uploaded"
    all_idr_files = common.encoded_get(server+idr_query, keypair)['@graph']

    for (i, experiment_id) in enumerate(ids):
        if experiment_id.startswith('#'):
            continue
        experiment_id = experiment_id.rstrip()
        experiment_uri = '/experiments/%s/' % (experiment_id)
        idr_files = \
            [f for f in all_idr_files if f['dataset'] == experiment_uri]
        idr_step_runs = set([f.get('step_run') for f in idr_files])
        if not len(idr_step_runs):
            if not args.all:
                logger.warning(
                    "%s: Found %d IDR step runs.  Skipping"
                    % (experiment_id, len(idr_step_runs)))
            continue

        idr_qc_uris = []
        assemblies = []
        for f in idr_files:
            quality_metrics = f.get('quality_metrics')
            if not len(quality_metrics) == 1:
                logger.error(
                    '%s: Expected one IDR quality metric for file %s. Found %d.'
                    % (experiment_id, f.get('accession'), len(quality_metrics)))
            idr_qc_uris.extend(quality_metrics)
            assembly = f.get('assembly')
            if not assembly:
                logger.error(
                    '%s: File %s has no assembly'
                    % (experiment_id, f.get('accession')))
            assemblies.append(assembly)
        idr_qc_uris = set(idr_qc_uris)
        if not len(idr_qc_uris) == 1:
            logger.error(
                '%s: Expected one unique IDR metric, found %d. Skipping.'
                % (experiment_id, len(idr_qc_uris)))
            continue
        assemblies = set(assemblies)
        if not len(assemblies) == 1:
            logger.error(
                '%s: Expected one unique assembly, found %d. Skipping.'
                % (experiment_id, len(assemblies)))
            continue
        assembly = next(iter(assemblies))

        idr_step_run_uri = next(iter(idr_step_runs))
        idr_step_run = common.encoded_get(server+idr_step_run_uri, keypair)
        try:
            dx_job_id_str = idr_step_run.get('dx_applet_details')[0].get('dx_job_id')
        except:
            logger.warning("Failed to get dx_job_id from step_run.dx_applet_details.dx_job_id")
            logger.debug(idr_step_run)
            dx_job_id_str = None #could try to pull it from alias
        dx_job_id = dx_job_id_str.rpartition(':')[2]
        dx_job = dxpy.DXJob(dx_job_id)
        job_desc = dx_job.describe()
        analysis_id = job_desc.get('analysis')

        logger.debug('%s' %(analysis_id))
        analysis = dxpy.DXAnalysis(analysis_id)
        desc = analysis.describe()
        project = desc.get('project')

        m = re.match('^(ENCSR[0-9]{3}[A-Z]{3}) Peaks', desc['name'])
        if m:
            experiment_accession = m.group(1)
        else:
            logger.error("No accession in %s, skipping." % (desc['name']))
            continue

        if args.all:  # we've already gotten all the experiment objects
            experiment = \
                next(e for e in all_experiments
                     if e['accession'] == experiment_accession)
        else:
            experiment = \
                common.encoded_get(urlparse.urljoin(
                    server,
                    '/experiments/%s' % (experiment_accession)), keypair)
        logger.debug('ENCODEd experiment %s' % (experiment['accession']))
        if args.lab and experiment['lab'].split('/')[2] not in args.lab:
            continue



        try:
            idr_stage = next(s['execution'] for s in desc['stages'] if s['execution']['name'] == "Final IDR peak calls")
        except:
            logger.error('Failed to find final IDR stage in %s' %(analysis_id))
        else:
            if idr_stage['state'] != 'done': #Final IDR peak calls stage not done, so loop through intermediate IDR stages to find errors
                Np = N1 = N2 = Nt = rescue_ratio = self_consistency_ratio = reproducibility_test = None
                notes = []
                #note this list contains a mis-spelled form of IDR Pooled Pseudoreplicates because until 11/13/15 the pipeline stage name was misspelled - need to be able to report on those runs
                idr_stage_names = ['IDR True Replicates', 'IDR Rep 1 Self-pseudoreplicates', 'IDR Rep 2 Self-pseudoreplicates', 'IDR Pooled Pseudoreplicates', 'IDR Pooled Pseudoeplicates']
                for stage_name in idr_stage_names:
                    try:
                        idr_stage = next(s['execution'] for s in desc['stages'] if s['execution']['name'] == stage_name)
                    except StopIteration:
                        continue
                    except:
                        raise
                    if idr_stage['state'] == 'failed':
                        try:
                            job_log = subprocess.check_output('dx watch %s' %(idr_stage['id']), shell=True, stderr=subprocess.STDOUT)
                        except subprocess.CalledProcessError as e:
                            job_log = e.output
                        else:
                            job_log = None
                        if job_log:
                            patterns = [r'Peak files must contain at least 20 peaks post-merge']
                            for p in patterns:
                                m = re.search(p,job_log)
                                if m:
                                    notes.append("%s: %s" %(stage_name,m.group(0)))
                        if not notes:
                            notes.append(idr_stage['failureMessage'])
                try:
                    done_time = next(transition['setAt'] for transition in desc['stateTransitions'] if transition['newState'] == "failed")
                except StopIteration:
                    done_time = "Not done or failed"
                except:
                    raise
            else:
                Np = idr_stage['output'].get('Np')
                N1 = idr_stage['output'].get('N1')
                N2 = idr_stage['output'].get('N2')
                Nt = idr_stage['output'].get('Nt')
                rescue_ratio = idr_stage['output'].get('rescue_ratio')
                self_consistency_ratio = idr_stage['output'].get('self_consistency_ratio')
                reproducibility_test = idr_stage['output'].get('reproducibility_test')
                notes = "IDR Complete"
                done_time = next(transition['setAt'] for transition in desc['stateTransitions'] if transition['newState'] == "done")

        if done_time:
            date = time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(done_time/1000))
        else:
            date = "Running"
        analysis_link = 'https://platform.dnanexus.com/projects/%s/monitor/analysis/%s' %(desc.get('project').split('-')[1], desc.get('id').split('-')[1])
        experiment_link = '%sexperiments/%s' %(server, experiment.get('accession'))
        row = {
            'date': date,
            'analysis':     analysis_link,
            'analysis id':  desc.get('id'),
            'experiment':   experiment_link,
            'target':       experiment['target'].split('/')[2],
            'biosample_term_name':  experiment.get('biosample_term_name'),
            'biosample_type':   experiment.get('biosample_type'),
            'lab':          experiment['lab'].split('/')[2],
            'rfa':          common.encoded_get(server+experiment.get('award'),keypair).get('rfa'),
            'assembly':     assembly,
            'Np':           Np,
            'N1':           N1,
            'N2':           N2,
            'Nt':           Nt,
            'rescue_ratio': rescue_ratio,
            'self_consistency_ratio': self_consistency_ratio,
            'reproducibility_test': reproducibility_test,
            'state':        desc.get('state'),
            'release':      experiment['status'],
            'total price':  desc.get('totalPrice')
        }

        if notes:
            row.update({'notes': '%s' %(notes)})
        else:
            row.update({'notes': '%s' %('OK')})
        #log = subprocess.check_output('dx watch %s' %(analysis.))
        writer.writerow(row)

if __name__ == '__main__':
    main()
