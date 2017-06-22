#!/usr/bin/env python2

import os.path
import sys
import subprocess
import logging
import re
import json
import urlparse
import requests
import csv
import StringIO
import common
import dxpy
import time
import pprint

EPILOG = '''Notes:

Examples:

    %(prog)s
'''
FILE_STATUSES_TO_MAP = ['in progress', 'released', 'uploading']
FILE_FORMATS_TO_MAP = ['fastq', 'fasta', 'sra']
#DEFAULT_APPLET_PROJECT = 'E3 ChIP-seq'
DEFAULT_APPLET_PROJECT = dxpy.WORKSPACE_ID
DEFAULT_OUTPUT_PROJECT = dxpy.WORKSPACE_ID
DEFAULT_OUTPUT_FOLDER = '/'

INPUT_SHIELD_APPLET_NAME = 'input_shield'
MAPPING_APPLET_NAME = 'encode_map'
FILTER_QC_APPLET_NAME = 'filter_qc'
XCOR_APPLET_NAME = 'xcor'
POOL_APPLET_NAME = 'pool'
ACCESSION_ANALYSIS_APPLET_NAME = 'accession_analysis'

REFERENCES = [
    {'assembly': 'GRCh38-minimal', 'organism': 'human', 'sex': 'male',   'file': 'ENCODE Reference Files:/Deprecated/GRCh38/GRCh38_minimal_XY.tar.gz'},
    {'assembly': 'GRCh38-minimal', 'organism': 'human', 'sex': 'female', 'file': 'ENCODE Reference Files:/Deprecated/GRCh38/GRCh38_minimal_X.tar.gz'},
    #warning these are not sex-specific yet
    {'assembly': 'GRCh38', 'organism': 'human', 'sex': 'male',   'file': 'ENCODE Reference Files:/GRCh38/ChIP-seq/GCA_000001405.15_GRCh38_no_alt_analysis_set.bwa.tar.gz'},
    {'assembly': 'GRCh38', 'organism': 'human', 'sex': 'female', 'file': 'ENCODE Reference Files:/GRCh38/ChIP-seq/GCA_000001405.15_GRCh38_no_alt_analysis_set.bwa.tar.gz'},
    #warning these are not sex-specific yet
    {'assembly': 'mm10-minimal',   'organism': 'mouse', 'sex': 'male',   'file': 'ENCODE Reference Files:/Deprecated/mm10/ChIP-seq/male.mm10.tar.gz'},
    {'assembly': 'mm10-minimal',   'organism': 'mouse', 'sex': 'female', 'file': 'ENCODE Reference Files:/Deprecated/mm10/ChIP-seq/female.mm10.tar.gz'},
    {'assembly': 'mm10',   'organism': 'mouse', 'sex': 'male', 'file': 'ENCODE Reference Files:/mm10/ChIP-seq/mm10_no_alt_analysis_set_ENCODE.tar.gz'},
    {'assembly': 'hg19',   'organism': 'human', 'sex': 'male',   'file': 'ENCODE Reference Files:/hg19/ChIP-seq/male.hg19.tar.gz'},
    {'assembly': 'hg19',   'organism': 'human', 'sex': 'female', 'file': 'ENCODE Reference Files:/hg19/ChIP-seq/female.hg19.tar.gz'}
    ]

APPLETS = {}


def get_args():
    import argparse
    parser = argparse.ArgumentParser(
        description=__doc__, epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    def t_or_f(arg):
        ua = str(arg).upper()
        if ua == 'TRUE'[:len(ua)]:
            return True
        elif ua == 'FALSE'[:len(ua)]:
            return False
        else:
            assert not (True or False), "Cannot parse %s to boolean" % (arg)

    parser.add_argument(
        'experiments',
        help='List of ENCSR accessions. Can be ENCSR,biorep_i,biorep_j,.. to restrict mapping to specific replicate(s).',
        nargs='*',
        default=None)

    parser.add_argument(
        '--infile',
        help='File containing ENCSR accessions',
        type=argparse.FileType('rU'),
        default=sys.stdin)

    parser.add_argument(
        '--assembly',
        help="Reference genome assembly, e.g. GRCh38, hg19, or mm10")

    parser.add_argument(
        '--sex_specific',
        help="Mapping should be to male or female reference.  Default male.",
        default=False, action='store_true')

    parser.add_argument(
        '--debug',
        help="Print debug messages",
        default=False, action='store_true')

    parser.add_argument(
        '--outp',
        help="Output project name or ID",
        default=DEFAULT_OUTPUT_PROJECT)

    parser.add_argument(
        '--outf',
        help="Output folder name or ID",
        default=DEFAULT_OUTPUT_FOLDER)

    parser.add_argument(
        '--applets',
        help="Name of project containing applets",
        default=DEFAULT_APPLET_PROJECT)

    parser.add_argument(
        '--key',
        help="The local keypair identifier from the local keyfile. Default is 'default'",
        default='default')

    parser.add_argument(
        '--keyfile',
        help="The local keypair filename.",
        default=os.path.expanduser("~/keypairs.json"))

    parser.add_argument(
        '--yes',
        help="Run the workflows created",
        default=False,
        action='store_true')

    parser.add_argument(
        '--raw',
        help="Produce only raw (unfiltered) bams",
        default=False,
        action='store_true')

    parser.add_argument(
        '--tag',
        help="String to add to the workflow name",
        default="")

    parser.add_argument(
        '--no_sfn_dupes',
        help="Disallow duplicte submitted filenames.  Otherwise warn but use files anyway.",
        default=False,
        action='store_true')

    parser.add_argument(
        '--force_se',
        help="Map only read1's of PE sequencing, and combine with SE data.",
        default=False,
        action='store_true')

    parser.add_argument(
        '--crop_length',
        help="Length to crop reads before mapping.  Default is native (no crop)",
        default='native')

    parser.add_argument(
        '--use_existing_folders',
        help="Reuse existing folders even if results have already been saved there",
        default=False,
        action='store_true')

    parser.add_argument(
        '--spp_version',
        help="Version of spp to use for the cross-correlation analysis",
        default='1.14')

    parser.add_argument('--accession', help="Accession the results to the ENCODE Portal", default=False, action='store_true')
    parser.add_argument('--fqcheck', help="If --accession, check that analysis is based on latest fastqs on ENCODEd", type=t_or_f, default=None)
    parser.add_argument('--force_patch', help="Force patching metadata for existing files", type=t_or_f, default=None)
    parser.add_argument('--scrub', help="Scrub bam files of sequence information", type=t_or_f, default=None)
    parser.add_argument('--encoded_check', help="Value is passed on to accession_analysis", type=t_or_f, default=None)

    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(
            format='%(levelname)s:%(message)s', level=logging.DEBUG)
        logging.debug("Logging level set to DEBUG")
    else:  # use the defaulf logging level
        logging.basicConfig(
            format='%(levelname)s:%(message)s')
        logging.info("Logging level set to detault")

    return args


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


def create_folder(project, folder_name):
    try:
        project.new_folder(folder_name, parents=True)
    except:
        logging.error(
            "Cannot create folder %s in project %s"
            % (folder_name, project.name))
        return None
    else:
        logging.info(
            "New folder %s created in project %s"
            % (folder_name, project.name))
        return folder_name


def resolve_folder(project, identifier):
    if not identifier.startswith('/'):
        identifier = '/' + identifier
    try:
        project.list_folder(identifier)
    except:
        return None
    else:
        return identifier


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

def files_to_map(exp_obj, server, keypair, no_sfn_dupes):
    if not exp_obj or not (exp_obj.get('files') or exp_obj.get('original_files')):
        logging.warning('Experiment %s or experiment has no files' %(exp_obj.get('accession')))
        files = []
    else:
        files = []
        for file_uri in exp_obj.get('original_files'):
            file_obj = common.encoded_get(urlparse.urljoin(server, file_uri), keypair=keypair)
            if file_obj.get('status') in FILE_STATUSES_TO_MAP and \
                    file_obj.get('output_type') == 'reads' and \
                    file_obj.get('file_format') in FILE_FORMATS_TO_MAP and \
                    file_obj.get('replicate'):
                if file_obj.get('submitted_file_name') in filenames_in(files):
                    if no_sfn_dupes:
                        logging.error('%s:%s Duplicate submitted_file_name found, skipping that file.' %(exp_obj.get('accession'),file_obj.get('accession')))
                    else:
                        logging.warning('%s:%s Duplicate submitted_file_name found, but allowing duplicates.' %(exp_obj.get('accession'),file_obj.get('accession')))
                        files.extend([file_obj])
                else:
                    files.extend([file_obj])
            elif file_obj.get('output_type') == 'reads' and \
                file_obj.get('file_format') in FILE_FORMATS_TO_MAP and not file_obj.get('replicate'):
                logging.error('%s: Reads file has no replicate' %(file_obj.get('accession')))
    print('returning from files_to_map with %s' % (pprint.pformat(files)))
    return files

def replicates_to_map(files, server, keypair, map_only_reps=[]):
    if not files:
        return []
    else:
        replicate_objects = []
        for f in files:
            replicate = common.encoded_get(urlparse.urljoin(server,f.get('replicate')),keypair)
            if not replicate in replicate_objects:
                if not map_only_reps or (map_only_reps and replicate['biological_replicate_number'] in map_only_reps):
                    replicate_objects.append(replicate)

        return replicate_objects

def choose_reference(experiment, biorep_n, server, keypair, sex_specific):

    replicates = [common.encoded_get(urlparse.urljoin(server,rep_uri), keypair, frame='embedded') for rep_uri in experiment['replicates']]
    replicate = next(rep for rep in replicates if rep.get('biological_replicate_number') == biorep_n)
    logging.debug('Replicate uuid %s' %(replicate.get('uuid')))
    organism_uri = replicate.get('library').get('biosample').get('organism')
    organism_obj = common.encoded_get(urlparse.urljoin(server,organism_uri), keypair)

    try:
        organism_name = organism_obj['name']
    except:
        logging.error('%s:rep%d Cannot determine organism.' %(experiment.get('accession'), biorep_n))
        raise
        return None
    else:
        logging.debug("Organism name %s" %(organism_name))

    if sex_specific:
        try:
            sex = replicate.get('library').get('biosample').get('sex')
            assert sex in ['male', 'female']
        except:
            logging.warning('%s:rep%d Sex is %s.  Mapping to male reference.' %(experiment.get('accession'), biorep_n, sex))
            sex = 'male'

        logging.debug('Organism %s sex %s' %(organism_name, sex))
    else:
        sex = 'male'
    
    genome_assembly = args.assembly

    reference = next((ref.get('file') for ref in REFERENCES if ref.get('organism') == organism_name and ref.get('sex') == sex and ref.get('assembly') == genome_assembly), None)
    logging.debug('Found reference %s' %(reference))
    return reference


def build_workflow(experiment, biorep_n, input_shield_stage_input, accession, use_existing_folders):

    output_project = resolve_project(args.outp, 'w')
    logging.debug('Found output project %s' % (output_project.name))

    applet_project = resolve_project(args.applets, 'r')
    logging.debug('Found applet project %s' % (applet_project.name))

    mapping_applet = \
        find_applet_by_name(MAPPING_APPLET_NAME, applet_project.get_id())
    logging.debug('Found applet %s' % (mapping_applet.name))

    input_shield_applet = \
        find_applet_by_name(INPUT_SHIELD_APPLET_NAME, applet_project.get_id())
    logging.debug('Found applet %s' % (input_shield_applet.name))

    folders = ['workflows', 'fastqs', 'raw_bams', 'bams']
    folder_paths = \
        ['/'.join([args.outf,
                   folder_name,
                   experiment.get('accession'),
                   'rep%d' % (biorep_n)])
         for folder_name in folders]
    paths_exist = \
        [resolve_folder(output_project, folder_path)
         for folder_path in folder_paths
         if resolve_folder(output_project, folder_path)]
    if any(paths_exist):
        msg = "%s: output paths already exist: %s" % (experiment.get('accession'), paths_exist)
        if use_existing_folders:
            logging.warning(msg)
        else:
            msg += "\nUse --use_existing_folders to supress but possibly create duplicate files"
            logging.error(msg)
            return None
    workflow_output_folder, fastq_output_folder, mapping_output_folder, final_output_folder = \
        tuple(create_folder(output_project, folder_path)
              for folder_path in folder_paths)

    if args.raw:
        workflow_title = \
            ('Map %s rep%d to %s (no filter)'
             % (experiment.get('accession'), biorep_n, args.assembly))
        workflow_name = 'ENCODE raw mapping pipeline'
    else:
        workflow_title = \
            ('Map %s rep%d to %s and filter'
             % (experiment.get('accession'), biorep_n, args.assembly))
        workflow_name = 'ENCODE mapping pipeline'

    if args.tag:
        workflow_title += ': %s' % (args.tag)

    workflow = dxpy.new_dxworkflow(
        title=workflow_title,
        name=workflow_name,
        project=output_project.get_id(),
        folder=workflow_output_folder
    )

    input_shield_stage_id = workflow.add_stage(
        input_shield_applet,
        name='Gather inputs %s rep%d' % (experiment.get('accession'), biorep_n),
        folder=fastq_output_folder,
        stage_input=input_shield_stage_input
    )

    input_names = \
        [name for name in ['reads1', 'reads2', 'crop_length', 'reference_tar',
         'bwa_version', 'bwa_aln_params', 'samtools_version', 'debug']
         if name in input_shield_stage_input]
    logging.debug('input_names: %s' % (input_names))
    mapping_stage_input = dict(zip(
        input_names,
        [dxpy.dxlink(
            {'stage': input_shield_stage_id, 'outputField': input_name})
         for input_name in input_names]))
    logging.debug('mapping_stage_input: %s' % (mapping_stage_input))
    mapping_stage_id = workflow.add_stage(
        mapping_applet,
        name='Map %s rep%d' % (experiment.get('accession'), biorep_n),
        folder=mapping_output_folder,
        stage_input=mapping_stage_input
    )

    if not args.raw:
        filter_qc_applet = \
            find_applet_by_name(FILTER_QC_APPLET_NAME, applet_project.get_id())
        logging.debug('Found applet %s' % (filter_qc_applet.name))

        filter_qc_stage_id = workflow.add_stage(
            filter_qc_applet,
            name='Filter and QC %s rep%d' % (experiment.get('accession'), biorep_n),
            folder=final_output_folder,
            stage_input={
                'input_bam': dxpy.dxlink({'stage': mapping_stage_id, 'outputField': 'mapped_reads'}),
                'paired_end': dxpy.dxlink({'stage': mapping_stage_id, 'outputField': 'paired_end'}),
                'scrub': args.scrub
            }
        )

        xcor_applet = find_applet_by_name(XCOR_APPLET_NAME, applet_project.get_id())
        logging.debug('Found applet %s' %(xcor_applet.name))

        xcor_stage_id = workflow.add_stage(
            xcor_applet,
            name='Calculate cross-correlation %s rep%d' %(experiment.get('accession'), biorep_n),
            folder=final_output_folder,
            stage_input={
                'input_bam': dxpy.dxlink({'stage': filter_qc_stage_id, 'outputField': 'filtered_bam'}),
                'paired_end': dxpy.dxlink({'stage': filter_qc_stage_id, 'outputField': 'paired_end'}),
                'spp_version': args.spp_version
            }
        )


    ''' This should all be done in the shield's postprocess entrypoint
    if args.accession_outputs:
        derived_from = input_shield_stage_input.get('reads1')
        if reads2:
            derived_from.append(reads2)
        files_json = {dxpy.dxlink({'stage': mapping_stage_id, 'outputField': 'mapped_reads'}) : {
            'notes': 'Biorep%d | Mapped to %s' %(biorep_n, input_shield_stage_input.get('reference_tar')),
            'lab': 'j-michael-cherry',
            'award': 'U41HG006992',
            'submitted_by': 'jseth@stanford.edu',
            'file_format': 'bam',
            'output_type': 'alignments',
            'derived_from': derived_from,
            'dataset': experiment.get('accession')}
        }
        output_shield_stage_id = workflow.add_stage(
            output_shield_applet,
            name='Accession outputs %s rep%d' %(experiment.get('accession'), biorep_n),
            folder=mapping_output_folder,
            stage_input={'files': [dxpy.dxlink({'stage': mapping_stage_id, 'outputField': 'mapped_reads'})],
                         'files_json': files_json,
                         'key': input_shield_stage_input.get('key')}
        )
    '''
    return workflow


def map_only(experiment, biorep_n, files, server, keypair, sex_specific,
             crop_length, accession, fqcheck, force_patch,
             use_existing_folders, encoded_check):

    pprint.pprint("in maponly with files %s" % (pprint.pformat(files)))
    if not files:
        logging.debug('%s:%s No files to map' %(experiment.get('accession'), biorep_n))
        return
    #look into the structure of files parameter to decide on pooling, paired end etc.

    workflows = []
    input_shield_stage_input = {}

    reference_tar = choose_reference(experiment, biorep_n, server, keypair, sex_specific)
    if not reference_tar:
        logging.warning('%s:%s Cannot determine reference' %(experiment.get('accession'), biorep_n))
        return

    input_shield_stage_input.update({
        'reference_tar': reference_tar,
        'debug': args.debug,
        'crop_length': crop_length
    })

    if all(isinstance(f, dict) for f in files): #single end
        input_shield_stage_input.update({'reads1': [f.get('accession') for f in files]})
        workflows.append(build_workflow(experiment, biorep_n, input_shield_stage_input, accession, use_existing_folders))
    elif all(isinstance(f, tuple) for f in files): #paired-end
        #launches separate mapping jobs for each readpair
        #TODO: upadte input_shield to take an array of read1/read2 PE pairs then pass that array from here
        input_shield_stage_input.update({'reads1': [], 'reads2': []})
        for readpair in files:
            try:
                input_shield_stage_input['reads1'].append(next(f.get('accession') for f in readpair if f.get('paired_end') == '1'))
                input_shield_stage_input['reads2'].append(next(f.get('accession') for f in readpair if f.get('paired_end') == '2'))
            except StopIteration:
                logging.error('%s rep %s: Unmatched read pairs' %(experiment.get('accession'),biorep_n))
                return []
        workflows.append(build_workflow(experiment, biorep_n, input_shield_stage_input, accession, use_existing_folders))
    else:
        logging.error('%s: List of files to map appears to be mixed single-end and paired-end: %s' %(experiment.get('accession'), files))

    analyses = []
    if args.yes:
        for workflow in [wf for wf in workflows if wf]:
            if args.debug:
                analysis = workflow.run(
                    {},
                    priority='high',
                    debug={'debugOn': ['AppInternalError', 'AppError']},
                    delay_workspace_destruction=True,
                    allow_ssh=['*'])
                analyses.append(analysis)
            else:
                analysis = workflow.run(
                    {},
                    priority='normal')
                analyses.append(analysis)
            if accession:
                applet_project = resolve_project(args.applets, 'r')
                accession_analysis_applet = find_applet_by_name(ACCESSION_ANALYSIS_APPLET_NAME, applet_project.get_id())
                accession_output_folder = '/' + accession_analysis_applet.name
                accession_job_input = {
                    'analysis_ids': [analysis.get_id()],
                    'wait_on_files': []
                }
                if fqcheck is not None:
                    accession_job_input.update({'fqcheck' : fqcheck})
                if force_patch is not None:
                    accession_job_input.update({'force_patch': force_patch})
                if encoded_check is not None:
                    accession_job_input.update({'encoded_check': encoded_check})
                time.sleep(5)
                max_retries = 10
                retries = max_retries
                while retries:
                    try:
                        accession_job = accession_analysis_applet.run(
                            accession_job_input,
                            name='Accession %s' % (analysis.name),
                            folder=accession_output_folder,
                            depends_on=analysis.describe()['dependsOn']
                        )
                    except Exception as e:
                        logging.error("%s launching auto-accession ... %d retries left" % (e, retries))
                        time.sleep(5)
                        retries -= 1
                        continue
                    else:
                        logging.info("Auto-accession will run as %s %s" % (accession_job.name, accession_job.get_id()))
                        break
                else:
                    logging.error("Auto-accession failed with %s" % ())

    return analyses


def main():
    global args
    args = get_args()

    authid, authpw, server = common.processkey(args.key, args.keyfile)
    keypair = (authid, authpw)

    if args.experiments:
        exp_ids = csv.reader(StringIO.StringIO('\n'.join([s.rstrip() for s in args.experiments])))
    else:
        exp_ids = csv.reader(args.infile)

    for row in exp_ids:
        if row[0].startswith('#'):
            continue
        exp_id = row[0].strip()
        if len(row) > 1:
            repns = []
            for s in row[1:]:
                repns.extend(s.split(','))
            map_only_reps = list(set([int(s) for s in repns]))
        else:
            map_only_reps = []
        outstrings = []
        encode_url = urlparse.urljoin(server, exp_id)
        experiment = common.encoded_get(encode_url, keypair)
        outstrings.append(exp_id)
        files = files_to_map(experiment, server, keypair, args.no_sfn_dupes)
        outstrings.append(str(len(files)))
        outstrings.append(str([f.get('accession') for f in files]))
        replicates = replicates_to_map(files, server, keypair, map_only_reps)
        biorep_numbers = \
            set([rep.get('biological_replicate_number') for rep in replicates])
        in_process = False
        print("files %s" % (pprint.pformat(files)))
        if files:
            for biorep_n in biorep_numbers:
                outstrings.append('rep%s' %(biorep_n))
                print("biorep n %s" % (biorep_n))
                biorep_files = [f for f in files if biorep_n in common.biorep_ns(f,server,keypair)]
                print("biorep_files %s" % (pprint.pformat(biorep_files)))
                paired_files = []
                unpaired_files = []
                print("biorep_files %s" % (pprint.pformat(biorep_files)))
                while biorep_files:
                    file_object = biorep_files.pop()
                    if file_object.get('paired_end') == None: # group all the unpaired reads for this biorep together
                        unpaired_files.append(file_object)
                    elif file_object.get('paired_end') in ['1','2']:
                        if file_object.get('paired_with'):
                            mate = next((f for f in biorep_files if f.get('@id') == file_object.get('paired_with')), None)
                        else: #have to find the file that is paired with this one
                            mate = next((f for f in biorep_files if f.get('paired_with') == file_object.get('@id')), None)
                        if mate:
                            biorep_files.remove(mate)
                        else:
                            logging.warning('%s:%s could not find mate' %(experiment.get('accession'), file_object.get('accession')))
                            mate = {}

                        # if mapping as SE, ignore the mate and just map the
                        # rep1 as SE with all the other SE for this rep, if any
                        if args.force_se:
                            unpaired_files.append(next(
                                f for f in [file_object, mate]
                                if f.get('paired_end') == '1'))
                        else:
                            paired_files.append((file_object, mate))

                if biorep_files:
                    logging.warning('%s: leftover file(s) %s' %(experiment.get('accession'), biorep_files))
                if paired_files:
                    pe_jobs = \
                        map_only(experiment, biorep_n, paired_files,
                                 server, keypair, args.sex_specific,
                                 args.crop_length, args.accession,
                                 args.fqcheck, args.force_patch,
                                 args.use_existing_folders, args.encoded_check)
                    in_process = True
                if unpaired_files:
                    se_jobs = \
                        map_only(experiment, biorep_n, unpaired_files,
                                 server, keypair, args.sex_specific,
                                 args.crop_length, args.accession,
                                 args.fqcheck, args.force_patch,
                                 args.use_existing_folders, args.encoded_check)
                    in_process = True
                if paired_files and pe_jobs:
                    outstrings.append('paired:%s' %([(a.get('accession'), b.get('accession')) for (a,b) in paired_files]))
                    outstrings.append('paired jobs:%s' %([j.get_id() for j in pe_jobs]))
                else:
                    outstrings.append('paired:%s' %(None))
                if unpaired_files and se_jobs:
                    outstrings.append('unpaired:%s' %([f.get('accession') for f in unpaired_files]))
                    outstrings.append('unpaired jobs:%s' %([j.get_id() for j in se_jobs]))
                else:
                    outstrings.append('unpaired:%s' %(None))
            if in_process:
                r = common.encoded_patch(encode_url, keypair, {"internal_status": "processing"}, return_response=True)
                try:
                    r.raise_for_status()
                except:
                    logging.error("Tried and failed to set internal_status")
                    logging.error(r.text)
            print('\t'.join(outstrings))
        else:  # no files
            if not replicates:
                logging.warning('%s: No files and no replicates' %experiment.get('accession'))
            else:
                logging.warning('%s: No files to map' %experiment.get('accession'))
        if files and not replicates:
            logging.warning('%s: Files but no replicates' %experiment.get('accession'))

if __name__ == '__main__':
    main()
