#!/usr/bin/env python2

import os
import subprocess
import logging
import dxpy
import common
import re
from pprint import pprint, pformat
import sys

EPILOG = '''Notes:

Examples:

    %(prog)s
'''

ASSEMBLY_METADATA = {
    'mm10': {
        'gsize': 'mm',
        'csizes': "ENCODE Reference Files:/mm10/mm10_no_alt.chrom.sizes",
        'blacklist': "ENCODE Reference Files:/mm10/blacklists/mm10.blacklist.bed.gz"
    },
    'mm10-minimal': {
        'gsize': 'mm',
        'csizes': "ENCODE Reference Files:/mm10/male.mm10.chrom.sizes",
        'blacklist': "ENCODE Reference Files:/mm10/blacklists/mm10.blacklist.bed.gz"
    },
    'GRCh38': {
        'gsize': 'hs',
        'csizes': "ENCODE Reference Files:/GRCh38/GRCh38_EBV.chrom.sizes",
        'blacklist': "ENCODE Reference Files:/GRCh38/blacklists/GRCh38.blacklist.bed.gz"
    },
    'hg19': {
        'gsize': 'hs',
        'csizes': "ENCODE Reference Files:/hg19/male.hg19.chrom.sizes",
        'blacklist': "ENCODE Reference Files:/hg19/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz"
    }
}


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

    parser.add_argument('experiments', help="Experiment accessions", nargs="*")
    parser.add_argument('--infile', help="File with experiment accessions")
    parser.add_argument('--debug', help="Print debug messages", default=False, action='store_true')
    parser.add_argument('--project', help="Project name or ID", default=dxpy.WORKSPACE_ID)
    parser.add_argument('--outf', help="Output folder name or ID", default="/")
    parser.add_argument('--inf', nargs='*', help="Folder(s) name or ID with tagAligns", default="/")
    parser.add_argument('--yes',   help="Run the workflows created",            default=False, action='store_true')
    parser.add_argument('--tag',   help="String to add to the workflow title")
    parser.add_argument('--key', help="The local keypair identifier from the keyfile.  Default is --key=default", default='default')
    parser.add_argument('--keyfile', default=os.path.expanduser("~/keypairs.json"), help="The local keypair file.  Default is --keyfile=%s" %(os.path.expanduser("~/keypairs.json")))
    parser.add_argument('--gsize', help="Genome size string for MACS2, e.g. mm or hs", default=None)
    parser.add_argument('--csizes', help="chrom.sizes file for bedtobigbed, e.g. ENCODE Reference Files:/mm10/male.mm10.chrom.sizes", default=None)
    parser.add_argument('--blacklist', help="Regions to exclude from final peaks list", default=None)
    parser.add_argument('--assembly', help="Genome assembly, e.g. hg19, mm10, GRCh38", required=True)
    # parser.add_argument('--idr', help="Run IDR. If not specified, run IDR for non-histone targets.", default=False, action='store_true')
    # parser.add_argument('--idrversion', help="IDR version (relevant only if --idr is specified", default="2")
    parser.add_argument('--dryrun', help="Formulate the run command, but don't actually run", default=False, action='store_true')
    parser.add_argument('--spp_version', help="spp version", default="1.14")
    parser.add_argument('--control', help="Use specified control tagAlign rather than infer one.", default=None)
    parser.add_argument('--accession', help="Accession the results to the ENCODE Portal", default=False, action='store_true')
    parser.add_argument('--fqcheck', help="If --accession, check that analysis is based on latest fastqs on ENCODEd", type=t_or_f, default=None)
    parser.add_argument('--skip_control', help="If --accession, accession no control files or metadata", type=t_or_f, default=None)
    parser.add_argument('--force_patch', help="Force patching metadata for existing files", type=t_or_f, default=None)

    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(
            format='%(levelname)s:%(message)s', level=logging.DEBUG)
        logging.debug("Debug logging ON")
    else:  # use the defaulf logging level
        logging.basicConfig(
            format='%(levelname)s:%(message)s')

    return args


def get_control_id(experiment):
    # url = server + '/experiments/%s/' %(exp_id)
    # experiment = encoded_get(url, keypair)
    # fastqs = common.encoded_get('')

    possible_controls = experiment.get('possible_controls')
    if not possible_controls or len(possible_controls) != 1:
        logging.error(
            "Tried to find one possible control, found %s"
            % (possible_controls))
        return None
    return possible_controls[0].get('accession')


def resolve_project(identifier, privs='r'):
    project = dxpy.find_one_project(
        name=identifier, level='VIEW', name_mode='exact', return_handler=True,
        zero_ok=True)
    if project is None:
        try:
            project = dxpy.get_handler(identifier)
        except:
            logging.error(
                'Could not find a unique project with name or id %s'
                % (identifier))
            raise ValueError(identifier)
    logging.debug(
        'Project %s access level is %s'
        % (project.name, project.describe()['level']))
    if privs == 'w' and project.describe()['level'] == 'VIEW':
        logging.error('Output project %s is read-only' % (identifier))
        raise ValueError(identifier)
    return project


def get_all_tas(experiment, default_project, ta_folders):
    logging.debug(
        'get_all_tas: enter with experiment %s default_project %s and ta_folders %s'
        % (experiment.get('accession'), default_project, ta_folders))
    exp_id = experiment['accession']
    possible_files = []
    for base_folder in ta_folders:
        if ':' in base_folder:
            project_name, base_path = base_folder.split(':')
            project = resolve_project(project_name)
            project = project.get_id()
            project_name += ":"
        else:
            project = default_project
            project_name = ""
            base_path = base_folder
        if not base_path.startswith('/'):
            base_path = '/' + base_path
        if not base_path.endswith('/'):
            base_path = base_path + '/'
        path = base_path + exp_id + '/'
        logging.debug(
            "get_all_tas: find_data objects in project %s project_name %s path %s"
            % (project, project_name, path))
        for dxfile in dxpy.find_data_objects(classname='file', state='closed', folder=path, describe=True, recurse=True, project=project):
            desc = dxfile.get('describe')
            logging.debug(
                "get_all_tas: checking object for match: folder %s name %s"
                % (desc.get('folder'), desc.get('name')))
            if exp_id in desc.get('folder') and '/bams' in desc.get('folder') and desc.get('name').endswith(('tagAlign', 'tagAlign.gz')):
                possible_files.append(desc)
    logging.debug(
        "get_all_tas: exit with possible_files %s" % (possible_files))
    return possible_files


def get_rep_ta(experiment, repn, default_project, ta_folders):
    logging.debug(
        'in get_rep_ta with experiment %s repn %s default_project %s and ta_folders %s'
        % (experiment.get('accession'), repn, default_project, ta_folders))
    exp_id = experiment['accession']

    possible_files = get_all_tas(experiment, default_project, ta_folders)
    print("%s %i possible files" % (exp_id, len(possible_files)))

    folders = [f.get('folder') for f in possible_files]
    rep_folder = \
        [folder for folder in folders if folder.endswith('rep%d' % (repn))]
    if len(rep_folder) != 1:
        logging.warning("Could not find folder rep%d" % (repn))
        return None
    rep_folder = rep_folder[0]
    rep_files = [f for f in possible_files if f.get('folder') == rep_folder]

    if len(rep_files) != 1:
        logging.warning(
            "Tried to find one rep%d ta, found %d" % (repn, len(rep_files)))
    if len(rep_files) > 0:
        if len(rep_files) > 1:
            logging.warning("Using first one found")
        return rep_files[0]
    else:
        return None


def get_possible_ctl_ta(experiment, repn, server, keypair, default_project,
                        ta_folders, used_control_ids):

    # Build a list of the possible_control experiments
    possible_control_experiments = []
    for uri in experiment.get('possible_controls'):
        possible_control_experiment = common.encoded_get(server+uri, keypair)
        target_uri = possible_control_experiment.get('target')
        # For now only use controls with no target or target "Control"
        # (i.e. not IgG)
        if not target_uri or target_uri.split('/')[2].startswith('Control'):
            possible_control_experiments.append(possible_control_experiment)
        elif 'control' in target_uri.split('/')[2]:
            logging.warning(
                '%s: possible control %s has target %s.  Such controls are allowed but deprecated.' % (
                    experiment.get('accession'),
                    possible_control_experiment.get('accession'),
                    target_uri))
            possible_control_experiments.append(possible_control_experiment)

        else:
            logging.warning(
                '%s: possible control %s has target %s, not "Control". Skipping.' % (
                    experiment.get('accession'),
                    possible_control_experiment.get('accession'),
                    target_uri))
    logging.debug(pformat(possible_control_experiments))
    try:
        matching_ta = \
            next(ta for ta in [get_rep_ta(e, repn, default_project, ta_folders) for e in possible_control_experiments] if ta and ta['id'] not in used_control_ids)
    except StopIteration:
        logging.warning('Failed to find control rep with matching repn')
        matching_ta = None
    except:
        raise
    else:
        return matching_ta

    tas = []
    for e in possible_control_experiments:
        unused_tas = [ta for ta in get_all_tas(e, default_project, ta_folders) if ta and ta['id'] not in used_control_ids]
        logging.debug(
            'get_possible_ctl_ta: experiment %s unused_tas %s'
            % (e.get('accession'), unused_tas))
        if unused_tas:
            tas.extend(unused_tas)
    if len(tas) > 1:
        logging.warning(
            'Found multiple tas %s, returning first one'
            % ([ta.get('project') + ':' + ta.get('folder') + '/' + ta.get('name') for ta in tas]))
    if tas:
        return tas[0]
    else:
        logging.warning('Failed to find any possible controls that have not already been used')

    tas = []
    for e in possible_control_experiments:
        all_tas = [ta for ta in get_all_tas(e, default_project, ta_folders) if ta]
        logging.debug(
            'get_possible_ctl_ta: experiment %s all_tas %s'
            % (e.get('accession'), all_tas))
        if all_tas:
            tas.extend(all_tas)
    if len(tas) > 1:
        logging.warning(
            'Found multiple tas %s, returning first one'
            % ([ta.get('project') + ':' + ta.get('folder') + '/' + ta.get('name') for ta in tas]))
    if tas:
        return tas[0]
    else:
        logging.error('Failed to find any possible control_tas')
        return None


def get_encffs(s):
    return re.findall('ENCFF[0-9]{3}[A-Z]{3}', s)


def get_ta_from_accessions(accessions, default_project, ta_folders):
    possible_files = []
    for base_folder in ta_folders:
        if ':' in base_folder:
            project_name, path = base_folder.split(':')
            project = resolve_project(project_name)
            project_id = project.get_id()
            project_name += ":"
        else:
            project_id = default_project
            project_name = ""
            path = base_folder
        if not path.startswith('/'):
            path = '/' + path
        if not path.endswith('/'):
            path += '/'
        logging.debug(
            "Looking for TA's in %s %s %s" % (project_id, project_name, path))
        for dxfile in dxpy.find_data_objects(
            classname='file',
            state='closed',
            folder=path + 'bams/',
            project=project_id,
            describe=True,
            recurse=True,
            name='*tagAlign.gz',
            name_mode='glob'
        ):
            possible_files.append(dxfile.get('describe'))
    matched_files = \
        [f for f in possible_files if all([acc in f['name'] for acc in accessions])]
    if not matched_files:
        logging.error(
            'Could not find tagAlign with accessions %s' % (accessions))
        return None
    elif len(matched_files) > 1:
        logging.warning(
            'Found multiple tagAligns that matched accessions %s'
            % (accessions))
        logging.warning(
            'Matched files %s'
            % ([(f['folder'], f['name']) for f in matched_files]))
        logging.warning('Using first one found')
        return matched_files[0]
    else:
        return matched_files[0]


def is_paired_end(tagAlign_file):
    job = dxpy.describe(tagAlign_file['createdBy']['job'])
    return job['output']['paired_end']


def get_tas(experiment, server, keypair, default_project, ta_folders, control_dxhandler):
    # tas = {
    #   'rep1_ta': {
    #       'file_id': "",
    #       'project_id': "",
    #       'folder': "",
    #       'name': "",
    #       'paired_end': False,
    #       'control_path': "",
    #       'enc_repn': 0
    # .for each ta_folder get list of TA's in /ta_folder/bams/ENCSR...
    # .from this list infer repns from the paths ../bams/ENCSR.../repn*
    # .from this list infer the ENCFF's for the fastqs that were used
    # for each repn go to the experiment and find all the fastqs for that rep
    # if there are different fastq's in the experiment, or different reps, warn
    # for each fastq found in the TA filename, find its controlled_by
    # if any have controlled_by, all must have controlled_by else error
    #    gather the list of controlled by and find a TA (anywhere in ta_folders) with those ENCFF's, else error
    # else get possible_controls and try to match the repn, else pick one (rememeber it)
    #   gather the list of fastqs in the possible_controls and find (one) TA with those ENCFF's, else error
    exp_id = experiment['accession']
    possible_files = []
    for base_folder in ta_folders:
        if ':' in base_folder:
            project_name, path = base_folder.split(':')
            project = resolve_project(project_name)
            project_id = project.get_id()
            project_name += ":"
        else:
            project_id = default_project
            project_name = ""
            path = base_folder
        if not path.startswith('/'):
            path = '/' + path
        if not path.endswith('/'):
            path += '/'
        logging.debug(
            "Looking for TA's in %s %s %s" % (project_id, project_name, path))
        for dxfile in dxpy.find_data_objects(
            classname='file',
            state='closed',
            folder=path + 'bams/%s/' % (exp_id),
            project=project_id,
            describe=True,
            recurse=True,
        ):
            desc = dxfile.get('describe')
            if desc.get('name').endswith(('tagAlign', 'tagAlign.gz')):
                possible_files.append(desc)
    logging.debug('Found %s possible files' % (len(possible_files)))
    logging.debug(
        '%s' % ([(f.get('folder'), f.get('name')) for f in possible_files]))
    repns = []
    files_to_ignore = []
    for f in possible_files:
        m = re.search('/rep(\d+)$', f['folder'])
        if m:
            repn = int(m.group(1))
            logging.debug("Matched rep%d" % (repn))
            if repn in repns:
                logging.warning(
                    "Ignoring additional rep%d bam, using first found"
                    % (repn))
                files_to_ignore.append(f)
            else:
                logging.debug("First time finding rep%d" % (repn))
                repns.append(repn)
        else:
            logging.error("Cannot parse rep number from %s" % (f['folder']))
            return None
    for f in files_to_ignore:
        possible_files.remove(f)
    logging.debug('Discovered repns %s' % (repns))
    if len(repns) not in [1, 2]:
        logging.error(
            "Required to have exactly 1 or 2 reps for %s.  Found %d: %s"
            % (exp_id, len(repns), repns))
        return None

    tas = {}
    used_controls = []
    for i, repn in enumerate(repns):
        if control_dxhandler:
            controlled_by_ta_name = control_dxhandler.name
            controlled_by_ta_id = control_dxhandler.get_id()
            controlled_by_accessions = None
        else:
            encode_files = \
                [common.encoded_get(server+'/files/%s/' % (f), keypair) for f in get_encffs(possible_files[i].get('name'))]
            controlled_by = \
                common.flat([f.get('controlled_by') for f in encode_files])
            if any(controlled_by):
                controlled_by_accessions = \
                    list(set([uri.split('/')[2] for uri in controlled_by if uri]))
                controlled_by_ta = \
                    get_ta_from_accessions(
                        controlled_by_accessions, default_project, ta_folders)
                if controlled_by_ta:
                    controlled_by_ta_name = controlled_by_ta.get('name')
                    controlled_by_ta_id = controlled_by_ta.get('id')
                else:
                    logging.error(
                        "%s: Could not find controlled_by_ta for accessions %s"
                        % (experiment.get('accession'), controlled_by_accessions))
                    controlled_by_ta_name = None
                    controlled_by_ta_id = None
            else:
                # evaluate possible controls
                controlled_by_accessions = None
                possible_controls = experiment.get('possible_controls')
                logging.warning(
                    '%s: No controlled_by for rep%d, attempting to infer from possible_controls %s'
                    % (experiment.get('accession'), repn, possible_controls))
                if not possible_controls or not any(possible_controls):
                    logging.error(
                        '%s: Could not find controlled_by or resolve possible_controls for rep%d'
                        % (experiment.get('accession'), repn))
                    controlled_by_ta_name = None
                    controlled_by_ta_id = None
                else:
                    control_ta = get_possible_ctl_ta(
                        experiment, repn, server, keypair, default_project,
                        ta_folders, used_controls)
                    controlled_by_ta_name = control_ta.get('name')
                    controlled_by_ta_id = control_ta.get('id')
        if not (controlled_by_ta_name and controlled_by_ta_id):
            logging.error(
                '%s: Failed to find controls' % (experiment.get('accession')))
            return None
        if controlled_by_ta_id and controlled_by_ta_id in used_controls:
            logging.warning(
                '%s: Using same control %s for multiple reps'
                % (controlled_by_ta_id, controlled_by_ta_name))
        used_controls.append(controlled_by_ta_id)
        # if encode repns are 1,2 then let the pipline input rep numbers (1 or 2) be the same.
        # Otherwise the mapping is arbitrary, but at least do it with smaller rep number first.
        if repn == min(repns):
            ta_index = 1
        else:
            ta_index = 2
        tas.update({
            'rep%d_ta' % (ta_index): {
                'file_id': possible_files[i].get('id'),
                'project_id': possible_files[i].get('project'),
                'folder': possible_files[i].get('folder'),
                'file_name': possible_files[i].get('name'),
                'enc_fqs': get_encffs(possible_files[i].get('name')),
                'controlled_by': controlled_by_accessions,
                'controlled_by_name': controlled_by_ta_name,
                'control_id': controlled_by_ta_id,
                'enc_repn': repn,
                'paired_end': is_paired_end(possible_files[i])
                }
            }
        )

    return tas

    #for each biorepn
    #   look for one TA in /ta_folder/bams/ENCSR.../repn else error
    #   get the list of controlled_by fastq's for all fastqs in that rep
    #   look for one TA in /ta_folder anywhere with (all) those fastqs in its name
    #   


    #   find exactly one tagAlign file anywhere in ta_folders that has each ENCFF in the filename (write a function for that) else error
    #   find the list of fastqs in controlled_by.
    #   If not empty
    #       find exactly one tagAlign anywhere in ta_folders with those ENCFFs.
    #   If empty
    #       Repeat the process above for each ENCSR in possible controls
    #
    #if none, error.  if more than one, error
    #each fastq in the rep should have the same controlled_by (or none)
    #for each rep

    #get the repn's from the experiment and the fastqs for each rep
    #for each repn look for the TA's in /ta_folder/bams/ENCSRxyzabc/repn
    #check to make sure there is exactly one TA with each fastq in its name (form ENCFF.....-ENCFF.....-ENCFF....*)
    #if not error
    #for each fastq, check to see if it has controlled_by.  If any do, all must, else error.
    #   if controlled-by,
    #       check that all fastqs in each biorep have same controlled_by else error
    #       

    #ta_folders assumed to have the structure /ta_folder/bams/ENCSRxyzabc/repn/*.tagAlign or tagAlign.gz
    #can only get folders in DNAnexus from files, so find all the ta paths in each /ta_folder/bams/ENCSRxyzabc/
    #from that list of files, infer the repn's from the file paths
    #for each repn, get the fastq's for that experiment's reps
    #if all the fastqs for the rep have controlled_by then figure out the control experiment accession and control_repn and find the path to that control
    #if none of the fastqs for the rep have controlled_by then get possible_controls and try to match the repn in the control, if no match, pick one
    #if there are multiple possible controls, pick one by heuristic (target starts with Control)
    #print that list of files to look through those file's tafor each folder in ta_folders look for ta_folder/bams/exp_id/repn
    #find the 


def resolve_dx_file(identifier):
    try:
        handler = dxpy.get_handler(identifier)
    except dxpy.DXError:
        try:
            handler = dxpy.find_one_data_object(
                classname='file',
                name=identifier,
                return_handler=True,
                zero_ok=False,
                more_ok=False)
        except dxpy.DXSearchError:
            logging.error('Failed to resolve control %s to unique dx object.  ID or name does not exist or multiple files of that name were found.' % (str(identifier)))
            return None
        else:
            return handler
    else:
        return handler


def main():
    args = get_args()
    authid, authpw, server = common.processkey(args.key, args.keyfile)
    keypair = (authid, authpw)

    experiments = []
    if args.experiments:
        experiments.extend(args.experiments)
    if args.infile:
        with open(args.infile, 'r') as fh:
            experiments.extend([e for e in fh])

    if args.control:
        control_dxhandler = resolve_dx_file(args.control)
    else:
        control_dxhandler = None

    for exp_id in experiments:
        if exp_id.startswith('#'):
            continue
        exp_id = exp_id.rstrip()
        print("Experiment %s" % (exp_id))
        experiment_url = server + '/experiments/%s/' % (exp_id)
        experiment = common.encoded_get(experiment_url, keypair)
        if experiment.get('target'):
            target_url = server + experiment.get('target')
            target = common.encoded_get(target_url, keypair)
        else:
            logging.error('Experiment has no target ... skipping')
            continue

        print(
            "%s %s %s"
            % (experiment['accession'], target.get('investigated_as'),
               experiment.get('description')))

        tas = get_tas(experiment, server, keypair, args.project, args.inf, control_dxhandler)
        if not tas:
            logging.error(
                'Failed to resolve all tagaligns for %s'
                % (experiment['accession']))
            continue
        if not tas.get('rep2_ta'):
            simplicate_experiment = True
            print("Simplicate experiment ta's:")
        else:
            simplicate_experiment = False
            print("Replicated experiment ta's:")
        pprint(tas)
        # sys.exit()
        # continue

        for key, value in tas.iteritems():
            if not value:
                logging.error('Missing %s ... skipping' % (key))
                continue

        workflow_title = '%s Peaks' % (exp_id)
        if args.tag:
            workflow_title += ' %s' % (args.tag)
        outf = args.outf

        if not outf.startswith('/') and outf != '/':
            outf = '/'+outf
        if not outf.endswith('/') and outf != '/':
            outf += '/'
        outf += '%s/peaks/' % (exp_id)

        try:
            investigated_as = target['investigated_as']
        except:
            logging.error(
                "%s: Failed to determine target type ... skipping" % (exp_id))
            continue
        else:
            print(investigated_as)

        rep1_pe = tas['rep1_ta']['paired_end']
        if not simplicate_experiment:
            rep2_pe = tas['rep2_ta']['paired_end']
        else:
            rep2_pe = None

        if simplicate_experiment and rep1_pe is None:
            logging.error(
                "%s: Cannot determine paired end: rep1 PE = %s... skipping"
                % (exp_id, rep1_pe))
            continue
        elif not simplicate_experiment and None in [rep1_pe, rep2_pe]:
            logging.error(
                "%s: Cannot determine paired end: rep1 PE = %s, rep2 PE = %s ... skipping"
                % (exp_id, rep1_pe, rep2_pe))
            continue

        if not simplicate_experiment and rep1_pe != rep2_pe:
            logging.error(
                "%s: rep1 PE %s differs from rep2 PE %s ... skipping"
                % (exp_id, rep1_pe, rep2_pe))
            continue

        if any('histone' in target_type for target_type in investigated_as):
            logging.info(
                "%s: Found to be histone.  No blacklist will be used."
                % (exp_id))
            wf_target = 'histone'
            blacklist = None
        else:
            logging.info("Assumed to be tf")
            wf_target = 'tf'
            if not args.blacklist:
                if args.assembly in ASSEMBLY_METADATA:
                    blacklist = ASSEMBLY_METADATA[args.assembly]['blacklist']
                else:
                    logging.warning(
                        "%s: No blacklist for assembly %s, proceeding with no blacklist"
                        % (exp_id, args.assembly))
                    blacklist = None

        if not args.gsize:
            if args.assembly in ASSEMBLY_METADATA:
                genomesize = ASSEMBLY_METADATA[args.assembly]['gsize']
            else:
                logging.error(
                    "%s: Must specify -gsize for assembly %s"
                    % (exp_id, args.assembly))
        else:
            genomesize = args.gsize

        if not args.csizes:
            if args.assembly in ASSEMBLY_METADATA:
                chrom_sizes = ASSEMBLY_METADATA[args.assembly]['csizes']
            else:
                logging.error(
                    "%s: Must specify -csizes for assembly %s"
                    % (exp_id, args.assembly))
        else:
            chrom_sizes = args.csizes

        command_strings = [
            '~/chip-seq-pipeline/dnanexus/chip_workflow.py',
            '--nomap --yes',
            '--target %s' % (wf_target),
            '--title "%s"' % (workflow_title),
            '--outf "%s"' % (outf),
            '--rep1pe %s' % (str(rep1_pe).lower()),
            '--rep1 %s' % (tas['rep1_ta'].get('file_id')),
            '--ctl1 %s' % (tas['rep1_ta'].get('control_id')),
            '--genomesize %s --chrom_sizes "%s"' % (genomesize, chrom_sizes),
            '--spp_version %s' % (args.spp_version)
        ]

        if not simplicate_experiment:
            command_strings.extend([
                '--rep2pe %s' % (str(rep2_pe).lower()),
                '--rep2 %s' % (tas['rep2_ta'].get('file_id')),
                '--ctl2 %s' % (tas['rep2_ta'].get('control_id')),
            ])

        if blacklist:
            command_strings.append('--blacklist "%s"' % (blacklist))
        if args.debug:
            command_strings.append('--debug')
        if args.accession:
            command_strings.append('--accession')
            if args.fqcheck is not None:
                command_strings.append('--fqcheck=%s' % (args.fqcheck))
            if args.skip_control is not None:
                command_strings.append('--skip_control=%s' % (args.skip_control))
            if args.force_patch is not None:
                command_strings.append('--force_patch=%s' % (args.force_patch))
        run_command = ' '.join(command_strings)
        print(run_command)

        if args.dryrun:
            logging.info('Dryrun')
        else:
            try:
                subprocess.check_call(run_command, shell=True)
            except subprocess.CalledProcessError as e:
                logging.error(
                    "%s: chip_workflow exited with non-zero code %d"
                    % (exp_id, e.returncode))
            else:
                print("%s workflow created" % (experiment['accession']))
                logging.debug(
                    "%s: patching internal_status to url %s"
                    % (exp_id, experiment_url))
                r = common.encoded_patch(
                    experiment_url, keypair, {'internal_status': 'processing'},
                    return_response=True)
                try:
                    r.raise_for_status()
                except:
                    logging.error(
                        "%s: Tried but failed to update experiment internal_status to processing"
                        % (exp_id))
                    logging.error(r.text)


if __name__ == '__main__':
    main()
