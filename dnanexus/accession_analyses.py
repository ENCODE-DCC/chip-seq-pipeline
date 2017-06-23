#!/usr/bin/env python2

import sys
import logging
import time
import subprocess
import shlex
import dxpy

ACCESSION_ANALYSIS_APPLET = '/ChIP-seq/applets/accession_analysis'

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())

timestring = time.strftime("%y%m%d%H%M%S")

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

    def t_or_f(arg):
        ua = str(arg).upper()
        if ua == 'TRUE'[:len(ua)]:
            return True
        elif ua == 'FALSE'[:len(ua)]:
            return False
        else:
            assert not (True or False), "Cannot parse %s to boolean" % (arg)

    parser.add_argument('analysis_ids', help='List of analysis IDs to accession', nargs='*', default=None)
    parser.add_argument('--infile', help='Local file containing analysis IDs', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('--outfile', help='DNAnexus file to save output summary table', default="accession_%s.csv" % (timestring))
    parser.add_argument('--destination', help='DNAnexus folder for output summary table', default='/accession_log/')
    parser.add_argument('--name', help='DNAnexus name for the accessioning run', default=None)
    parser.add_argument('--watch', help="Watch the run log", default=False, action='store_true')
    # Any applet arguments with default None will not be passed to the DNAnexus
    # applet, so that its defaults will be used
    # Otherwise the string will be cast to boolean
    parser.add_argument('--project', help='DNAnexus project in which to run', default=None)
    parser.add_argument('--pipeline', help='Over-ride automatic pipeline determination', default=None)
    parser.add_argument('--key', help="The local keypair identifier from the keyfile.", default=None)
    parser.add_argument('--keyfile', help="The local keyfile.", default=None)
    parser.add_argument('--debug', help="Print debug messages", type=t_or_f, default=None)
    parser.add_argument('--dryrun', help="Set up runs but don't change anything.", type=t_or_f, default=None)
    parser.add_argument('--force_patch', help="Force patching metadata for existing files", type=t_or_f, default=None)
    parser.add_argument('--force_upload', help="Force re-uploading for existing files. Files not in status uploading are skipped", type=t_or_f, default=None)
    parser.add_argument('--fqcheck', help="Check that analysis is based on latest fastqs on ENCODEd", type=t_or_f, default=None)
    parser.add_argument('--accession_raw', help="Accession unfiltered bams", type=t_or_f, default=None)
    parser.add_argument('--signal_only', help="Accession through signal generation only", type=t_or_f, default=None)
    parser.add_argument('--skip_control', help="Accession no control files or metadata", type=t_or_f, default=None)
    parser.add_argument('--encoded_check', help="Check if ENCODE server is ready or not", type=t_or_f, default=None)
    args = parser.parse_args()

    return args


def main():

    args = get_args()
    if args.debug:
        logger.setLevel(logging.DEBUG)
        logger.debug("Logging level set to DEBUG")
    else:
        logger.setLevel(logging.INFO)

    if args.analysis_ids:
        ids = [i for i in args.analysis_ids if not i.startswith('#')]
    elif args.infile:
        ids = [i for i in args.infile if not i.startswith('#')]
    else:
        # never reached because inile defaults to stdin
        raise InputError("Must supply analysis id's in arguments or --infile")

    if not args.name:
        if len(ids) > 1:
            job_name = "batch_%s" % (timestring)
        else:
            analysis = dxpy.DXAnalysis(ids[0])
            job_name = "Accession %s" % (analysis.name)
    else:
        job_name = args.name

    tokens = [
        'dx run %s' % (ACCESSION_ANALYSIS_APPLET),
        '-i "outfn=%s"' % (args.outfile),
        '--destination "%s"' % (args.destination),
        '--name "%s"' % (job_name),
        '--yes'
    ]
    if args.watch:
        tokens.append('--watch')
    if args.project is not None:
        tokens.append('-i "project=%s"' % (args.project))
    if args.pipeline is not None:
        tokens.append('-i "pipeline=%s"' % (args.pipeline))
    if args.key is not None:
        tokens.append('-i "key=%s"' % (args.key))
    # if args.keyfile is not None:
    #     tokens.append('-i "keyfile=%s"' % (args.keyfile))
    if args.debug is not None:
        tokens.append('-i "debug=%s"' % (args.debug))
    if args.dryrun is not None:
        tokens.append('-i "dryrun=%s"' % (args.dryrun))
    if args.force_patch is not None:
        tokens.append('-i "force_patch=%s"' % (args.force_patch))
    if args.force_upload is not None:
        tokens.append('-i "force_upload=%s"' % (args.force_upload))
    if args.fqcheck is not None:
        tokens.append('-i "fqcheck=%s"' % (args.fqcheck))
    if args.accession_raw is not None:
        tokens.append('-i "accession_raw=%s"' % (args.accession_raw))
    if args.signal_only is not None:
        tokens.append('-i "signal_only=%s"' % (args.signal_only))
    if args.skip_control is not None:
        tokens.append('-i "skip_control=%s"' % (args.skip_control))
    if args.encoded_check is not None:
        tokens.append('-i "encoded_check=%s"' % (args.encoded_check))

    for (i, analysis_id) in enumerate(ids):
        if analysis_id.startswith('#'):
            continue
        analysis_id = analysis_id.rstrip()
        logger.debug('%s' % (analysis_id))
        tokens.append('-i "analysis_ids=%s"' % (analysis_id))

    command_string = ' '.join(tokens)
    logger.debug(command_string)
    subprocess.check_call(shlex.split(command_string))

if __name__ == '__main__':
    main()
