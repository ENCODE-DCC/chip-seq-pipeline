#!/usr/bin/env python

import os, sys, subprocess, logging, dxpy, json, re, socket, getpass, urlparse
from posixpath import basename, dirname
import common

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

    parser.add_argument('infile',       help="analysis ID's", nargs='?', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('outfile',      nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('--debug',      help="Print debug messages",                default=False, action='store_true')
    parser.add_argument('--project',    help="Project name or ID",          default=dxpy.WORKSPACE_ID)
    parser.add_argument('--nodownload', help="Don't transfer data files, only make the hub", default=False, action='store_true')
    parser.add_argument('--dxf',        help="Generate DNAnexus URL's. Implies --nodownload", default=False, action='store_true')
    parser.add_argument('--truncate',   help="Replace existing trackDb file", default=False, action='store_true')
    parser.add_argument('--key',        help="The keypair identifier from the keyfile.", default='www')
    parser.add_argument('--ddir',       help="The local directory to store data files", default=os.path.expanduser('~/tracks'))
    parser.add_argument('--tdbpath',    help="The local path to the trackhub trackDb", default=os.path.expanduser('~/tracks/E3_ChIP_hub/mm10/trackDb.txt'))
    parser.add_argument('--turl',       help="The base URL to the tracks", default='http://'+socket.getfqdn()+'/'+getpass.getuser()+'/tracks/')
    parser.add_argument('--tag',        help="A short string to add to the composite track longLabel")
    parser.add_argument('--lowpass',    help="Add replicated peak tracks with peaks less than this(these) width(s)", nargs='*', type=int)
    parser.add_argument('--pipeline',   help="tf or histone.  If omitted, try to determine automatically", default=None)

    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
    else: #use the defaulf logging level
        logging.basicConfig(format='%(levelname)s:%(message)s')

    if args.dxf:
        args.nodownload = True

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
    url = urlparse.urljoin(url,'?format=json&frame=embedded')
    if keypair:
        response = requests.get(url, auth=keypair, headers=HEADERS)
    else:
        response = requests.get(url, headers=HEADERS)
    return response.json()

def pprint_json(JSON_obj):
    print json.dumps(JSON_obj, sort_keys=True, indent=4, separators=(',', ': '))

def composite_stanza(accession, longLabel):
    return(
        "track %s\n" %(accession) + \
        "compositeTrack on\n" + 
        "shortLabel %s\n" %(accession) + \
        "longLabel %s\n" %(longLabel) + \
        "type bed 3\n" + \
        "visibility full\n" + \
        "subGroup1 view Views PK=Peaks SIG=Signals\n\n")

def viewpeaks_stanza(accession):
    return(
        "\ttrack %sviewpeaks\n" %(accession) + \
        "\tparent %s on\n" %(accession) + \
        "\tshortLabel Peaks\n" + \
        "\tlongLabel Peaks\n" + \
        "\tview PK\n" + \
        "\tvisibility dense\n" + \
        "\ttype bigBed 6 +\n" + \
        "\tscoreFilter 0\n" + \
        "\tscoreFilterLimits 0:1000\n" + \
        "\tviewUi on\n\n")

def peaks_stanza(accession, url, name, n, tracktype='bigBed 6 +', lowpass=[], dx=None):
    return_string = \
        "\t\ttrack %s%d\n" %(accession,n) + \
        "\t\tbigDataUrl %s\n" %(url) + \
        "\t\tshortLabel %s\n" %(name[:17]) + \
        "\t\tparent %sviewpeaks on\n" %(accession) + \
        "\t\ttype %s\n" %(tracktype) + \
        "\t\tvisibility dense\n" + \
        "\t\tview PK\n" + \
        "\t\tpriority %d\n\n" %(n)
    n_stanzas = 1
    if not lowpass:
        lowpass = []
    if isinstance(lowpass,int):
        lowpass = [lowpass]
    extra_stanza_count = 0
    for (i, cutoff) in enumerate(lowpass,start=1):
        fn = dx.get_id()
        if not os.path.isfile(fn):
            dxpy.download_dxfile(dx.get_id(),fn)
        cutoffstr = '-lt%d' %(cutoff)
        outfn = fn + cutoffstr
        print fn, os.path.getsize(fn), subprocess.check_output('wc -l %s' %(fn), shell=True).split()[0]
        bed_fn = fn + '.bed'
        common.block_on('bigBedToBed %s %s' %(fn, bed_fn))
        common.run_pipe([
            'cat %s' %(bed_fn),
            r"""awk 'BEGIN{FS="\t";OFS="\t"}{if (($3-$2) < %d) {print $0}}'""" %(cutoff)], outfn)
        print outfn, os.path.getsize(outfn), subprocess.check_output('wc -l %s' %(outfn), shell=True).split()[0]
        if tracktype =='bigBed 6 +':
            as_file = 'narrowPeak.as'
        elif tracktype == 'bigBed 12 +':
            as_file = 'gappedPeak.as'
        else:
            print "Cannot match tracktype %s to any .as file" %(tracktype)
        bb_fn = common.bed2bb(outfn,'mm10.chrom.sizes',as_file)
        newdx = dxpy.upload_local_file(filename=bb_fn, folder="/tracks", wait_on_close=True)
        new_url, headers = newdx.get_download_url(duration=sys.maxint, preauthenticated=True)

        new_lines = [
            "\t\ttrack %s%d" %(accession,n+i),
            "\t\tbigDataUrl %s" %(new_url),
            "\t\tshortLabel %s" %(name[:17-len(cutoffstr)] + cutoffstr),
            "\t\tparent %sviewpeaks on" %(accession),
            "\t\ttype %s" %(tracktype),
            "\t\tvisibility dense",
            "\t\tview PK",
            "\t\tpriority %d\n\n" %(n+i)]
        new_stanza = '\n'.join(new_lines)
        return_string += new_stanza
        n_stanzas += 1
        os.remove(bed_fn)
        os.remove(bb_fn)
        os.remove(outfn)
        os.remove(fn)

    return(return_string, n_stanzas)

def viewsignal_stanza(accession):
    return(
        "\ttrack %sviewsignals\n" %(accession) + \
        "\tparent %s on\n" %(accession) + \
        "\tshortLabel Signals\n" + \
        "\tlongLabel Signals\n" + \
        "\tview SIG\n" + \
        "\tvisibility dense\n" + \
        "\ttype bigWig\n" + \
        "\tviewUi on\n\n")

def signal_stanza(accession, url, name, n, tracktype='bigWig'):
    return(
        "\t\ttrack %s%d\n" %(accession,n) + \
        "\t\tbigDataUrl %s\n" %(url) + \
        "\t\tshortLabel %s\n" %(name[:17]) + \
        "\t\tparent %sviewsignals on\n" %(accession) + \
        "\t\ttype %s\n" %(tracktype) + \
        "\t\tview SIG\n" + \
        "\t\tvisibility dense\n" + \
        "\t\tviewLimits 1:10\n" + \
        "\t\tmaxHeightPixels 127:64:2\n" + \
        "\t\tpriority %d\n\n" %(n))

def tf(args, analysis, experiment_accession, first_analysis):
    authid, authpw, server = processkey(args.key)
    keypair = (authid,authpw)

    stages = analysis.get('stages')
    peaks_stage =   next(stage for stage in stages if stage['execution']['name'] == "SPP Peaks")['execution']
    signals_stage = next(stage for stage in stages if stage['execution']['name'] == "ENCODE Peaks")['execution']
    idr_stage =     next(stage for stage in stages if stage['execution']['name'] == "Final IDR peak calls")['execution']

    peaks_output_names = [
        'rep1_peaks_bb',
        'rep2_peaks_bb',
        'pooled_peaks_bb'
    ]
    signals_output_names = [
        'rep1_pvalue_signal',
        'rep2_pvalue_signal',
        'pooled_pvalue_signal',
        'rep1_fc_signal',
        'rep2_fc_signal',
        'pooled_fc_signal'
    ]
    idr_output_names = [
        'conservative_set_bb',
        'optimal_set_bb'
    ]

    output_names = peaks_output_names + signals_output_names + idr_output_names
    outputs = {}
    outputs.update(dict(zip(peaks_output_names,[{'dx': dxpy.DXFile(peaks_stage['output'][output_name])} for output_name in peaks_output_names])))
    outputs.update(dict(zip(signals_output_names,[{'dx': dxpy.DXFile(signals_stage['output'][output_name])} for output_name in signals_output_names])))
    outputs.update(dict(zip(idr_output_names,[{'dx': dxpy.DXFile(idr_stage['output'][output_name])} for output_name in idr_output_names])))

    track_directory = os.path.join(args.ddir, experiment_accession)
    url_base = urlparse.urljoin(args.turl, experiment_accession+'/')
    #print "url_base %s" %(url_base)
    if not args.nodownload and not os.path.exists(track_directory):
        os.makedirs(track_directory)
    if first_analysis:
        if os.path.exists(args.tdbpath):
            if args.truncate:
                trackDb = open(args.tdbpath,'w')
            else:
                trackDb = open(args.tdbpath,'a')
        else:
            if not os.path.exists(os.path.dirname(args.tdbpath)):
                os.makedirs(os.path.dirname(args.tdbpath))
            trackDb = open(args.tdbpath, 'w')
    else:
        trackDb = open(args.tdbpath,'a')            

    for (output_name, output) in outputs.iteritems():
        local_path = os.path.join(track_directory, output['dx'].name)
        print output_name, output['dx'].get_id(), local_path
        if not args.nodownload:
            dxpy.download_dxfile(output['dx'].get_id(), local_path)
        outputs[output_name].update({'local_path' : local_path})
        #print "Joining %s and %s" %(url_base, os.path.basename(local_path))
        if args.dxf:
            url, headers = output['dx'].get_download_url(duration=sys.maxint, preauthenticated=True)
            outputs[output_name].update({'url': url})
        else:
            outputs[output_name].update({'url': urlparse.urljoin(url_base,os.path.basename(local_path))})
        #print outputs[output_name]['url']

    experiment = encoded_get(urlparse.urljoin(server,'/experiments/%s' %(experiment_accession)), keypair)
    description = '%s %s' %(
        experiment['target']['label'],
        experiment['replicates'][0]['library']['biosample']['biosample_term_name'])
    longLabel = 'E3 TF ChIP - %s - %s' %(experiment_accession, description)
    if args.tag:
        longLabel += ' - %s' %(args.tag)
    trackDb.write(composite_stanza(experiment_accession, longLabel))

    first_peaks = True
    first_signal = True
    priority = 1
    for (n, output_name) in enumerate(output_names,start=1):
        if output_name.endswith('peaks_bb') or output_name.endswith('set_bb'):
            if first_peaks:
                trackDb.write(viewpeaks_stanza(experiment_accession))
                first_peaks = False
            stanzas, n_stanzas = peaks_stanza(experiment_accession, outputs[output_name]['url'], output_name, priority, tracktype="bigBed 6 +", lowpass=args.lowpass, dx=outputs[output_name]['dx'])
            trackDb.write(stanzas)
            priority += n_stanzas
        elif output_name.endswith('_signal'):
            if first_signal:
                trackDb.write(viewsignal_stanza(experiment_accession))
                first_signal = False
            trackDb.write(signal_stanza(experiment_accession, outputs[output_name]['url'], output_name, priority, tracktype="bigWig"))
            priority += 1

    trackDb.close()


def histone(args, analysis, experiment_accession, first_analysis):
    authid, authpw, server = processkey(args.key)
    keypair = (authid,authpw)

    stages = analysis.get('stages')
    peaks_stage = next(stage for stage in stages if stage['execution']['name'] == "ENCODE Peaks")['execution']
    replicated_stages = [stage['execution'] for stage in stages if 'Final' in stage['execution']['name']]

    output_names = [
        'rep1_narrowpeaks_bb',
        'rep2_narrowpeaks_bb',
        'pooled_narrowpeaks_bb',
        'rep1_gappedpeaks_bb',
        'rep2_gappedpeaks_bb',
        'pooled_gappedpeaks_bb',
        'rep1_pvalue_signal',
        'rep2_pvalue_signal',
        'pooled_pvalue_signal',
        'rep1_fc_signal',
        'rep2_fc_signal',
        'pooled_fc_signal']

    outputs = dict(zip(output_names,[{'dx': dxpy.DXFile(peaks_stage['output'][output_name])} for output_name in output_names]))

    output_names.insert(3,'replicated_narrowpeaks_bb')
    outputs.update({'replicated_narrowpeaks_bb' : {'dx': dxpy.DXFile(next(stage['execution']['output']['overlapping_peaks_bb'] for stage in stages if stage['execution']['name'] == 'Final narrowpeaks'))}})
    output_names.insert(7,'replicated_gappedpeaks_bb')
    outputs.update({'replicated_gappedpeaks_bb' : {'dx': dxpy.DXFile(next(stage['execution']['output']['overlapping_peaks_bb'] for stage in stages if stage['execution']['name'] == 'Final gappedpeaks'))}})

    track_directory = os.path.join(args.ddir, experiment_accession)
    url_base = urlparse.urljoin(args.turl, experiment_accession+'/')
    #print "url_base %s" %(url_base)
    if not args.nodownload and not os.path.exists(track_directory):
        os.makedirs(track_directory)
    if first_analysis:
        if os.path.exists(args.tdbpath):
            if args.truncate:
                trackDb = open(args.tdbpath,'w')
            else:
                trackDb = open(args.tdbpath,'a')
        else:
            if not os.path.exists(os.path.dirname(args.tdbpath)):
                os.makedirs(os.path.dirname(args.tdbpath))
            trackDb = open(args.tdbpath, 'w')
    else:
        trackDb = open(args.tdbpath,'a')            

    for (output_name, output) in outputs.iteritems():
        local_path = os.path.join(track_directory, output['dx'].name)
        print output_name, output['dx'].get_id(), local_path
        if not args.nodownload:
            dxpy.download_dxfile(output['dx'].get_id(), local_path)
        outputs[output_name].update({'local_path' : local_path})
        #print "Joining %s and %s" %(url_base, os.path.basename(local_path))
        if args.dxf:
            url, headers = output['dx'].get_download_url(duration=sys.maxint, preauthenticated=True)
            outputs[output_name].update({'url': url})
        else:
            outputs[output_name].update({'url': urlparse.urljoin(url_base,os.path.basename(local_path))})
        #print outputs[output_name]['url']

    experiment = encoded_get(urlparse.urljoin(server,'/experiments/%s' %(experiment_accession)), keypair)
    description = '%s %s %s %s' % (
        experiment['target']['label'],
        experiment['replicates'][0]['library']['biosample']['biosample_term_name'],
        experiment['replicates'][0]['library']['biosample'].get('life_stage'),
        experiment['replicates'][0]['library']['biosample'].get('age_display'))
    longLabel = 'E3 Histone ChIP - %s - %s' %(experiment_accession, description)
    if args.tag:
        longLabel += ' - %s' %(args.tag)
    trackDb.write(composite_stanza(experiment_accession, longLabel))

    first_peaks = True
    first_signal = True
    priority = 1
    for (n, output_name) in enumerate(output_names,start=1):
        if output_name.endswith('narrowpeaks_bb'):
            if first_peaks:
                trackDb.write(viewpeaks_stanza(experiment_accession))
                first_peaks = False
            stanzas, n_stanzas = peaks_stanza(experiment_accession, outputs[output_name]['url'], output_name, priority, tracktype="bigBed 6 +", lowpass=args.lowpass, dx=outputs[output_name]['dx'])
            trackDb.write(stanzas)
            priority += n_stanzas
        elif output_name.endswith('gappedpeaks_bb'):
            if first_peaks:
                trackDb.write(viewpeaks_stanza(experiment_accession))
                first_peaks = False
            stanzas, n_stanzas = peaks_stanza(experiment_accession, outputs[output_name]['url'], output_name, priority, tracktype="bigBed 12 +", lowpass=args.lowpass, dx=outputs[output_name]['dx'])
            trackDb.write(stanzas)
            priority += n_stanzas
        elif output_name.endswith('_signal'):
            if first_signal:
                trackDb.write(viewsignal_stanza(experiment_accession))
                first_signal = False
            trackDb.write(signal_stanza(experiment_accession, outputs[output_name]['url'], output_name, priority, tracktype="bigWig"))
            priority += 1

    trackDb.close()

def main():
    args = get_args()
    first_analysis = True

    for (i, analysis_id) in enumerate(args.infile):
        analysis_id = analysis_id.strip()
        try:
            analysis = dxpy.describe(analysis_id)
        except:
            print "Invalid analysis ID %s. Skipping." % (analysis_id)
            continue

        experiment_m = re.match('^(ENCSR[0-9]{3}[A-Z]{3})', analysis['name'])
        if not experiment_m:
            print "No accession in %s, skipping." % (analysis['name'])
            continue
            # print "Temporary hack"
            # experiment_accession = "ENCSR048KZD"
        else:
            experiment_accession = experiment_m.group(1)

        if args.pipeline:
            pipeline = args.pipeline
        elif analysis['executableName'] == 'histone_chip_seq':
            pipeline = 'histone'
        elif analysis['executableName'] == 'tf_chip_seq':
            pipeline = 'tf'

        if pipeline == 'histone':
            histone(args, analysis, experiment_accession, first_analysis)
        elif pipeline == 'tf':
            tf(args, analysis, experiment_accession, first_analysis)
        else:
            print "Unrecognized pipeline: %s, skipping." % (pipeline)
            continue

        first_analysis = False

if __name__ == '__main__':
    main()
