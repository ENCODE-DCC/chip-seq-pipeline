#!/usr/bin/env python

import os, sys, logging, urlparse, requests, csv, StringIO, re, copy
import common

logger = logging.getLogger(__name__)

EPILOG = '''Notes:

Examples:

    %(prog)s
'''

def get_args():
    import argparse
    parser = argparse.ArgumentParser(
        description=__doc__, epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('experiments',  help='List of ENCSR accessions to report on', nargs='*', default=None)
    parser.add_argument('--infile',     help='File containing ENCSR accessions', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('--outfile',    help='tsv table of files with metadata', type=argparse.FileType('wb'), default=sys.stdout)
    parser.add_argument('--assembly',   help='Genome assembly like hg19 or mm10', required=True)
    parser.add_argument('--debug',      help="Print debug messages", default=False, action='store_true')
    parser.add_argument('--key',        help="The keypair identifier from the keyfile.", default='www')
    parser.add_argument('--keyfile',    help="The keyfile.", default=os.path.expanduser("~/keypairs.json"))

    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
    else: #use the defaulf logging level
        logging.basicConfig(format='%(levelname)s:%(message)s')

    return args

def biorep_ns(file_accession,server,keypair):
    m = re.match('^/?(files)?/?(\w*)', file_accession)
    if m:
        acc = m.group(2)
    else:
        return
    url = urlparse.urljoin(server, '/files/%s' %(acc))
    file_object = common.encoded_get(url, keypair)
    if file_object.get('derived_from'):
        for f in file_object.get('derived_from'):
            for repnum in biorep_ns(f,server,keypair):
                yield repnum
    else:
        url = urlparse.urljoin(server, '%s' %(file_object.get('replicate')))
        replicate_object = common.encoded_get(url, keypair)
        yield replicate_object.get('biological_replicate_number')

def biorep_ages(file_accession,server,keypair):
    m = re.match('^/?(files)?/?(\w*)', file_accession)
    if m:
        acc = m.group(2)
    else:
        return
    url = urlparse.urljoin(server, '/files/%s' %(acc))
    file_object = common.encoded_get(url, keypair)
    if file_object.get('derived_from'):
        for f in file_object.get('derived_from'):
            for bioage in biorep_ages(f,server,keypair):
                yield bioage
    else:
        url = urlparse.urljoin(server, '%s' %(file_object.get('replicate')))
        replicate_object = common.encoded_get(url, keypair)
        url = urlparse.urljoin(server, '%s' %(replicate_object.get('library')))
        library_object = common.encoded_get(url, keypair)
        url = urlparse.urljoin(server, '%s' %(library_object.get('biosample')))
        biosample_object = common.encoded_get(url, keypair)
        yield biosample_object.get('age_display')


def main():

    args = get_args()
    if args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    authid, authpw, server = common.processkey(args.key, args.keyfile)
    keypair = (authid,authpw)

    if args.experiments:
        exp_ids = args.experiments
    else:
        exp_ids = args.infile

    for (i, exp_id) in enumerate(exp_ids):
        exp_id = exp_id.rstrip()
        logger.info('%s' %(exp_id))
        url = urlparse.urljoin(server, 'metadata/type=experiment&accession=%s/metadata.tsv' %(exp_id))
        r = requests.get(url, auth=keypair)
        try:
            r.raise_for_status()
        except:
            logger.error('%s failed to get metadata.  GET returned %s' %(exp_id, r.return_code))
            logger.debug('%s' %(r.text))
            logger.error('Skipping ...')
            continue

        reader = csv.DictReader(StringIO.StringIO(r.text), delimiter='\t')
        fieldnames = copy.copy(reader.fieldnames)
        # fieldnames.remove('Biological replicate(s)')
        # fieldnames.insert(4,'Biological replicate(s)')
        # fieldnames.remove('Biosample Age')
        # fieldnames.insert(10,'Biosample Age')
        fieldnames.append('Derived from')
        writer = csv.DictWriter(args.outfile,fieldnames, delimiter='\t')
        writer.writeheader()
        for file_metadata in reader:
            file_accession = file_metadata.get('File accession')
            url = urlparse.urljoin(server, 'files/%s' %(file_accession))
            file_object = common.encoded_get(url, keypair)
            
            # bio_reps = sorted(list(set(biorep_ns(file_accession, server, keypair))))
            # file_metadata['Biological replicate(s)'] = ",".join([str(n) for n in bio_reps])

            # bio_ages = sorted(list(set(biorep_ages(file_accession, server, keypair)))) or ""
            # file_metadata.update({'Biosample Age': ",".join(bio_ages)})
            
            if file_object.get('derived_from'):
                derived_from = ",".join([str(f.split('/')[2]) for f in file_object.get('derived_from')])
            else:
                derived_from = None
            file_metadata.update({'Derived from': derived_from})

            #print file_metadata
            writer.writerow(file_metadata)

if __name__ == '__main__':
    main()
