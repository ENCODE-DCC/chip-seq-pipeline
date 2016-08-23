#!/usr/bin/env python2

import sys
import os
import subprocess
import shlex
import logging
import re
import urlparse
import dateutil.parser
from time import time, sleep

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
logger.propagate = True


def test():
    print "In common.test"


def flat(l):
    result = []
    for el in l:
        if hasattr(el, "__iter__") and not isinstance(el, basestring):
            result.extend(flat(el))
        else:
            result.append(el)
    return result


def rstrips(string, substring):
    if not string.endswith(substring):
        return string
    else:
        return string[:len(string)-len(substring)]

def touch(fname, times=None):
    with open(fname, 'a'):
        os.utime(fname, times)

def block_on(command):
    process = subprocess.Popen(shlex.split(command), stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
    for line in iter(process.stdout.readline, ''):
        sys.stdout.write(line)
    process.wait()
    return process.returncode

def run_pipe(steps, outfile=None):
    #break this out into a recursive function
    #TODO:  capture stderr
    from subprocess import Popen, PIPE
    p = None
    p_next = None
    first_step_n = 1
    last_step_n = len(steps)
    for n,step in enumerate(steps, start=first_step_n):
        print "step %d: %s" %(n,step)
        if n == first_step_n:
            if n == last_step_n and outfile: #one-step pipeline with outfile
                with open(outfile, 'w') as fh:
                    print "one step shlex: %s to file: %s" %(shlex.split(step), outfile)
                    p = Popen(shlex.split(step), stdout=fh)
                break
            print "first step shlex to stdout: %s" %(shlex.split(step))
            p = Popen(shlex.split(step), stdout=PIPE)
            #need to close p.stdout here?
        elif n == last_step_n and outfile: #only treat the last step specially if you're sending stdout to a file
            with open(outfile, 'w') as fh:
                print "last step shlex: %s to file: %s" %(shlex.split(step), outfile)
                p_last = Popen(shlex.split(step), stdin=p.stdout, stdout=fh)
                p.stdout.close()
                p = p_last
        else: #handles intermediate steps and, in the case of a pipe to stdout, the last step
            print "intermediate step %d shlex to stdout: %s" %(n,shlex.split(step))
            p_next = Popen(shlex.split(step), stdin=p.stdout, stdout=PIPE)
            p.stdout.close()
            p = p_next
    out,err = p.communicate()
    return out,err

def uncompress(filename):
    #leaves compressed file intact
    m = re.match('(.*)(\.((gz)|(Z)|(bz)|(bz2)))',filename)
    if m:
        basename = m.group(1)
        logging.info(subprocess.check_output(shlex.split('ls -l %s' %(filename))))
        logging.info("Decompressing %s" %(filename))
        #logging.info(subprocess.check_output(shlex.split('gzip -dc %s' %(filename))))
        out,err = run_pipe([
            'gzip -dc %s' %(filename)],
            basename)
        logging.info(subprocess.check_output(shlex.split('ls -l %s' %(basename))))
        return basename
    else:
        return filename

def compress(filename):
    #leaves uncompressed file intact
    if re.match('(.*)(\.((gz)|(Z)|(bz)|(bz2)))',filename):
        return filename
    else:
        logging.info(subprocess.check_output(shlex.split('cp %s tmp' %(filename))))
        logging.info(subprocess.check_output(shlex.split('ls -l %s' %(filename))))
        logging.info("Compressing %s" %(filename))
        logging.info(subprocess.check_output(shlex.split('gzip %s' %(filename))))
        new_filename = filename + '.gz'
        logging.info(subprocess.check_output(shlex.split('cp tmp %s' %(filename))))
        logging.info(subprocess.check_output(shlex.split('ls -l %s' %(new_filename))))
        return new_filename

def count_lines(fname):
    wc_output = subprocess.check_output(shlex.split('wc -l %s' %(fname)))
    lines = wc_output.split()[0]
    return int(lines)

def bed2bb(bed_filename, chrom_sizes, as_file, bed_type='bed6+4'):
    if bed_filename.endswith('.bed'):
        bb_filename = bed_filename[:-4] + '.bb'
    else:
        bb_filename = bed_filename + '.bb'
    bed_filename_sorted = bed_filename + ".sorted"

    logging.debug("In bed2bb with bed_filename=%s, chrom_sizes=%s, as_file=%s" %(bed_filename, chrom_sizes, as_file))

    print "Sorting"
    print subprocess.check_output(shlex.split("sort -k1,1 -k2,2n -o %s %s" %(bed_filename_sorted, bed_filename)), shell=False, stderr=subprocess.STDOUT)

    for fn in [bed_filename, bed_filename_sorted, chrom_sizes, as_file]:
        print "head %s" %(fn)
        print subprocess.check_output('head %s' %(fn), shell=True, stderr=subprocess.STDOUT)

    command = "bedToBigBed -type=%s -as=%s %s %s %s" %(bed_type, as_file, bed_filename_sorted, chrom_sizes, bb_filename)
    print command
    try:
        process = subprocess.Popen(shlex.split(command), stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
        for line in iter(process.stdout.readline, ''):
            sys.stdout.write(line)
        process.wait()
        returncode = process.returncode
        if returncode != 0:
            raise subprocess.CalledProcessError
    except:
        e = sys.exc_info()[0]
        sys.stderr.write('%s: bedToBigBed failed. Skipping bb creation.' %(e))
        return None

    #print subprocess.check_output('ls -l', shell=True, stderr=subprocess.STDOUT)

    #this is necessary in case bedToBegBed failes to create the bb file but doesn't return a non-zero returncode
    try:
        os.remove(bed_filename_sorted)
    except:
        pass
    if not os.path.isfile(bb_filename):
        bb_filename = None

    print "Returning bb file %s" %(bb_filename)
    return bb_filename

def rescale_scores(fn, scores_col, new_min=10, new_max=1000):
    n_peaks = count_lines(fn)
    sorted_fn = '%s-sorted' %(fn)
    rescaled_fn = '%s-rescaled' %(fn)
    out,err = run_pipe([
        'sort -k %dgr,%dgr %s' %(scores_col, scores_col, fn),
        r"""awk 'BEGIN{FS="\t";OFS="\t"}{if (NF != 0) print $0}'"""],
        sorted_fn)
    out, err = run_pipe([
        'head -n 1 %s' %(sorted_fn),
        'cut -f %s' %(scores_col)])
    max_score = float(out.strip())
    out, err = run_pipe([
        'tail -n 1 %s' %(sorted_fn),
        'cut -f %s' %(scores_col)])
    min_score = float(out.strip())
    out,err = run_pipe([
        'cat %s' %(sorted_fn),
        r"""awk 'BEGIN{OFS="\t"}{n=$%d;a=%d;b=%d;x=%d;y=%d}""" %(scores_col, min_score, max_score, new_min, new_max) + \
        r"""{$%d=int(((n-a)*(y-x)/(b-a))+x) ; print $0}'""" %(scores_col)],
        rescaled_fn)
    return rescaled_fn


def slop_clip(filename, chrom_sizes, bed_type='bed'):
    assert bed_type in ['bed', 'gappedPeak'], \
        "slop_clip: unsupported bed_type %s" % (bed_type)
    clipped_fn = '%s-clipped' % (filename)
    # Remove coordinates outside chromosome sizes
    pipe = ['slopBed -i %s -g %s -b 0' % (filename, chrom_sizes),
            'bedClip stdin %s %s' % (chrom_sizes, clipped_fn)]
    print pipe
    out, err = run_pipe(pipe)
    if bed_type == 'bed':
        return clipped_fn
    elif bed_type == 'gappedPeak':
        clipped_gappedPeaks_fn = '%s-gapclipped' % (clipped_fn)
        import csv
        import copy
        with open(clipped_fn, 'rb') as in_fh, open(clipped_gappedPeaks_fn, 'wb') as out_fh:
            fieldnames = [
                'chrom', 'chromStart', 'chromEnd', 'name', 'score',
                'strand', 'thickStart', 'thickEnd', 'reserved',
                'blockCount', 'blockSizes', 'blockStarts',
                'signalValue', 'pValue', 'qValue']
            reader = \
                csv.DictReader(in_fh, fieldnames=fieldnames, delimiter='\t')
            writer = \
                csv.DictWriter(out_fh, fieldnames=fieldnames, delimiter='\t', lineterminator='\n')
            for line in reader:
                peak = dict(zip(fieldnames, [
                    line['chrom'],
                    int(line['chromStart']),
                    int(line['chromEnd']),
                    line['name'],
                    line['score'],
                    line['strand'],
                    int(line['thickStart']),
                    int(line['thickEnd']),
                    line['reserved'],
                    int(line['blockCount']),
                    [int(s) for s in line['blockSizes'].split(',')],
                    [int(s) for s in line['blockStarts'].split(',')],
                    line['signalValue'],
                    line['pValue'],
                    line['qValue']]))
                chromStart = peak['chromStart']
                chromEnd = peak['chromEnd']
                if peak['thickStart'] < chromStart:
                    peak['thickStart'] = chromStart
                if peak['thickEnd'] > chromEnd:
                    peak['thickEnd'] = chromEnd
                blocks = [
                    dict(zip(
                        ['start', 'end'],
                        [chromStart+peak['blockStarts'][i], chromStart+(peak['blockStarts'][i] + peak['blockSizes'][i])]))
                    for i in range(peak['blockCount'])]
                newblocks = []
                for block in blocks:
                    if block['start'] < chromStart and block['end'] < chromStart:
                        continue
                    elif block['start'] > chromEnd and block['end'] > chromEnd:
                        continue
                    elif block['end'] > chromEnd:
                        block['end'] = chromEnd
                    elif block['start'] < block['start']:
                        block['start'] = chromStart
                    newblocks.append(block)
                if not [b for b in blocks if b['start'] == chromStart]:
                    newblocks.insert(
                        0,
                        {'start': chromStart, 'end': chromStart+1})
                if not [b for b in blocks if b['end'] == chromEnd]:
                    newblocks.append(
                        {'start': chromEnd-1, 'end': chromEnd})

                peak['blockCount'] = len(newblocks)
                peak['blockSizes'] = \
                    [(block['end']-block['start']) for block in newblocks]
                peak['blockStarts'] = \
                    [block['start']-chromStart for block in newblocks]
                if peak['blockCount']:
                    peak['blockSizes'] = \
                        ','.join([str(x) for x in peak['blockSizes']])
                    peak['blockStarts'] = \
                        ','.join([str(x) for x in peak['blockStarts']])
                    writer.writerow(peak)

        return clipped_gappedPeaks_fn


def processkey(key=None, keyfile=None):

    import json

    if not (key or keyfile) and os.getenv('ENCODE_AUTHID',None) and os.getenv('ENCODE_AUTHPW',None) and os.getenv('ENCODE_SERVER',None):
        authid = os.getenv('ENCODE_AUTHID',None)
        authpw = os.getenv('ENCODE_AUTHPW',None)
        server = os.getenv('ENCODE_SERVER',None)
    else:
        if not keyfile:
            if 'KEYFILE' in globals(): #this is to support scripts where KEYFILE is a global
                keyfile = KEYFILE
            elif os.path.isfile(os.path.expanduser('~/keypairs.json')):
                keyfile = os.path.expanduser('~/keypairs.json')
            else:
                logging.error("Keyfile must be specified, in ~/keypairs.json or in global KEYFILE.")
                return None
        if key:
            try:
                keysf = open(keyfile,'r')
            except IOError as e:
                logging.error("Failed to open keyfile %s" %(keyfile))
                logging.error("e.")
                return None
            except:
                raise
            keys_json_string = keysf.read()
            keysf.close()
            try:
                keys = json.loads(keys_json_string)
            except ValueError as e:
                logging.error(e.message)
                logging.error("Keyfile %s not in parseable JSON" %(keyfile))
                return None
            except:
                raise
            try:
                key_dict = keys[key]
            except ValueError:
                logging.error(e.message)
                logging.error("Keyfile %s has no key named %s" %(keyfile,key))
                return None
            except:
                raise
        else:
            key_dict = {}

        if key_dict:
            authid = key_dict.get('key')
            authpw = key_dict.get('secret')
            server = key_dict.get('server')
        else:
            return None

    if not server.endswith("/"):
        server += "/"

    return (authid,authpw,server)


def encoded_get(url, keypair=None, frame='object', return_response=False):
    import urlparse
    import urllib
    import requests
    #it is not strictly necessary to include both the accept header, and format=json, but we do
    #so as to get exactly the same URL as one would use in a web browser

    RETRY_CODES = [500]
    RETRY_EXCEPTIONS = (requests.exceptions.ConnectionError, requests.exceptions.SSLError)
    HEADERS = {'accept': 'application/json'}

    url_obj = urlparse.urlsplit(url)
    new_url_list = list(url_obj)
    query = urlparse.parse_qs(url_obj.query)
    if 'format' not in query:
        new_url_list[3] += "&format=json"
    if 'frame' not in query:
        new_url_list[3] += "&frame=%s" %(frame)
    if 'limit' not in query:
        new_url_list[3] += "&limit=all"
    if new_url_list[3].startswith('&'):
        new_url_list[3] = new_url_list[3].replace('&','',1)
    get_url = urlparse.urlunsplit(new_url_list)
    logging.debug('encoded_get: %s' %(get_url))
    max_retries = 10
    max_sleep = 10
    while max_retries:
        try:
            if keypair:
                response = requests.get(get_url, auth=keypair, headers=HEADERS)
            else:
                response = requests.get(get_url, headers=HEADERS)
        except RETRY_EXCEPTIONS as e:
            logging.warning(
                "%s ... %d retries left." % (e, max_retries))
            sleep(max_sleep - max_retries)
            max_retries -= 1
            continue
        except Exception as e:
            logging.error("%s" % (e))
            if return_response:
                return response
            else:
                return None
        else:
            if response.status_code in RETRY_CODES:
                logging.warning(
                    "%d %s ... %d retries left."
                    % (response.status_code, response.text, max_retries))
                sleep(max_sleep - max_retries)
                max_retries -= 1
                continue
            if return_response:
                return response
            else:
                try:
                    return response.json()
                except:
                    return response.text
    logging.error("Max retries exhausted.")
    if return_response:
        return response
    else:
        return None


def encoded_update(method, url, keypair, payload, return_response):
    import urlparse, urllib, requests, json
    if method == 'patch':
        request_method = requests.patch
    elif method == 'post':
        request_method = requests.post
    elif method == 'put':
        request_method = requests.put
    else:
        logging.error('Invalid HTTP method: %s' %(method))
        return

    RETRY_CODES = [500]
    RETRY_EXCEPTIONS = (requests.exceptions.ConnectionError, requests.exceptions.SSLError)
    HEADERS = {'accept': 'application/json', 'content-type': 'application/json'}

    max_retries = 10
    max_sleep = 10
    while max_retries:
        try:
            response = request_method(
                url, auth=keypair, headers=HEADERS, data=json.dumps(payload))
        except RETRY_EXCEPTIONS as e:
            logging.warning(
                "%s ... %d retries left." % (e, max_retries))
            sleep(max_sleep - max_retries)
            max_retries -= 1
            continue
        except Exception as e:
            logging.error("%s" % (e))
            return None
        else:
            if response.status_code in RETRY_CODES:
                logging.warning(
                    "%d %s ... %d retries left."
                    % (response.status_code, response.text, max_retries))
                sleep(max_sleep - max_retries)
                max_retries -= 1
                continue
            if return_response:
                return response
            else:
                try:
                    return response.json()
                except:
                    return response.text
    logging.error("Max retries exhausted.")
    if return_response:
        return response
    else:
        return None

def encoded_patch(url, keypair, payload, return_response=False):
    return encoded_update('patch', url, keypair, payload, return_response)

def encoded_post(url, keypair, payload, return_response=False):
    return encoded_update('post', url, keypair, payload, return_response)

def encoded_put(url, keypair, payload, return_response=False):
    return encoded_update('put', url, keypair, payload, return_response)

def pprint_json(JSON_obj):
    import json
    print json.dumps(JSON_obj, sort_keys=True, indent=4, separators=(',', ': '))

def merge_dicts(*dict_args):
    '''
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    '''
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result

def md5(fn):
    if 'md5_command' not in globals():
        global md5_command
        try:
            subprocess.check_call('which md5', shell=True)
        except:
            try:
                subprocess.check_call('which md5sum', shell=True)
            except:
                md5_command = None
            else:
                md5_command = 'md5sum'
        else:
            md5_command = 'md5 -q'

    md5_output = subprocess.check_output(' '.join([md5_command, fn]), shell=True)
    return md5_output.partition(' ')[0].rstrip()

def after(date1, date2):
    try:
        result = dateutil.parser.parse(date1) > dateutil.parser.parse(date2)
    except TypeError:
        if not re.search('\+.*$', date1):
            date1 += 'T00:00:00-07:00'
        if not re.search('\+.*$', date2):
            date1 += 'T00:00:00-07:00'
    try:
        result = dateutil.parser.parse(date1) > dateutil.parser.parse(date2)
    except Exception as e:
        logger.error("%s Cannot compare %s with %s" %(e, date1, date2))
        raise
    else:
        return result


def biorep_ns_generator(f, server, keypair):
    if isinstance(f, dict):
        acc = f.get('accession')
    else:
        m = re.match('^/?(files)?/?(\w*)', f)
        if m:
            acc = m.group(2)
        else:
            acc = re.search('ENCFF[0-9]{3}[A-Z]{3}', f).group(0)
    if not acc:
        return
    url = urlparse.urljoin(server, '/files/%s' % (acc))
    file_object = encoded_get(url, keypair)
    if file_object.get('derived_from'):
        for derived_from in file_object.get('derived_from'):
            for repnum in biorep_ns_generator(derived_from, server, keypair):
                yield repnum
    else:
        url = urlparse.urljoin(server, '%s' % (file_object.get('replicate')))
        replicate_object = encoded_get(url, keypair)
        yield replicate_object.get('biological_replicate_number')


def biorep_ns(f, server, keypair):
    return [n for n in set(biorep_ns_generator(f, server, keypair)) if n is not None]


def derived_from_references_generator(f, server, keypair):
    if isinstance(f, dict):
        acc = f.get('accession')
    else:
        m = re.match('^/?(files)?/?(\w*)', f)
        if m:
            acc = m.group(2)
        else:
            acc = re.search('ENCFF[0-9]{3}[A-Z]{3}', f).group(0)
    if not acc:
        return
    url = urlparse.urljoin(server, '/files/%s' % (acc))
    file_object = encoded_get(url, keypair)

    if not file_object.get('derived_from'):
        return
    else:
        for derived_from_uri in file_object.get('derived_from', []):
            derived_from_url = urlparse.urljoin(server, derived_from_uri)
            derived_from_file = encoded_get(derived_from_url, keypair)
            if derived_from_file.get('output_category') == "reference":
                yield derived_from_file.get('@id')
            else:
                for derived_from_reference in derived_from_references_generator(derived_from_file, server, keypair):
                    yield derived_from_reference


def derived_from_references(f, server, keypair):
    return [n for n in set(derived_from_references_generator(f, server, keypair)) if n is not None]


def s3_cp(file_object, local_fname, server, keypair, overwrite=False):
    # TODO check overwrite and regenerate credential if necessary
    if 'upload_credentials' not in file_object:
        url = server + '/files/%s/upload/' % (file_object['accession'])
        response_json = encoded_get(url, keypair)
        logger.debug('s3_cp: Got %s response_json %s' % (url, response_json))
        creds = response_json['@graph'][0]['upload_credentials']
    else:
        creds = file_object['upload_credentials']
    logger.debug('s3_cp: Got creds %s' % (creds))
    env = os.environ.copy()
    env.update({
        'AWS_ACCESS_KEY_ID': creds['access_key'],
        'AWS_SECRET_ACCESS_KEY': creds['secret_key'],
        'AWS_SECURITY_TOKEN': creds['session_token'],
    })
    logger.info("Uploading file.")
    start = time()
    logger.debug('accession_file local_fname %s' % (local_fname))
    logger.debug('accession_file upload_url %s' % (creds['upload_url']))
    try:
        return_code = subprocess.check_call(
            ['aws', 's3', 'cp', local_fname, creds['upload_url'], '--quiet'],
            env=env)
    except subprocess.CalledProcessError as e:
        # The aws command returns a non-zero exit code on error.
        logger.error("Upload failed with exit code %d" % e.returncode)
        return e.returncode
    else:
        end = time()
        duration = end - start
        logger.info("Uploaded in %.2f seconds" % duration)
        return return_code
