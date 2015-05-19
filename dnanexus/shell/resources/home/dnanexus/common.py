#!/usr/bin/env python

import sys, os, subprocess, shlex, logging, re, urlparse
import dateutil.parser

def test():
	print "In common.test"

def rstrips(string, substring):
	if not string.endswith(substring):
		return string
	else:
		return string[:len(string)-len(substring)]

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

	print "In bed2bb with bed_filename=%s, chrom_sizes=%s, as_file=%s" %(bed_filename, chrom_sizes, as_file)

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

	print subprocess.check_output('ls -l', shell=True, stderr=subprocess.STDOUT)

	#this is necessary in case bedToBegBed failes to create the bb file but doesn't return a non-zero returncode
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

def processkey(key, keyfile=None):

	import json

	if not keyfile:
		try:
			keyfile = KEYFILE
		except:
			print >> sys.stderr, "No keyfile was specified."
			raise
	if key:
		keysf = open(keyfile,'r')
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

def encoded_get(url, keypair=None, frame='object'):
	import urlparse, requests
	HEADERS = {'content-type': 'application/json'}
	url = urlparse.urljoin(url,'?format=json&frame=%s' %(frame))
	if keypair:
		response = requests.get(url, auth=keypair, headers=HEADERS)
	else:
		response = requests.get(url, headers=HEADERS)
	return response.json()

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
		if not subprocess.call('which md5', shell=True):
			md5_command = 'md5 -q'
		elif not subprocess.call('which md5sum', shell=True):
			md5_command = 'md5sum'
		else:
			md5_command = ''
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
		logger.error("%s Cannot compare bam date %s with fastq date %s" %(e, bam.get('date_created'), f.get('date_created')))
		raise
	
	return result

def biorep_ns_generator(file_accession,server,keypair):
	m = re.match('^/?(files)?/?(\w*)', file_accession)
	if m:
		acc = m.group(2)
	else:
		return
	url = urlparse.urljoin(server, '/files/%s' %(acc))
	file_object = encoded_get(url, keypair)
	if file_object.get('derived_from'):
		for f in file_object.get('derived_from'):
			for repnum in biorep_ns_generator(f,server,keypair):
				yield repnum
	else:
		url = urlparse.urljoin(server, '%s' %(file_object.get('replicate')))
		replicate_object = encoded_get(url, keypair)
		yield replicate_object.get('biological_replicate_number')

def biorep_ns(file_accession,server,keypair):
	return list(set(biorep_ns_generator(file_accession,server,keypair)))


