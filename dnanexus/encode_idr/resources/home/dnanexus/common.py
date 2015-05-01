#!/usr/bin/env python

import sys, os, subprocess, shlex

def test():
	print "In common.test"

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
