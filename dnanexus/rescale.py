#!/usr/bin/env python

import os, time, common, sys
import dxpy

def rescale_scores(fn, scores_col, new_min=10, new_max=1000):
	n_peaks = common.count_lines(fn)
	sorted_fn = 'sorted-%s' %(fn)
	rescaled_fn = 'rescaled-%s' %(fn)
	out,err = common.run_pipe([
		'sort -k %dgr,%dgr %s' %(scores_col, scores_col, fn),
		r"""awk 'BEGIN{FS="\t";OFS="\t"}{if (NF != 0) print $0}'"""],
		sorted_fn)
	out, err = common.run_pipe([
		'head -n 1 %s' %(sorted_fn),
		'cut -f %s' %(scores_col)])
	max_score = float(out.strip())
	out, err = common.run_pipe([
		'tail -n 1 %s' %(sorted_fn),
		'cut -f %s' %(scores_col)])
	min_score = float(out.strip())
	out,err = common.run_pipe([
		'cat %s' %(sorted_fn),
		r"""awk 'BEGIN{OFS="\t"}{n=$%d;a=%d;b=%d;x=%d;y=%d}""" %(scores_col, min_score, max_score, new_min, new_max) + \
		r"""{$%d=int(((n-a)*(y-x)/(b-a))+x) ; print $0}'""" %(scores_col)],
		rescaled_fn)
	return rescaled_fn

def main():

	rescaled_fn = rescale_scores(sys.argv[1],5)

if __name__ == '__main__':
	main()
