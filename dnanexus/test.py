#!/usr/bin/env python

import pdb
import os.path, sys, subprocess, logging, re, json, urlparse, requests
import dxpy

project = dxpy.find_one_project(name='input_mapping', return_handler=True)
test_applet = dxpy.find_one_data_object(classname='applet', name='test', project=project.get_id(), return_handler=True)
test2_applet = dxpy.find_one_data_object(classname='applet', name='test2', project=project.get_id(), return_handler=True)

def main():
	wf = dxpy.new_dxworkflow(
		title="Test",
		name="test",
		project=project.get_id(),
		folder='/testing')

	stage1 = wf.add_stage(
		test_applet,
		name="test",
		folder="/testing"
	)

	stage2 = wf.add_stage(
		test2_applet,
		name="test2",
		folder='/testing',
		stage_input={'inhash': dxpy.dxlink({'stage': stage1, 'outputField': 'testout'})}
	)

	wf.run({})

if __name__ == '__main__':
	main()
