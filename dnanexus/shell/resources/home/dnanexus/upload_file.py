#!/usr/bin/env python
# -*- coding: latin-1 -*-

import os, sys, time, subprocess, json, requests

HEADERS = {
    'Content-type': 'application/json',
    'Accept': 'application/json',
}

path = 'test.fastq'
FILE_URL = 'http://test.encodedcc.org/TSTFF867178/upload/'
ENCODED_KEY = '...'
ENCODED_SECRET_KEY = '...'

response = requests.get(FILE_URL, headers=HEADERS, auth=(ENCODED_KEY, ENCODED_SECRET_KEY))
try:
    response.raise_for_status()
except:
    print('File object GET failed')
    raise
item = response.json()['@graph'][0]
print(json.dumps(item, indent=4, sort_keys=True))

creds = item['upload_credentials']
env = os.environ.copy()
env.update({
    'AWS_ACCESS_KEY_ID': creds['access_key'],
    'AWS_SECRET_ACCESS_KEY': creds['secret_key'],
    'AWS_SECURITY_TOKEN': creds['session_token'],
})

# ~10s/GB from Stanford - AWS Oregon
# ~12-15s/GB from AWS Ireland - AWS Oregon
print("Uploading file.")
start = time.time()
try:
    subprocess.check_call(['aws', 's3', 'cp', path, creds['upload_url']], env=env)
except subprocess.CalledProcessError as e:
    # The aws command returns a non-zero exit code on error.
    print("Upload failed with exit code %d" % e.returncode)
    sys.exit(e.returncode)
else:
    end = time.time()
    duration = end - start
    print("Uploaded in %.2f seconds" % duration)
