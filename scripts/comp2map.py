#!/usr/bin/python

import argparse

parser = argparse.ArgumentParser(description='Convert NanoSim compatibility TSV to read-to-model map')
parser.add_argument('fpath', metavar='TSV', help='compatibility file path')

args = parser.parse_args()

print('read_id', 'transcript_id', sep='\t')
with open(args.fpath) as fh:
    for line in fh:
        cols = line.split('\t')
        if len(cols) >= 5:
            readname = cols[0]
            if len(cols[2]) > 0:
                for t in cols[2].split(','):
                    print(readname, t, sep='\t')
