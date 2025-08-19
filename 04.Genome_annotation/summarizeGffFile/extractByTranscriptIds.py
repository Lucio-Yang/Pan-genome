#!python

import GffFile
import re
import argparse
import MyUtil

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='extract by trqanscript list')
    parser.add_argument('-f', '--gff', dest='gff_file', type=str, help="The gff file", required=True)
    parser.add_argument('-t', '--transcriptList', dest='transcript_list', type=str, help="The gene list to be extracted", required=True)
    args = parser.parse_args()

    transcriptIdsSet = set()
    with open(args.transcript_list) as f:
        for line in f:
            line = line.strip()
            if not line.startswith("#"):
                transcriptIdsSet.add(line)
            #print(line)
    f.close()
    MyUtil.extractByTranscripts(args.gff_file, transcriptIdsSet)

