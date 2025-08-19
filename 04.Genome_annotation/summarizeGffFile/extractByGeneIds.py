#!python

import GffFile
import re
import argparse
import MyUtil

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='extract by gene list')
    parser.add_argument('-f', '--gff', dest='gff_file', type=str, help="The gff file", required=True)
    parser.add_argument('-g', '--geneList', dest='gene_list', type=str, help="The gene list to be extracted", required=True)
    args = parser.parse_args()

    geneIdsSet = set()
    with open(args.gene_list) as f:
        for line in f:
            line = line.strip()
            if not line.startswith("#"):
                geneIdsSet.add(line)
            #print(line)
    f.close()
    MyUtil.extractByGeneids(args.gff_file, geneIdsSet)

