#!python

import GffFile
import re
import argparse
import MyUtil

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='exclude by gene list')
    parser.add_argument('-f', '--gff', dest='gff_file', type=str, help="The gff file", required=True)
    parser.add_argument('-g', '--geneList', dest='gene_list', type=str, help="The gene list to be extracted", required=True)
    args = parser.parse_args()

    excludeGeneIdsSet = set()
    with open(args.gene_list) as f:
        for line in f:
            line = line.strip()
            if not line.startswith("#"):
                excludeGeneIdsSet.add(line)
            #print(line)
    f.close()

    geneIdsSet=set()
    chromosome_gene_dict, chromosome_gene_list, transcript_gene_map, gene_chr_map = GffFile.readGff(args.gff_file)
    for chr in chromosome_gene_dict:
        for gene in chromosome_gene_dict[chr]:
            if gene not in excludeGeneIdsSet:
                geneIdsSet.add(gene)

    MyUtil.extractByGeneids(args.gff_file, geneIdsSet)

