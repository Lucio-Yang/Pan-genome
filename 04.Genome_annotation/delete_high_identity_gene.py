#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to remove shorter sequences in alignments when their identity exceeds a specified threshold, then extract corresponding annotation information from GFF file for the remaining sequences.
Note: Alignments are generated using blastp with the parameters: `-outfmt 6 -evalue 1e-5`
"""
__author__ = "Bin Chen"
__date__ = "2025-2-25"
__version__ = "1.0"

import argparse
import sys
from collections import defaultdict
import re

def gff_parser(gff_file):
    '''
    Read a GFF3 file and calculate the length of the coding region (CDS) for each mRNA
    '''
    cds_lengths = defaultdict(int)  # Store CDS lengths for each mRNA
    mRNA_to_gene = {}  # Map mRNA IDs to gene IDs

    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            
            feature_type = fields[2]
            start, end = int(fields[3]), int(fields[4])
            attributes = {attr.split("=")[0]: attr.split("=")[1] for attr in fields[8].split(";") if "=" in attr}

            if feature_type == "mRNA":
                mRNA_to_gene[attributes.get("ID", "")] = attributes.get("Parent", "")

            elif feature_type == "CDS":
                parent_id = attributes.get("Parent", "")
                cds_lengths[parent_id] += (end - start + 1)

    #print(f'mRNA_to_gene: {mRNA_to_gene}')
    #print(f'cds_lengths: {cds_lengths}')
    return mRNA_to_gene, cds_lengths


def get_redundant_seqid(blastp_m8, identity_threshold):
    '''
    Get shorter sequence ID in an alignment when their identity exceeds a specified threshold
    m8 header: ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    '''
    with open(blastp_m8, 'r') as m8:
        redundant_seqids = set()
        for line in m8:
            fields = line.strip().split("\t")   
            if fields[0] == fields[1]:
                continue
            if float(fields[2]) < identity_threshold*100:
                continue
            
            if cds_lengths[fields[0]] > cds_lengths[fields[1]]:
                redundant_seqids.add(fields[1])
            else:
                redundant_seqids.add(fields[0])
        
    #print(f'redundant_seqids: {redundant_seqids}')
    return redundant_seqids


def filter_gff(gff_file, outfile):
    '''
    Remove annotation information for the redundant sequences from the GFF3 file.
    '''
    redundant_seqgids = {mRNA_to_gene[k] for k in mRNA_to_gene if k in redundant_seqids}  # Get redundant Gene IDs
    ID_pat = re.compile(r'ID=([^;]+)')
    Parent_pat = re.compile(r'Parent=([^;]+)')
    
    with open(gff_file, 'r') as f, open(outfile, 'w') as out:
        for line in f:
            if line.startswith("#"):
                out.write(line)
                continue

            fields = line.strip().split("\t")
            feature, attributes = fields[2], fields[8]
            if feature == "gene":
                gene_ID = ID_pat.search(attributes).group(1)
                if gene_ID in redundant_seqgids:
                    continue
            
            if ID_pat.search(attributes).group(1) in redundant_seqids or Parent_pat.search(attributes).group(1) in redundant_seqids:
                continue

            out.write(line)

    print(f"Filtered GFF3 saved as {out}", flush=True)
    

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Remove shorter sequences in alignments when their identity exceeds a specified threshold, then extract corresponding annotation information from GFF file for the remaining sequences.")
        print("Usage: python3 delete_high_identity_gene.py <gff_file> <blastp_m8> <identity_threshold> <outfile>")
        print("Note: identity_threshold must be of type FLOAT")
        sys.exit(1)

    gff_file = sys.argv[1]
    blastp_m8 = sys.argv[2]
    identity_threshold = float(sys.argv[3])
    outfile = sys.argv[4]

    mRNA_to_gene, cds_lengths = gff_parser(gff_file)
    redundant_seqids = get_redundant_seqid(blastp_m8, identity_threshold)
    filter_gff(gff_file, outfile)
