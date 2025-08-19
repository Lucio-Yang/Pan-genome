#!python
import numpy as np

import FastaFile
import GffFile
import re
# song@mpipz.mpg.de


def overlap_with_certain_gene(position, chromosome_name, strand, chromosome_gene_dict):
    for gene_name in chromosome_gene_dict[chromosome_name]:
        if chromosome_gene_dict[chromosome_name][gene_name].start > position:
            return None
        if chromosome_gene_dict[chromosome_name][gene_name].end < position:
            continue
        if (chromosome_gene_dict[chromosome_name][gene_name].start < position) & (chromosome_gene_dict[chromosome_name][gene_name].end > position) & (chromosome_gene_dict[chromosome_name][gene_name].strand == strand):
            return gene_name
    return None


def overlap_certain_gene_with_genes(gene, chromosome_name, chromosome_gene_dict):
    overlapped_genes = set()
    for gene_name in chromosome_gene_dict[chromosome_name]:
        # print(gene_name)
        # print(chromosome_gene_dict[chromosome_name][gene_name].start)
        # print(chromosome_gene_dict[chromosome_name][gene_name].end)
        # print(chromosome_gene_dict[chromosome_name][gene_name].strand)
        if chromosome_gene_dict[chromosome_name][gene_name].start > gene.end:
            continue
        if chromosome_gene_dict[chromosome_name][gene_name].end < gene.start:
            continue
        if (chromosome_gene_dict[chromosome_name][gene_name].start <= gene.start) and (gene.start <= chromosome_gene_dict[chromosome_name][gene_name].end ) and (chromosome_gene_dict[chromosome_name][gene_name].strand == gene.strand):
            overlapped_genes.add(gene_name)
            continue
        if (chromosome_gene_dict[chromosome_name][gene_name].start <= gene.end) and (gene.end <= chromosome_gene_dict[chromosome_name][gene_name].end) and (chromosome_gene_dict[chromosome_name][gene_name].strand == gene.strand):
            overlapped_genes.add(gene_name)
            continue
    return overlapped_genes


def extractByGeneidsAndTranscripts(gffFile, geneIdsSet, transcriptiDSet):
    with open(gffFile) as f:
        for line in f:
            line = line.strip()
            m = re.search(r'^(\S+)\t(\S+)\t\S+\t(\d+)\t+(\d+)\t+(\S+)\t+(\S+)\t+(\S+)\t+.*Parent=([a-zA-Z0-9._^+-]+)', line)
            if( m != None ):
                transcript_name = m.group(8)
                if transcript_name in transcriptiDSet:
                    print(line)
            m1 = re.search(r'^(\S+)\t(\S+)\tmRNA\t(\d+)\t+(\d+)\t+(\S+)\t+(\S+)\t+(\S+)\t+.*ID=([a-zA-Z0-9._^+-]+);', line)
            if (m1 != None):
                transcript_name = m1.group(8)
                if transcript_name in transcriptiDSet:
                    print(line)
            else:
                m2 = re.search(r'^(\S+)\t(\S+)\tgene\t(\d+)\t+(\d+)\t+(\S+)\t+(\S+)\t+(\S+)\t+.*ID=([a-zA-Z0-9._^+-]+)$', line)
                if (m2 != None):
                    gene_name = m2.group(8)
                    if gene_name in geneIdsSet:
                        print(line)
                else:
                    m3 = re.search(r'^(\S+)\t(\S+)\tgene\t(\d+)\t+(\d+)\t+(\S+)\t+(\S+)\t+(\S+)\t+.*ID=([a-zA-Z0-9._^+-]+);', line)
                    if (m3 != None):
                        gene_name = m3.group(8)
                        if gene_name in geneIdsSet:
                            print(line)


def extractByGeneids(gffFile, geneIdsSet):
    chromosome_gene_dict, chromosome_gene_list, transcript_gene_map, gene_chr_map = GffFile.readGff( gffFile )
    transcriptiDSet = set()
    for chr in chromosome_gene_dict:
        for geneName in chromosome_gene_dict[chr]:
            if geneName in geneIdsSet:
                for transcript in chromosome_gene_dict[chr][geneName].transcripts:
                    transcriptiDSet.add(transcript.name)
#                    print(transcript.name)
    extractByGeneidsAndTranscripts(gffFile, geneIdsSet, transcriptiDSet)

def extractByTranscripts(gffFile, transcriptiDSet):
    chromosome_gene_dict, chromosome_gene_list, transcript_gene_map, gene_chr_map = GffFile.readGff(gffFile)
    geneIdsSet = set()
    for transcript in transcriptiDSet:
        if transcript in transcript_gene_map:
            geneIdsSet.add(transcript_gene_map[transcript])

    extractByGeneids(gffFile, geneIdsSet)
