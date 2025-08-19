# -*- coding: utf-8 -*-
# @Author: wjq
# @Date:   2022-02-07 17:58:46
# @Last Modified by:   wjq
# @Last Modified time: 2023-02-06 14:48:01

import glob
import os
import re
from Bio import SeqIO


def extract(name, genes, peps, dir):
    sf = open(f'{dir}/{name}.fa', 'w')
    for gene in genes:
        if gene == '*':
            continue
        seq = peps[gene]
        sf.write(f'>{gene}\n{seq}\n')


if __name__=='__main__':


    peps = {}
    files = glob.glob('cds/*.prot')
    for file in files:
        for seq_record in SeqIO.parse(file, 'fasta'):
            peps[seq_record.id] = seq_record.seq

    single_seqdir = 'Single_Copy_cds'
    for path in [single_seqdir]:
        if not os.path.isdir(path):
            os.mkdir(path)

    single = [x.strip() for x in open('SingleCopySGs.txt')]
    rows = open('sample.list.SG.pan.modified').readlines()
    for i in range(0, len(rows)):
        # lis = rows[i].split()
        lis = re.split('[\\s,]+', rows[i].strip())
        if lis[0] in single:
            extract(lis[0], lis[1:], peps, single_seqdir)



