# -*- coding: utf-8 -*-


import glob
import os
import re
from Bio import SeqIO
import subprocess
from multiprocessing import Pool

def run_iqtree(name):
    subprocess.call(f'mafft --auto Single_Copy_cds/{name}.fa  > mafft_out/{name}.mafft.fa', shell=True)
    # subprocess.call(f'FastTree mafft_out/{name}.mafft.fa > Single_tree/{name}_AA.nwk', shell=True)   
    #subprocess.call(f'trimal -in mafft_out/{name}.mafft.fa -out trimal_aln/{name}.mafft.trimal.aln -gt 0.3', shell=True)
    # subprocess.call(f'pal2nal.pl mafft_out/{name}.mafft.fa Single_Copy_cds/{name}.fa -output fasta > mafft_out/{name}.cds.aln', shell=True)
    # subprocess.call(f'pxclsq -p 0.2 -s mafft_out/{name}.cds.aln -o mafft_out/{name}.cds.aln-0.2cln', shell=True)
    subprocess.call(f'iqtree2 -s mafft_out/{name}.mafft.fa -m MFP -bb 1000 -nt 10 -redo -quiet --prefix Single_tree/{name}_AA --seed 1', shell=True)
    subprocess.call(f'nw_topology -b -I Single_tree/{name}_AA.contree > deal_bootstrap_tree/{name}_deal.nwk', shell=True)
    
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


    p = Pool(5)
    tree_dir = 'Single_tree'
    alignment_dir = 'mafft_out'
    single_seqdir = 'Single_Copy_cds'

    for path in [tree_dir, alignment_dir, single_seqdir]:
        if not os.path.isdir(path):
            os.mkdir(path)

    # rows = open('SingleCopySGs.txt').readlines()
    single = [x.strip() for x in open('SingleCopySGs.txt')]
    rows = open('sample.list.SG.pan.modified').readlines()
    for i in range(0, len(rows)):
        if 'Orthogroup' in rows[i]:continue
        lis = re.split('[\\s,]+', rows[i].strip())
        if lis[0] in single:
            genes = [ge for x in lis[1:] for ge in x.split(', ')]
            extract(lis[0], genes, peps, single_seqdir)
            # extract(lis[0], lis[1:], peps, single_seqdir)
            p.apply_async(run_iqtree, args=(lis[0],))
    p.close()
    p.join()



