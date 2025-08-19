# -*- coding: utf-8 -*-

import re
import os
import glob
import subprocess
tree_dir = 'Single_tree'  # Single_tree
change_name_dir = 'newname-Single_tree'
reroot_dir = 'reroot_tree'
deal_length_dir = 'deal_bootstrap_tree'

for path in [change_name_dir, reroot_dir, deal_length_dir]:
    if not os.path.isdir(path):
        os.mkdir(path)
single = [x.strip() for x in open('SingleCopySGs.txt')]
rows = open('sample.list.SG.pan.modified').readlines()

names = []
for i in range(1, len(rows)):
    lis = re.split('[\\s,]+', rows[i].strip())
    if not os.path.isfile(f'{tree_dir}/{lis[0]}_AA.nwk'): continue
    if lis[0] not in single: continue
    tr = open(f'{tree_dir}/{lis[0]}_AA.nwk').read()
    sf = open(f'{change_name_dir}/{lis[0]}.tre', 'w')
    # sf = open(f'ch_name_tree/{lis[0]}.tre', 'w')
    names.append(lis[0])
    for ge in lis[1:]:
        if 'Hv' in ge[:2]:
            name = 'Hv'
        else:
            name = ge.split('_')[0]
        # name = ge[:4]
        tr = tr.replace(ge, name)
    sf.write(tr)
# rows = open('contain14-species.single-copy_groups.tsv').readlines()
    
for na in names:
    os.system(f'nw_reroot {change_name_dir}/{na}.tre Hv > {reroot_dir}/{na}.tree')
    os.system(f'nw_topology {reroot_dir}/{na}.tree > {deal_length_dir}/{na}_deal.tree')

subprocess.call(f'cat {reroot_dir}/*.tree > all.reroot.tre', shell=True)


subprocess.call(f'cat {change_name_dir}/*.tre > {change_name_dir}_gene_trees.tree', shell=True)
