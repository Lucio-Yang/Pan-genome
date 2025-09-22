#!/usr/bin/python3

import sys,gzip
import numpy as np
from scipy.stats import fisher_exact

freq1 = sys.argv[1]
freq2 = sys.argv[2]
group = sys.argv[3]

pop1 = {}
print("SV_id","POP1_ref","POP1_alt","POP1_freq","POP2_ref","POP2_alt","POP2_freq","P-value","Group",sep = '\t')

for line in open(freq1, 'r').readlines()[1:]:

	mes = line.rstrip().split()
	if not line.startswith('#'):
		f_alt = float(mes[-2])	### AC
		f_ref = 1-f_alt	### AN-AC
		sv_id = mes[1]
		num = int(mes[-1])
		n_alt = int(f_alt * num)
		n_ref = int(f_ref * num)
#		print(sv_id,ref,alt)
		pop1[sv_id] = {}
		pop1[sv_id]['ref'] = n_ref
		pop1[sv_id]['alt'] = n_alt
		pop1[sv_id]['n'] = num

for l in open(freq2, 'r').readlines()[1:]:
	
	m = l.rstrip().split()
	if not l.startswith('#'):
		sv_id2 = m[1]
		f_a = float(m[-2])        ### AC
		f_r = 1-f_a       ### AN-AC
		n = int(m[-1])
		n_a = int(f_a * n)
		n_r = int(f_r * n)
		data = np.array([[pop1[sv_id2]['ref'],n_r],[pop1[sv_id2]['alt'],n_a]])
		odds_ratio, p_value = fisher_exact(data)
		freq1 = pop1[sv_id2]['alt']/pop1[sv_id2]['n']
		freq2 = f_a
		print(sv_id2,pop1[sv_id2]['ref'],pop1[sv_id2]['alt'],freq1,n_r,n_a,freq2,p_value,group,sep = '\t')


