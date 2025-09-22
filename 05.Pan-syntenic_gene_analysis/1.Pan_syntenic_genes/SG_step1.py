#!/usr/bin/python3
import sys,random

sample = sys.argv[1]
sg = sys.argv[2]
num = int(sys.argv[3])

t = {}

count = 1
for l in open(sample,'r'):
	m = l.rstrip()
	t[m] = count
	count += 1

result = random.sample(range(1,13),num)

for line in open(sg,'r'):
	mes = line.rstrip().split('\t')
	mes2 = mes[4:]
	t2 = []
	for i in mes2:
		s = i.split('_')[0]
		if s:
			s_id = t[s]
			if s_id in result:
				t2.append(i)
	if t2:
		print('\t'.join(t2))

		
