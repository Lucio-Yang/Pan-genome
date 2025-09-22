#!/usr/bin/python3
import sys

tmp = sys.argv[1]
num = int(sys.argv[2])

core=0
pan=0

for line in open(tmp,'r'):

	mes = line.rstrip().split('\t')
	n = len(mes)
	if n == num:
		core += 1
	else:
		pan += 1

print('core = ',core,sep = '')
print('pan = ',pan,sep = '')
