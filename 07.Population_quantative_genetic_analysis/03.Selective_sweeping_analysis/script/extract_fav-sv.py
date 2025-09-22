#!/usr/bin/python3

import sys


for line in open(sys.argv[1],'r'):
    if line.startswith('S'):
        print(line.rstrip(),'fav-SV',sep = '\t')
    else:
        mes = line.split('\t')
        f_w = float(mes[3])
        f_d = float(mes[6])
        f_f = float(mes[9])
        p1 = float(mes[-2])
        p2 = float(mes[-1])
        if 'Domes' in line or 'Improve' in line:
            if f_f >= f_d >= f_w:
                print(line.rstrip(),'Up',sep = '\t')
            elif f_f <= f_d <= f_w:
                print(line.rstrip(),'Down',sep = '\t')
            else:
                print(line.rstrip(),'None',sep = '\t')
