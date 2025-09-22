#!/usr/bin/python3

import sys


for line in open(sys.argv[1],'r'):

    if line.startswith('S'):
        print(line.rstrip(),'Type',sep = '\t')
    else:
        mes = line.split('\t')
        diff1 = abs(float(mes[6]) - float(mes[3]))
        diff2 = abs(float(mes[9]) - float(mes[6]))
        f_w = float(mes[3])
        f_d = float(mes[6])
        f_f = float(mes[9])
        p1 = float(mes[-2])
        p2 = float(mes[-1])
        t = []
        n1 = 0; n2 = 0
        t.append(f_w); t.append(f_d); t.append(f_f)
        for i in t:
            if i > 0.5:
                n1 += 1
            elif i < 0.5:
                n2 += 1
        if p1 <= 0.05 and p2 <= 0.05:
            if n1 < 3 and n2 < 3:
                print(line.rstrip(),'Domes-Improve',sep = '\t')
            else:
                print(line.rstrip(),'None-None',sep = '\t')
        elif p1 <= 0.05 and p2 > 0.05:
            if f_w > 0.5 and f_d < 0.5:
                print(line.rstrip(),'Domes-None',sep = '\t')
            elif f_w < 0.5 and f_d > 0.5:
                print(line.rstrip(),'Domes-None',sep = '\t')
            else:
                print(line.rstrip(),'None-None',sep = '\t')
        elif p1 > 0.05 and p2 <= 0.05:
            if f_d > 0.5 and f_f < 0.5:
                print(line.rstrip(),'None-Improve',sep = '\t')
            elif f_d < 0.5 and f_f > 0.5:
                print(line.rstrip(),'None-Improve',sep = '\t')
            else:
                print(line.rstrip(),'None-None',sep = '\t')
        elif p1 > 0.05 and p2 > 0.05:
            print(line.rstrip(),'None-None',sep = '\t')
#        else:
#            print(line.rstrip(),'None-None',sep = '\t')
