## Genome preparing
**Spliting the chromsome of reference genome fasta and annotation file**

```
#!/usr/bin/bash

#####  split_fasta.bed1 #####
1A	1
1B	2
2A	3
2B	4
3A	5
3B	6
4A	7
4B	8
5A	9
5B	10
6A	11
6B	12
7A	13
7B	14
#####  split_fasta.bed2  #####
1	chr1	chr1_part1	1A	1	300801990
10	chr5	chr5_part2	3A	380730203	761460405

11	chr6	chr6_part1	3B	1	433323414
12	chr6	chr6_part2	3B	433323415	866646828
13	chr7	chr7_part1	4A	1	386494760
14	chr7	chr7_part2	4A	386494761	772976832
15	chr8	chr8_part1	4B	1	351312746
16	chr8	chr8_part2	4B	351312747	702625493
17	chr9	chr9_part1	5A	1	362006400
18	chr9	chr9_part2	5A	362006401	724012800
19	chr10	chr10_part1	5B	1	369959980
2	chr1	chr1_part2	1A	300801991	601435485
20	chr10	chr10_part2	5B	369959981	739917088
21	chr11	chr11_part1	6A	1	314263251
22	chr11	chr11_part2	6A	314263252	628526502
23	chr12	chr12_part1	6B	1	363350500
24	chr12	chr12_part2	6B	363350501	726696819
25	chr13	chr13_part1	7A	1	376919097
26	chr13	chr13_part2	7A	376919098	753838194
27	chr14	chr14_part1	7B	1	384208404
28	chr14	chr14_part2	7B	384208405	768416808
3	chr2	chr2_part1	1B	1	359947670
4	chr2	chr2_part2	1B	359947671	719895341
5	chr3	chr3_part1	2A	1	399021356
6	chr3	chr3_part2	2A	399021357	798042713
7	chr4	chr4_part1	2B	1	414443095
8	chr4	chr4_part2	2B	414443096	828886191
9	chr5	chr5_part1	3A	1	380730202
#############################

cat split_fasta.bed1 | while read f
do
	c1=$(echo ${f} | awk '{print$1}')
	c2=$(echo ${f} | awk '{print$2}')
	echo '>'${c2} > Kronos.chr${c2}.fasta
	samtools faidx Kronos.fasta ${c1} | sed -n '2,$p' >> Kronos.chr${c2}.fasta
done

cat split_fasta.bed2 | while read f
do
	c1=$(echo ${f} | awk '{print$1}')
	c2=$(echo ${f} | awk '{print$4}')
	l1=$(echo ${f} | awk '{print$5}')
	l2=$(echo ${f} | awk '{print$6}')
	echo '>'${c1} >> Kronos_split.fasta
	samtools faidx Kronos.fasta ${c2}:${l1}-${l2} | sed -n '2,$p' >> Kronos_split.fasta
done
samtools faidx Kronos_split.fasta
awk '{print$1}' Kronos_split.fasta.fai | while read f
do
	samtools faidx Kronos_split.fasta ${f} > TW_t2_split.chr${f}.fasta
	samtools faidx Kronos_split.chr${f}.fasta
done
```