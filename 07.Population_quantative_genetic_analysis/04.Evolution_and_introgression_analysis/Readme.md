### Population analysis

**Network analysis**

Network-based approach to explore the species history and community structure. 
```text
##### random SNP sites selection
python evolution_network_random_SNP.py

cat id | while read f
do
	shuf -n 2000 ${f} > ${f}.shuf
	awk '{print$2}' ${f}.shuf > ${f}.shuf.list
	plink --bfile TW_HW_merged.SNP.Hv.nochr.imput.maf0.05.indep-pairwisefilter --extract ${f}.shuf.list --recode --out ${f}
	vcftools --vcf ${f}.vcf --plink --out ${f}
	perl ped2fa.pl ${f}.ped ${f}.fasta
done

#construction the maximum likelihood trees
raxmlHPC-PTHREADS-SSE3 -s random1.fas -n random1.out -m GTRGAMMA --JC69 --asc-corr=felsenstein -T 24 -f a -N 1000 -x 12345 -p 12345

#MST construction
python MST-Kruskal.py
```

#### ABBA-BABA analysis

Inferring the gene-flows among the tetraploid wheats.

```text
##### pop_list.sample.txt #####
A1890_A1890	WEW
A1897_A1897	WEW
A1901_A1901	WEW
A1912_A1912	WEW
A3108_A3108	DW
A3109_A3109	DW
A3113_A3113	DW
A3316_A3316	DEW
A3317_A3317	DEW
A3319_A3319	DEW
A3320_A3320	DEW
SRR12640144.man_SRR12640144.man	Rivet
SRR12640145.man_SRR12640145.man	Rivet
SRR12640146.man_SRR12640146.man	Rivet
SRR12640147.man_SRR12640147.man	Rivet
#############################

example commond:
python genomics_general-master/ABBABABAwindows.py -g my_data/chr1_output.geno.gz -f phased -o results_Ispahanicum/chr1.wheat.ABBABABA_Ispahanicum_DW_DEW_out.w100m50.csv.gz -P1 pop1 -P2 pop2 -P3 pop3 -O out --popsFile my_data/pop_list.sample513.txt --windType sites  -w 100 -s 50 -m 3 --T 128
```

Calculate the proportion of genome introgression (PGI) among tetraploid wheats.

```text
bash calculate_PGI.sh
```

