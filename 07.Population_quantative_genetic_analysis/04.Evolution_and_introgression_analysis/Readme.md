### Population analysis

#### Network analysis

The network-based approach used for [bread wheat](https://www.nature.com/articles/s41588-019-0393-z) was used to reconstruct the species history and community structure of tetraploid wheat.Based on a high-quality VCF file, we divided the genome into 1,000 parts and performed repeated random-haplotype sampling at heterozygous sites, randomly sampling 2,000 SNP sites each time to infer 1,000 maximum-likelihood tree topologies with the GTRGAMMA model and JC69 distances in RAxML (asc-corr = felsenstein). While these repeated random haplotype sampling procedures including heterozygous sites (RRHS) trees were also analysed in the form of conventional consensus topologies and Densitree visualizations to infer taxonomic clades, we used the Kruskal algorithm130 implemented in Python to analyse the evolutionary distances between each pair of accessions using 1,000 minimum spanning trees (MST).

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

The D statistic was used to reflect introgression proportions at the genome-wide level, and the fd statistic was used to estimate introgression proportions within a specific size window. To explore the introgression regions among the tetraploid wheat accessions, we calculated D and fd values across the whole genome using [Python code](https://github.com/simonhmartin/genomics_general), with the window size = 100 and step size = 50 SNPs. The minimum number of loci in each window was set to 3. The fd statistic ranges from 0 to 1, and the D statistic ranges from –1 to 1. When windows with D < 0, or D > 0, but fd > 1, the fd statistic was biologically meaningless or noisy; therefore, fd was set to zero.

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

