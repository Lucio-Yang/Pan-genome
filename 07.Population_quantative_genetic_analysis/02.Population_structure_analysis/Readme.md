## Population structure analyis
** Preparing for input files**
```
#!/bin/bash
#SBATCH -o job.run.%j.out
#SBATCH -p m256Partition
#SBATCH -J run
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

#### Merging SNP and InDel vcf
tabix -p vcf -C TW.sample742.filtered.indel.nochr.miss0.8.maf0.05.bg.header.vcf.gz
tabix -p vcf -C TW.sample742.filtered.snp.nochr.miss0.8.maf0.05.bg.header.vcf.gz
tabix -p vcf -C TWGG_sample742_SV.vcf.gz
bcftools concat TW.sample742.filtered.snp.nochr.miss0.8.maf0.05.bg.header.vcf.gz TW.sample742.filtered.indel.nochr.miss0.8.maf0.05.bg.header.vcf.gz -a --threads 100 -Oz > TW.sample742.filtered.snp-indel.nochr.miss0.8.maf0.05.bg.header.vcf.gz
bcftools concat TW.sample742.filtered.snp-indel.nochr.miss0.8.maf0.05.bg.header.sorted.vcf.gz TWGG_sample742_SV.vcf.gz -a --threads 100 -Oz > TW.sample742.filtered.snp-indel-sv.nochr.miss0.8.maf0.05.bg.header.changegeno.vcf.gz
bcftools sort TW.sample742.filtered.snp-indel.nochr.miss0.8.maf0.05.bg.header.vcf.gz -m 250G -T ./tmp -Oz > TW.sample742.filtered.snp-indel.nochr.miss0.8.maf0.05.bg.header.sorted.vcf.gz
bcftools sort TW.sample742.filtered.snp-indel-sv.nochr.miss0.8.maf0.05.bg.header.vcf.gz -m 250G -T ./tmp -Oz > TW.sample742.filtered.snp-indel-sv.nochr.miss0.8.maf0.05.bg.header.sorted.vcf.gz
tabix -p vcf -C TW.sample742.filtered.snp-indel.nochr.miss0.8.maf0.05.bg.header.sorted.vcf.gz
tabix -p vcf -C TW.sample742.filtered.snp-indel-sv.nochr.miss0.8.maf0.05.bg.header.sorted.vcf.gz

#### Converting plink format
sample=TW.sample742.filtered.snp-indel-sv.nochr.miss0.8.maf0.05.bg.header.sorted
vcftools --gzvcf ${sample}.vcf.gz --plink --out ${sample}
plink --file $sample --make-bed --out $sample
plink --bfile $sample --indep-pairwise 50 10 0.2 --out ${sample}.indep-pairwisefilter
plink --bfile $sample --extract ${sample}.indep-pairwisefilter.prune.in --make-bed --out ${sample}.indep-pairwisefilter
plink --bfile ${sample}.indep-pairwisefilter --recode --out ${sample}.indep-pairwisefilter

#### Extracting A genome vcf and B genome vcf for LD calculation
plink --bfile TW.sample742.filtered.snp.nochr.miss0.8.maf0.05.bg.header --chr 2 4 6 8 10 12 14 --make-bed --out TW.sample742.noout.filtered.snp.nochr.miss0.8.maf0.05.bg.header.Bgenomes
plink --bfile TW.sample742.filtered.snp.nochr.miss0.8.maf0.05.bg.header --chr 1 3 5 7 9 11 13 --make-bed --out TW.sample742.noout.filtered.snp.nochr.miss0.8.maf0.05.bg.header.Agenomes
plink --bfile TW.sample742.filtered.indel.nochr.miss0.8.maf0.05.bg.header --chr 2 4 6 8 10 12 14 --make-bed --out TW.sample742.noout.filtered.indel.nochr.miss0.8.maf0.05.bg.header.Bgenomes
plink --bfile TW.sample742.filtered.indel.nochr.miss0.8.maf0.05.bg.header --chr 1 3 5 7 9 11 13 --make-bed --out TW.sample742.noout.filtered.indel.nochr.miss0.8.maf0.05.bg.header.Agenomes
plink --bfile TW.sample742.filtered.snp-indel.nochr.miss0.8.maf0.05.bg.header.sorted --chr 2 4 6 8 10 12 14 --make-bed --out TW.sample742.filtered.snp-indel.nochr.miss0.8.maf0.05.bg.header.sorted.changeindel.Bgenomes
plink --bfile TW.sample742.filtered.snp-indel.nochr.miss0.8.maf0.05.bg.header.sorted --chr 1 3 5 7 9 11 13 --make-bed --out TW.sample742.filtered.snp-indel.nochr.miss0.8.maf0.05.bg.header.sorted.changeindel.Agenomes
plink --bfile TW.sample742.filtered.snp-indel-sv.nochr.miss0.8.maf0.05.bg.header.sorted --chr 2 4 6 8 10 12 14 --keep-allele-order --make-bed --out TW.sample742.filtered.snp-indel-sv.nochr.miss0.8.maf0.05.bg.header.sorted.Bgenomes
plink --bfile TW.sample742.filtered.snp-indel-sv.nochr.miss0.8.maf0.05.bg.header.sorted --chr 1 3 5 7 9 11 13 --keep-allele-order --make-bed --out TW.sample742.filtered.snp-indel-sv.nochr.miss0.8.maf0.05.bg.header.sorted.Agenomes
```
**Construction of phylogenetic trees**
```
#### Extracting for random variants
python snp_number.py TW.sample742.filtered.snp.nochr.miss0.8.maf0.05.bg.header.indep-pairwisefilter.map

#### Running raxml
raxmlHPC-PTHREADS-AVX2 -s random_10w_snp.fasta -n random_10w_snp -f a -m GTRGAMMA -T 128 -p 123 -x 123 -# 100
```

**Running ADMIXTURE**
```
for K in `seq 2 20`
do	
	admixture --cv TW.sample742.filtered.snp.nochr.miss0.8.maf0.05.bg.header.indep-pairwisefilter.bed $K -j52 | tee log${K}.out
done
```

**Running PCA**
```
gcta --bfile TW.sample742.filtered.snp.nochr.miss0.8.maf0.05.bg.header.indep-pairwisefilter --autosome --make-grm --out TW.sample742.filtered.snp.nochr.miss0.8.maf0.05.bg.header.indep-pairwisefilter --thread-num 128
```

**Calculating LD decay**
```
for p in WEW DTW FTW
do
	for i in Agenomes Bgenomes
    do
		PopLDdecay -InVCF snp_sample742_${i}_${p}.vcf.gz -MaxDist 3000 -OutStat SNP_sample742_${i}_${p}
    	PopLDdecay -InVCF indel_sample742_${i}_${p}.vcf.gz -MaxDist 3000 -OutStat SNP_sample742_${i}_${p}
    	PopLDdecay -InVCF sv_sample742_${i}_${p}.vcf.gz -MaxDist 3000 -OutStat SNP_sample742_${i}_${p}
    	PopLDdecay -InVCF snp-indel-sv_sample742_${i}_${p}.vcf.gz -MaxDist 3000 -OutStat snp-indel-sv_sample742_${i}_${p}
	done
done
```
