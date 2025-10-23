## Heritability estimating
We estimated heritability for 32 traits based on SNP, InDel and SV respectively according to [previous study](https://www.nature.com/articles/s41586-022-04808-9). The following example is merely exemplified by SNP variants.
```
ldak --bfile snp --window-prune 0.98 --window-kb 100 --cut-weights snp_weights
awk < snp_weights/thin.in '{print $1, 1}' > snp_weights.thin

for f in `seq 1 32`
do
mkdir -p pheno${f}
	echo "#!/bin/bash
	#SBATCH -o job.${f}_snp_ldak.%j.out
	#SBATCH -p tcum256c128Partition
	#SBATCH -J ${f}_snp_ldak
	#SBATCH --nodes=1
	#SBATCH --ntasks-per-node=50" > pheno${f}/run_${f}_snp_ldak.sh
	echo "ldak --linear ${f}_snp_quant --bfile ../snp --pheno ../../gcta_pheno/total.pheno --mpheno ${f} --max-threads 50" >> pheno${f}/run_${f}_snp_ldak.sh
    
	echo "for j in {-20..10}; do
	alpha=\`echo \$j | awk '{print -\$1/20}'\`; echo \$alpha
	ldak --calc-tagging ${f}_snp.\$j --bfile ../snp --power \$alpha --window-kb 100 --max-threads 50
	done" >> pheno${f}/run_${f}_snp_ldak.sh
    
	echo "for j in {-20..10}; do
	ldak --sum-hers ${f}_snp_alpha.\$j --summary ${f}_snp_quant.summaries --tagfile ${f}_snp.\$j.tagging --check-sums NO --max-threads 50
	done" >> pheno${f}/run_${f}_snp_ldak.sh
    
	echo "for j in {-20..10}; do
	alpha=\`echo \$j | awk '{print -\$1/20}'\`;
	grep Alt_logl ${f}_snp_alpha.\$j.extra | awk -v alpha=\$alpha '{print alpha, \$2}'
	done > ${f}_snp_alpha.likes" >> pheno${f}/run_${f}_snp_ldak.sh
    
	echo "ldak --find-gaussian ${f}_snp_alpha --likelihoods ${f}_snp_alpha.likes --max-threads 50" >> pheno${f}/run_${f}_snp_ldak.sh
    
	echo "ldak --calc-kins-direct LDAK-Thin_${f}_snp --bfile ../snp --weights ../snp_weights.thin --power -.5 --max-threads 50" >> pheno${f}/run_${f}_snp_ldak.sh
    
	echo "ldak --pheno ../../gcta_pheno/total.pheno --mpheno ${f}  --grm LDAK-Thin_${f}_snp --covar ../TW.gwas.snp.eigenvec.pca10 --reml ${f}_snp --constrain YES --max-threads 50" >> pheno${f}/run_${f}_snp_ldak.sh
done

```

## Genome wide association analysis
In this study, we performed GWAS for five subgroups, including all accessions, wild emmer wheat accessions, domesticated emmer wheat accessions, durum accessions, cultivated accessions (domesticated emmer wheat + durum). And also, we performed GWAS for SNPs, InDels and SVs respectively. The following example is merely exemplified by SNP variants.

**Preparing for covar and kinship files**
```
##### Extracting GWAS samples and filtering MAF
plink --bfile TW.sample742.filtered.snp.nochr.miss0.8.maf0.05.bg.header --noweb --keep gwas_sample.list --maf 0.05 --recode --make-bed --out TW.gwas.snp
plink --bfile TW.sample742.filtered.snp.nochr.miss0.8.maf0.05.bg.header --noweb --keep gwas_sample_WEW.list --maf 0.05 --recode --make-bed --out TW.gwas_wew.snp
plink --bfile TW.sample742.filtered.snp.nochr.miss0.8.maf0.05.bg.header --noweb --keep gwas_sample_DEW.list --maf 0.05 --recode --make-bed --out TW.gwas_dew.snp
plink --bfile TW.sample742.filtered.snp.nochr.miss0.8.maf0.05.bg.header --noweb --keep gwas_sample_DW.list --maf 0.05 --recode --make-bed --out TW.gwas_dw.snp
plink --bfile TW.sample742.filtered.snp.nochr.miss0.8.maf0.05.bg.header --noweb --keep gwas_sample_DEW-DW.list --maf 0.05 --recode --make-bed --out TW.gwas_dew-dw.snp
plink --bfile TW.gwas.snp --noweb --keep gwas_sample260.list --maf 0.05 --recode --make-bed --out TW.gwas_sample260.snp

##### Calculating the kinship matrix
gemma -bfile TW.gwas.snp -gk 2 -o kin.snp
gemma -bfile TW.gwas_WEW.snp -gk 2 -o kin.snp.WEW
gemma -bfile TW.gwas_DEW.snp -gk 2 -o kin.snp.DEW
gemma -bfile TW.gwas_DW.snp -gk 2 -o kin.snp.DW
gemma -bfile TW.gwas_DEW-DW.snp -gk 2 -o kin.snp.DEW-DW

##### Calculating the covar
for bfile in TW.gwas.snp TW.gwas_WEW.snp TW.gwas_DEW.snp TW.gwas_DW.snp TW.gwas_DEW-DW.snp
do
gcta --bfile ${bfile} --autosome --make-grm --out ${bfile} --thread-num 128
gcta  --grm ${bfile} --pca 20  --out ${bfile} --thread-num 128
awk '{print"1",$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' ${bfile}.eigenvec > ${bfile}.eigenvec.pca10.txt
done
```
**Performing GWAS**
```
path=~/01.SNP
bfile=TW.gwas.snp
kin=~/kin.snp.sXX.txt
pca=~/TW.gwas.snp.eigenvec.pca10.txt

cat traits.id | while read f
do
mkdir -p ${f}
cd ${f}
        for i in `seq 1 14`
        do
        	mkdir -p chr${i}
        	cd chr${i}
        	ln -s ${path}/${bfile}.chr${i}.bed
        	ln -s ${path}/${bfile}.chr${i}.bim
        	python3 make_pheno_fam.py ${f}.gwas_input.txt ${path}/${bfile}.chr${i}.fam > ${bfile}.chr${i}.fam
        	/usr/bin/time -v gemma -bfile ${bfile}.chr${i} -k ${kin} -lmm 4 -c ${pca} -miss 1 -o ${f}
        cd ../
        done
cd ../
done
```