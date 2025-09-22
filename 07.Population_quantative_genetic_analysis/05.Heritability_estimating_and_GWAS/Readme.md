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
