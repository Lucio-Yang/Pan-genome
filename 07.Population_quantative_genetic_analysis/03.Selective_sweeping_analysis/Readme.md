##Identification of domestication and improvement region for tetraploid wheat<br>
For selective analysis, we identified the selective signature based on SNP&InDel using Fst, π(Nucleotide diversity) and XPCLR value. And also, we calculated the SV frequency between three subpopulations (WEW, Wild emmer wheat; DTW, domesticated tetraploid wheat; FTW, free-threshing tetraploid wheat) to identift the selective SVs.

**<em>F</em>st analysis**
```
########### Fst for pairs of each subspecies
#!/usr/bin/bash
vcf=TW.sample742.filtered.snp-indel.nochr.miss0.8.maf0.05.bg.header.sorted.vcf.gz
for i in `seq 1 45`
do
pop1=$(sed -n "${i}p" pop_vs_pop.txt | awk '{print$1}')
pop2=$(sed -n "${i}p" pop_vs_pop.txt | awk '{print$2}')
vcftools --gzvcf ${vcf} --weir-fst-pop pop_${pop1}.list --weir-fst-pop pop_${pop2}.list --out fst_${pop1}_${pop2} --fst-window-size 500000 --fst-window-step 50000
done

########### Fst for three main subgroup
########### 01. Getting Fst value of each subgroup
vcf=TW.sample742.filtered.snp-indel.nochr.miss0.8.maf0.05.bg.header.sorted.vcf.gz

pop1=WEW
pop2=DTW
pop3=FTW

vcftools --gzvcf ${vcf} --weir-fst-pop pop_${pop1}.list --weir-fst-pop pop_${pop2}.list --out fst_${pop1}_${pop2} --fst-window-size 500000 --fst-window-step 50000
vcftools --gzvcf ${vcf} --weir-fst-pop pop_${pop1}.list --weir-fst-pop pop_${pop3}.list --out fst_${pop1}_${pop3} --fst-window-size 500000 --fst-window-step 50000
vcftools --gzvcf ${vcf} --weir-fst-pop pop_${pop2}.list --weir-fst-pop pop_${pop3}.list --out fst_${pop2}_${pop3} --fst-window-size 500000 --fst-window-step 50000

########### 02. Extracting each chrom
for i in `seq 1 14`
do
awk '$1 == '''$i''' {print $0}' fst_WEW_DTW.windowed.weir.fst > fst_WEW_DTW.fst.chr${i}.txt
awk '$1 == '''$i''' {print $0}' fst_DTW_FTW.windowed.weir.fst > fst_DTW_FTW.fst.chr${i}.txt
done

```
**Nucleotide diversity analysis**
```
########### Nucleotide diversity for each subspecies
vcf=TW.sample742.filtered.snp-indel.nochr.miss0.8.maf0.05.bg.header.sorted.vcf.gz
cat pop.list| while read f
do
bcftools view -S pop_${f}.list ${vcf} --threads 50 -Oz > variant_${f}.vcf.gz
plink --file ${bfile} --noweb --keep pop_${f}.list --recode vcf-iid bgz --keep-allele-order --out variant_${f}
vcftools --gzvcf variant_${f}.vcf.gz --out pi_${f} --window-pi 500000 --window-pi-step 50000
for i in 1 3 5 7 9 11 13; do awk '$1 == '''${i}'''{print$0}' pi_${f}.windowed.pi >> pi_${f}.windowed.Agenome.pi; done
for i in 2 4 6 8 10 12 14; do awk '$1 == '''${i}'''{print$0}' pi_${f}.windowed.pi >> pi_${f}.windowed.Bgenome.pi; done
done

########### Nucleotide diversity for three main subgroups
########### 01. Getting pi value
vcf=TW.sample742.filtered.snp-indel.nochr.miss0.8.maf0.05.bg.header.sorted.vcf.gz

for pop in WEW DEW DTW
do
bcftools view -S pop_${pop}.list ${vcf} --threads 50 -Oz > variant_${pop}.vcf.gz
tabix -p vcf -C variant_${pop}.vcf.gz
done

vcftools --gzvcf variant_WEW.vcf.gz --out pi_WEW --window-pi 500000 --window-pi-step 50000
vcftools --gzvcf variant_DTW.vcf.gz --out pi_DTW --window-pi 500000 --window-pi-step 50000
vcftools --gzvcf variant_FTW.vcf.gz --out pi_FTW --window-pi 500000 --window-pi-step 50000

############ 02. Calculating pi ratio
python3 cal_piratio.py pi_WEW.windowed.pi pi_DTW.windowed.pi > piratio_WEW_DTW.txt
python3 cal_piratio.py pi_DTW.windowed.pi pi_FTW.windowed.pi > piratio_DTW_FTW.txt

############ 03. Extracting each chrom
for i in `seq 1 14`
do
awk '$1 == '''$i''' {print$0}' piratio_WEW_DTW.txt > piratio_WEW_DTW.chr${i}
awk '$1 == '''$i''' {print$0}' piratio_DTW_FTW.txt > piratio_DTW_FTW.chr${i}
done

```
**XPCLR analysis**
```
pop1=pop_WEW.list
pop2=pop_DTW.list
pop3=pop_FTW.list

for f in `seq 1 14`
do
	vcf=chr${f}.vcf.gz
	xpclr --format vcf --input ${vcf} --out ./chr${f}_WEW_DTW_sweep.txt --samplesA ${pop1} --samplesB ${pop2} --size 500000 --step 50000 --maxsnps 500 --chr ${f} -V 10
	xpclr --format vcf --input ${vcf} --out ./chr${f}_DTW_FTW_sweep.txt --samplesA ${pop2} --samplesB ${pop3} --size 500000 --step 50000 --maxsnps 500 --chr ${f} -V 10
done
```

**Running GenWin for potential selective region**
```
for i in `seq 1 14`; do Rscript genwin.r chr${i}_WEW_DTW_sweep.txt ${i}; done
for i in `seq 1 14`; do Rscript genwin.r chr${i}_DTW_FTW_sweep.txt ${i}; done
for i in `seq 1 14`; do Rscript genwin.r fst_WEW_DTW.fst.chr${i}.txt ${i}; done
for i in `seq 1 14`; do Rscript genwin.r fst_DTW_FTW.fst.chr${i}.txt ${i}; done
for i in `seq 1 14`; do Rscript genwin.r piratio_WEW_DTW.chr${i} ${i}; done
for i in `seq 1 14`; do Rscript genwin.r piratio_DTW_FTW.chr${i} ${i}; done

head -1 chr1_WEW_DTW_sweep.txt_genwin.txt > XPCLR_domes_genwin.txt
for i in `seq 1 14`; do sed -n '2,$p' chr${i}_WEW_DTW_sweep.txt_genwin.txt >> XPCLR_domes_genwin.txt; done
head -1 XPCLR_domes_genwin.txt > XPCLR_domes_genwin.top0.05.txt
x=$(wc -l < XPCLR_domes_genwin.txt)
y=$(echo "$x - 1" | awk '{printf "%.0f", $1 * 0.05}')
sort -gk5 -r XPCLR_domes_genwin.txt | head -$y >> XPCLR_domes_genwin.top0.05.txt

head -1 chr1_DTW_FTW_sweep.txt_genwin.txt > XPCLR_improve_genwin.txt
for i in `seq 1 14`; do sed -n '2,$p' chr${i}_dew_dw_sweep.txt_genwin.txt >> XPCLR_improve_genwin.txt; done
head -1 XPCLR_improve_genwin.txt > XPCLR_improve_genwin.top0.05.txt
x=$(wc -l < XPCLR_improve_genwin.txt)
y=$(echo "$x - 1" | awk '{printf "%.0f", $1 * 0.05}')
sort -gk5 -r XPCLR_improve_genwin.txt | head -$y >> XPCLR_improve_genwin.top0.05.txt

head -1 fst_WEW_DTW.fst.chr1.txt_genwin.txt > fst_domes_genwin.txt
for i in `seq 1 14`; do sed -n '2,$p' fst_WEW_DTW.fst.chr${i}.txt_genwin.txt >> fst_domes_genwin.txt; done
head -1 fst_domes_genwin.txt > fst_domes_genwin.top0.05.txt
x=$(wc -l < fst_domes_genwin.txt)
y=$(echo "$x - 1" | awk '{printf "%.0f", $1 * 0.05}')
sort -gk5 -r fst_domes_genwin.txt | head -$y >> fst_domes_genwin.top0.05.txt

head -1 fst_DTW_FTW.fst.chr1.txt_genwin.txt > fst_improve_genwin.txt
for i in `seq 1 14`; do sed -n '2,$p' fst_DTW_FTW.fst.chr${i}.txt_genwin.txt >> fst_improve_genwin.txt; done
head -1 fst_improve_genwin.txt > fst_improve_genwin.top0.05.txt
x=$(wc -l < fst_improve_genwin.txt)
y=$(echo "$x - 1" | awk '{printf "%.0f", $1 * 0.05}')
sort -gk5 -r fst_improve_genwin.txt | head -$y >> fst_improve_genwin.top0.05.txt

head -1 piratio_WEW_DTW.chr1_genwin.txt > piratio_domes_genwin.txt
for i in `seq 1 14`; do sed -n '2,$p' piratio_WEW_DTW.chr${i}_genwin.txt >> piratio_domes_genwin.txt; done
head -1 piratio_domes_genwin.txt > piratio_domes_genwin.top0.05.txt
x=$(wc -l < piratio_domes_genwin.txt)
y=$(echo "$x - 1" | awk '{printf "%.0f", $1 * 0.05}')
sort -gk5 -r piratio_domes_genwin.txt | head -$y >> piratio_domes_genwin.top0.05.txt

head -1 piratio_DTW_FTW.chr1_genwin.txt > piratio_improve_genwin.txt
for i in `seq 1 14`; do sed -n '2,$p' piratio_DTW_FTW.chr${i}_genwin.txt >> piratio_improve_genwin.txt; done
head -1 piratio_improve_genwin.txt > piratio_improve_genwin.top0.05.txt
x=$(wc -l < piratio_improve_genwin.txt)
y=$(echo "$x - 1" | awk '{printf "%.0f", $1 * 0.05}')
sort -gk5 -r piratio_improve_genwin.txt | head -$y >> piratio_improve_genwin.top0.05.txt

awk -v OFS='\t' '{print$1,$2,$3}' fst_domes_genwin.top0.05.bed > domes_genwin.top0.05.bed
awk -v OFS='\t' '{print$1,$2,$3}' piratio_domes_genwin.top0.05.bed >> domes_genwin.top0.05.bed
awk -v OFS='\t' '{print$1,$2,$3}' XPCLR_domes_genwin.top0.05.bed >> domes_genwin.top0.05.bed
sort -V domes_genwin.top0.05.bed > domes_genwin.top0.05.sorted.bed
bedtools merge -i domes_genwin.top0.05.sorted.bed > domes_genwin.top0.05.sorted.merged.bed

awk -v OFS='\t' '{print$1,$2,$3}' fst_improve_genwin.top0.05.bed > improve_genwin.top0.05.bed
awk -v OFS='\t' '{print$1,$2,$3}' piratio_improve_genwin.top0.05.bed >> improve_genwin.top0.05.bed
awk -v OFS='\t' '{print$1,$2,$3}' XPCLR_improve_genwin.top0.05.bed >> improve_genwin.top0.05.bed
sort -V improve_genwin.top0.05.bed > improve_genwin.top0.05.sorted.bed
bedtools merge -i improve_genwin.top0.05.sorted.bed > improve_genwin.top0.05.sorted.merged.bed

```

**Identification of selective SVs based on Fisher's exact test**
```
bcftools view -S pop_WEW.list TWGG_sample742_SV.vcf.gz --threads 50 -Oz > TWGG_sample742_SV.WEW.vcf.gz
bcftools view -S pop_FTW.list TWGG_sample742_SV.vcf.gz --threads 50 -Oz > TWGG_sample742_SV.FTW.vcf.gz
bcftools view -S pop_DTW.list TWGG_sample742_SV.vcf.gz --threads 50 -Oz > TWGG_sample742_SV.DTW.vcf.gz

plink --vcf TWGG_sample742_SV.WEW.vcf.gz --keep-allele-order --freq --out TWGG_sample742_SV.WEW
plink --vcf TWGG_sample742_SV.FTW.vcf.gz --keep-allele-order --freq --out TWGG_sample742_SV.FTW
plink --vcf TWGG_sample742_SV.DTW.vcf.gz --keep-allele-order --freq --out TWGG_sample742_SV.DTW

python3 cal_fisher.py TWGG_sample742_SV.WEW.frq TWGG_sample742_SV.DTW.frq WEW-DTW > sv_WEW-DTW.freq
python3 cal_fisher.py TWGG_sample742_SV.DTW.frq TWGG_sample742_SV.FTW.frq DTW-FTW > sv_DTW-FTW.freq

python3 merge_domes-improve.py sv_WEW-DTW.freq sv_DTW-FTW.freq > sv_WEW-DTW-FTW.freq
python3 filter_selective_sv.py sv_WEW-DTW-FTW.freq > sv_WEW-DTW-FTW.freq.txt
python3 extract_fav-sv.py sv_WEW-DTW-FTW.freq.txt > sv_fav.freq.txt
```