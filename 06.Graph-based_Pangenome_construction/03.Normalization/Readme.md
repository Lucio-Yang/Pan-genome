## Converting MAF to VCF
**Converting the output of AnchorWave to the gVCF format using TASSEL**
```
#!/usr/bin/bash

for chrom_id in `seq 1 14`
do
	cat cat specieslist | while read query
	do
	run_pipeline.pl -Xmx100G -debug -MAFToGVCFPlugin -referenceFasta /data/yangg/project/pangenome/03.pangenome/00.genome/02.each_chrom/Kronos/Kronos.chr${chrom_id}.fasta -mafFile Kronos_${query}_chr${chrom_id}.maf -sampleName Kronos_${query}_chr${chrom_id}_anchorwave -gvcfOutput ./Kronos_${query}_chr${chrom_id}.gvcf -fillGaps false -bgzipAndIndex false > ./Kronos_${query}_chr${chrom_id}_outputMafToGVCF.txt
	done
done
```
**Converting gVCF to VCF format**
```
cat sample.list | while read f; do python3 convert_gvcf2vcf.py ${f}.gvcf > ${f}.vcf; sed -i 's/1:0/1\/1:0/g' ${f}.vcf; python3 filter_length.py ${f}.vcf > ${f}_longerthan50bp.vcf; bgzip -@ 50 ${f}_longerthan50bp.vcf; tabix -C -p vcf ${f}_longerthan50bp.vcf.gz; done
```

## Variants merge
**Merging VCF files using BCFtools**
```
#!/usr/bin/bash

for i in `seq 1 14`
do
	bcftools merge ../02.gvcf2vcf/Kronos_Svevo_chr${i}_longerthan50bp.vcf.gz ../02.gvcf2vcf/Kronos_PI294478_chr${i}_longerthan50bp.vcf.gz ../02.gvcf2vcf/Kronos_NU00021_chr${i}_longerthan50bp.vcf.gz ../02.gvcf2vcf/Kronos_IG77365_chr${i}_longerthan50bp.vcf.gz ../02.gvcf2vcf/Kronos_IG99236_chr${i}_longerthan50bp.vcf.gz ../02.gvcf2vcf/Kronos_NU01905_chr${i}_longerthan50bp.vcf.gz ../02.gvcf2vcf/Kronos_NU01954_chr${i}_longerthan50bp.vcf.gz ../02.gvcf2vcf/Kronos_XM001097_chr${i}_longerthan50bp.vcf.gz ../02.gvcf2vcf/Kronos_Zavitan_chr${i}_longerthan50bp.vcf.gz | bcftools norm -m -any -N | bcftools norm -d none --fasta-ref ../../01.genome/Kronos.chr${i}.fasta | bcftools sort > chr${i}_longerthan50bp.merged.vcf
	python3 re-chrom_position.py chr${i}_longerthan50bp.merged.vcf > chr${i}_longerthan50bp.merged.split.vcf
	bcftools annotate --rename-chrs change_chrom.txt chr${i}_longerthan50bp.merged.split.vcf --force -Ov > chr${i}_longerthan50bp.merged.split.changechr.vcf
	bgzip -@ 50 chr${i}_longerthan50bp.merged.split.changechr.vcf
	tabix -C -p vcf chr${i}_longerthan50bp.merged.split.changechr.vcf.gz
done
```
**Spliting chromosomes**
```
##### Bed file #####
1	1
1	2
2	3
2	4
3	5
3	6
4	7
4	8
5	9
5	10
6	11
6	12
7	13
7	14
8	15
8	16
9	17
9	18
10	19
10	20
11	21
11	22
12	23
12	24
13	25
13	26
14	27
14	28
####################
#!/usr/bin/bash

cat bed | while read f
do
	c1=$(echo ${f} | awk '{print$1}')
	c2=$(echo ${f} | awk '{print$2}')
	bcftools view -r ${c2} chr${c1}_longerthan50bp.merged.split.changechr.vcf.gz > chr${c2}.vcf
done
```

## Merging and optimizing the multiallelic SVs using Panpop
```
#!/usr/bin/bash

for i in `seq 1 28`
do
	mkdir panpop_chr${i}_output
	perl PART_run.pl -i chr${i}.vcf -o panpop_chr${i}_output -r Kronos_split.chr${i}.fasta -t 50 --tmpdir ./tmp
done
```
**Sorting chromosomes ID according to ASCII**
```
#!/usr/bin/bash

zcat panpop_chr1_output/3.final.vcf.gz | head -16  > TWGG.vcf

for i in 10 11 12 13 14 15 16 17 18 19 2 20 21 22 23 24 25 26 27 28 3 4 5 6 7 8 9
do
zcat panpop_chr${i}_output/3.final.vcf.gz | grep "contig" >> TWGG.vcf
done

zcat panpop_chr1_output/3.final.vcf.gz | sed -n '17,28p' >> TWGG.vcf
for i in 1 10 11 12 13 14 15 16 17 18 19 2 20 21 22 23 24 25 26 27 28 3 4 5 6 7 8 9
do
zcat -c panpop_chr${i}_output/3.final.vcf.gz | sed -n '29,$p' >> TWGG.vcf
done
bgzip -@ 100 -c TWGG.vcf > TWGG.vcf.gz
tabix -p vcf TWGG.vcf.gz
```
