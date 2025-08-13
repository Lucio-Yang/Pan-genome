# TWGG: Tetraploid wheat graph-based pangenome



## Genome preparing
**Spliting the chromsome of reference genome fasta and annotation file**

```
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
	echo '>'${c2} > TW_t2.chr${c2}.fasta
	samtools faidx TW_t2.fasta ${c1} | sed -n '2,$p' >> TW_t2.chr${c2}.fasta
done

cat split_fasta.bed2 | while read f
do
	c1=$(echo ${f} | awk '{print$1}')
	c2=$(echo ${f} | awk '{print$4}')
	l1=$(echo ${f} | awk '{print$5}')
	l2=$(echo ${f} | awk '{print$6}')
	echo '>'${c1} >> TW_t2_split.fasta
	samtools faidx TW_t2.fasta ${c2}:${l1}-${l2} | sed -n '2,$p' >> TW_t2_split.fasta
done
samtools faidx TW_t2_split.fasta
awk '{print$1}' TW_t2_split.fasta.fai | while read f
do
	samtools faidx TW_t2_split.fasta ${f} > TW_t2_split.chr${f}.fasta
	samtools faidx TW_t2_split.chr${f}.fasta
done
```

## Alignments by [AnchorWave](https://github.com/baoxingsong/AnchorWave)
 ```
 #!/bin/bash
 chrom_number=14
 ref=TW_t2
 result_path=/home/yangg/my_data/project/pangenome/03.pangenome/01.anchorwave
 genome_path=/data/yangg/project/pangenome/03.pangenome/00.genome/02.each_chrom
 gff_path=/data/yangg/project/pangenome/03.pangenome/00.genome/03.gff
 specieslist=/data/yangg/project/pangenome/03.pangenome/01.anchorwave/genome.list


 for chrom_id in `seq 1 14`
 do
 	mkdir -p chr${chrom_id}
    cd chr${chrom_id}
    
    header1="#!/bin/bash\n#SBATCH -o job.run_anchorwave_ref.%j.out\n#SBATCH -p tcum256c128Partition\n#SBATCH -J TWGG\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=50\n"
    echo -e ${header1} > run_${ref}_chr${chrom_id}.sh
    echo -e "anchorwave gff2seq -r ${genome_path}/${ref}/${ref}.chr${chrom_id}.fasta -i ${gff_path}/${ref}/${ref}.chr${chrom_id}.gff -o ${ref}_chr${chrom_id}_cds.fa" >> run_${ref}_chr${chrom_id}.sh
    echo -e "minimap2 -x splice -t 15 -k 12 -a -p 0.4 -N 20 --cs ${genome_path}/${ref}/${ref}.chr${chrom_id}.fasta ${ref}_chr${chrom_id}_cds.fa > ${ref}_chr${chrom_id}.sam" >> run_${ref}_chr${chrom_id}.sh
	cat specieslist | while read query
	do
    	header2="#!/bin/bash\n#SBATCH -o job.${ref}_${query}_${chrom_id}.%j.out\n#SBATCH -p tcum256c128Partition\n#SBATCH -J TWGG\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=50\n"
    	echo -e ${header2} > run_${ref}_${query}_chr${chrom_id}.sh
    	echo "echo Start Time is \`date\`" >> run_${ref}_${query}_chr${chrom_id}.sh
    	echo "time1=\`date +\"%Y-%m-%d %H:%M:%S\"\`" >> run_${ref}_${query}_chr${chrom_id}.sh
    	echo "#---------------------------------------------------------------------------" >> run_${ref}_${query}_chr${chrom_id}.sh
    	echo -e "###ref is ${ref}" >> run_${ref}_${query}_chr${chrom_id}.sh
    	echo -e "minimap2 -x splice -t 50 -k 12 -a -p 0.4 -N 20 --cs ${genome_path}/${query}/${query}.chr${chrom_id}.fasta ${ref}_chr${chrom_id}_cds.fa > ${query}_chr${chrom_id}.sam\n" >> run_${ref}_${query}_chr${chrom_id}.sh
    	echo -e "anchorwave proali -r ${genome_path}/${ref}/${ref}.chr${chrom_id}.fasta -i ${gff_path}/${ref}/${ref}.chr${chrom_id}.gff -ar ${ref}_chr${chrom_id}.sam -as ${ref}_chr${chrom_id}_cds.fa -a ${query}_chr${chrom_id}.sam -s ${genome_path}/${query}/${query}.chr${chrom_id}.fasta -n ${ref}_${query}_chr${chrom_id}.anchor -o ${ref}_${query}_chr${chrom_id}.maf -f ${ref}_${query}_chr${chrom_id}.fragment.maf -R 1 -Q 1 -t 50 \n" >> run_${ref}_${query}_chr${chrom_id}.sh
    	echo "#---------------------------------------------------------------------------" >> run_${ref}_${query}_chr${chrom_id}.sh
    	echo "echo End Time is \`date\`" >> run_${ref}_${query}_chr${chrom_id}.sh
    	echo "time2=\`date +\"%Y-%m-%d %H:%M:%S\"\`" >> run_${ref}_${query}_chr${chrom_id}.sh
    	echo "timerun1=\$((\$(date +%s -d \"\$time2\") - \$(date +%s -d \"\$time1\")))" >> run_${ref}_${query}_chr${chrom_id}.sh
    	echo "echo \$timerun1 | awk '{print \"Running time is \" \$1*1/3600 \" hours(\" \$1*1/60  \" mins)\"}'" >> run_${ref}_${query}_chr${chrom_id}.sh
    	chmod 755 run_${ref}_chr${chrom_id}.sh
    	chmod 755 run_${ref}_${query}_chr${chrom_id}.sh
	done
	cd ../
done
```

## Converting MAF to VCF
**Converting the output of AnchorWave to the gVCF format using TASSEL**
```
for chrom_id in `seq 1 14`
do
	cat cat specieslist | while read query
	do
	run_pipeline.pl -Xmx100G -debug -MAFToGVCFPlugin -referenceFasta /data/yangg/project/pangenome/03.pangenome/00.genome/02.each_chrom/TW_t2/TW_t2.chr${chrom_id}.fasta -mafFile TW_t2_${query}_chr${chrom_id}.maf -sampleName TW_t2_${query}_chr${chrom_id}_anchorwave -gvcfOutput ./TW_t2_${query}_chr${chrom_id}.gvcf -fillGaps false -bgzipAndIndex false > ./TW_t2_${query}_chr${chrom_id}_outputMafToGVCF.txt
	done
done
```
** Converting gVCF to VCF format **
```
cat sample.list | while read f; do python3 convert_gvcf2vcf.py ${f}.gvcf > ${f}.vcf; sed -i 's/1:0/1\/1:0/g' ${f}.vcf; python3 filter_length.py ${f}.vcf > ${f}_longerthan10bp.vcf; bgzip -@ 50 ${f}_longerthan10bp.vcf; tabix -C -p vcf ${f}_longerthan10bp.vcf.gz; done
```

## Variants merge
** Merging VCF files using BCFtools **
```
for i in `seq 1 14`
do
	bcftools merge ../02.gvcf2vcf/TW_t2_TW_svevo_chr${i}_longerthan10bp.vcf.gz ../02.gvcf2vcf/TW_t2_TW_t15_chr${i}_longerthan10bp.vcf.gz ../02.gvcf2vcf/TW_t2_TW_t21_chr${i}_longerthan10bp.vcf.gz ../02.gvcf2vcf/TW_t2_TW_t5_chr${i}_longerthan10bp.vcf.gz ../02.gvcf2vcf/TW_t2_TW_t9_chr${i}_longerthan10bp.vcf.gz ../02.gvcf2vcf/TW_t2_TW_1905_chr${i}_longerthan10bp.vcf.gz ../02.gvcf2vcf/TW_t2_TW_1954_chr${i}_longerthan10bp.vcf.gz ../02.gvcf2vcf/TW_t2_TW_1097_chr${i}_longerthan10bp.vcf.gz ../02.gvcf2vcf/TW_t2_TW_zavitan_chr${i}_longerthan10bp.vcf.gz | bcftools norm -m -any -N | bcftools norm -d none --fasta-ref ../../01.genome/TW_t2.chr${i}.fasta | bcftools sort > chr${i}_longerthan10bp.merged.vcf
	python3 re-chrom_position.py chr${i}_longerthan10bp.merged.vcf > chr${i}_longerthan10bp.merged.split.vcf
	bcftools annotate --rename-chrs change_chrom.txt chr${i}_longerthan10bp.merged.split.vcf --force -Ov > chr${i}_longerthan10bp.merged.split.changechr.vcf
	bgzip -@ 50 chr${i}_longerthan10bp.merged.split.changechr.vcf
	tabix -C -p vcf chr${i}_longerthan10bp.merged.split.changechr.vcf.gz
done
```
** Spliting chromosomes **
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
cat bed | while read f
do
	c1=$(echo ${f} | awk '{print$1}')
	c2=$(echo ${f} | awk '{print$2}')
	bcftools view -r ${c2} chr${c1}_longerthan10bp.merged.split.changechr.vcf.gz > chr${c2}.vcf
done
```

## Merging and optimizing the multiallelic SVs using Panpop
```
for i in `seq 1 28`
do
	mkdir panpop_chr${i}_output
	perl PART_run.pl -i chr${i}.vcf -o panpop_chr${i}_output -r TW_t2_split.chr${i}.fasta -t 50 --tmpdir ./tmp
done
```
** Sorting chromosomes ID according to ASCII**
```
gunzip -c panpop_chr1_output/3.final.vcf.gz | head -16  > TWGG.vcf
for i in 10 11 12 13 14 15 16 17 18 19 2 20 21 22 23 24 25 26 27 28 3 4 5 6 7 8 9
do
gunzip -c panpop_chr${i}_output/3.final.vcf.gz | grep "contig" >> TWGG.vcf
done
gunzip -c panpop_chr1_output/3.final.vcf.gz | sed -n '17,28p' >> TWGG.vcf

for i in 1 10 11 12 13 14 15 16 17 18 19 2 20 21 22 23 24 25 26 27 28 3 4 5 6 7 8 9
do
gunzip -c panpop_chr${i}_output/3.final.vcf.gz | sed -n '29,$p' >> TWGG.vcf
done
bgzip -@ 100 -c TWGG.vcf > TWGG.vcf.gz
tabix -p vcf TWGG.vcf.gz
```

## Construction of tetraploid wheat graph genome (TWGG)

vg version = 1.56.0

**VG indexing**
`vg autoindex --workflow giraffe -t 128 -r TW_t2_split.fasta -v TWGG.vcf.gz -p TWGG`

`vg snarls --threads 52 -T TWGG.giraffe.gbz > TWGG.snarls`

**Short reads mapping**
```
path=/home/yangg/my_data/project/pangenome/00.data/03.subspecies_resequencing_data/01.data_fastp/fastp_fastq
cat sample.list | while read f
do
mkdir -p ${f}
echo "#!/bin/bash
#SBATCH -o job.vg_${f}.%j.out
#SBATCH -p tcum256c128Partition
#SBATCH -J ${f}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=50" > ${f}/vg_call.sh

echo "ulimit -u 40960" >> ${f}/vg_call.sh
echo -e "vg giraffe -t 128 -Z /data/yangg/project/pangenome/03.pangenome/06.vg2/TWGG.giraffe.gbz -m /data/yangg/project/pangenome/03.pangenome/06.vg2/TWGG.min -d /data/yangg/project/pangenome/03.pangenome/06.vg2/TWGG.dist -f ${path}/${f}_1.fastp.fq.gz -f ${path}/${f}_2.fastp.fq.gz -o bam > ${f}.mapped.bam" >> ${f}/vg_call.sh
done
```

## SV calling
```
cat sample.list | while read f
do
	samtools addreplacerg -r '@RG\tID:\${f}\tSM:\${f}' ${f}.mapped.bam -o ${f}.mapped.readgroup.bam
	samtools sort ${f}.mapped.readgroup.bam > ${f}.mapped.readgroup.sorted.bam
	samtools index ${f}.mapped.readgroup.sorted.bam
    echo "id,path,depth,read length" > ${f}_for_paragraph.txt
    echo ${f},${f}".mapped.readgroup.sorted.bam,15,150" >> ${f}_for_paragraph.txt

	mkdir -p paragraph_${f}_output
    find /tmp -name "*.vcf.gz" | xargs rm
    find /tmp -name "*.csi" | xargs rm
    multigrmpy.py -m ${f}_for_paragraph.txt -i TWGG.vcf -r TW_t2_split.fasta -o paragraph_${f}_output -t 80 -M 15
done
```




