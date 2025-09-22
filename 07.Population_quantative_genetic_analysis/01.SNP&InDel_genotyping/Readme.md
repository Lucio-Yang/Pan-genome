## WGS mapping<br>
**Mapping the WGS data of each accession to Kronos reference**
```
ref=kronos.final.split.fa

cat sample.id | while read f
do

mkdir -p ${f}

echo -e "#!/bin/bash
#SBATCH -o job.lush_${f}.%j.out
#SBATCH -J ${f}_lush
#SBATCH --nodes=1
#SBATCH -p tcum256c128Partition
#SBATCH --ntasks-per-node=128" > ${f}/run_${f}.sh

echo -e "echo Start Time is \`date\`"  >> ${f}/run_${f}.sh

echo -e "./bwa-mem2-2.2.1_x64-linux/bwa-mem2.avx2 mem -t 80 -R '@RG\\\tID:${f}\\\tPL:Illumina\\\tSM:${f}' ${ref} ${f}_1.fastp.fq.gz ${f}_2.fastp.fq.gz > ${f}.sam" >> ${f}/run_${f}.sh
echo -e "samtools view -@ 128 -bS ${f}.sam > ${f}.bam" >> ${f}/run_${f}.sh
echo -e "samtools collate -@ 128 -O -u ${f}.bam -T ${f}.collate | samtools fixmate -@ 128 -u -m - - | samtools sort -@ 128 -o ${f}.sorted.bam" >> ${f}/run_${f}.sh
echo -e "samtools index -c -@ 128 ${f}.sorted.bam" >> ${f}/run_${f}.sh
echo -e "samtools markdup -@ 128 -r -s -f ${f}.sorted.rmdup.metrics -O BAM ${f}.sorted.bam ${f}.sorted.rmdup.bam" >> ${f}/run_${f}.sh
echo -e "samtools index -@ 128 -c ${f}.sorted.rmdup.bam" >> ${f}/run_${f}.sh
echo -e "/home/xud/software/LUSH-DNASeq-pipeline-main/bin/LUSH_toolkit-HC/lush_hc HaplotypeCaller --nthreads 50 -R ${ref} -I ${f}.sorted.rmdup.bam --pcr-indel-model NONE --base-quality-score-threshold 20 --emit-ref-confidence GVCF -O ${f}_variant.g.vcf.gz" >> ${f}/run_${f}.sh

echo -e "echo End Time is \`date\`" >> ${f}/run_${f}.sh

done
```
**Genotyping for all samples**
```
#!/usr/bin/bash
echo "#!/bin/bash
#SBATCH -o job.gtx.%j.out
#SBATCH -J gtx
#SBATCH --nodes=1
#SBATCH -p fat1t128c
#SBATCH --ntasks-per-node=128" > run_gtx_gvcf2vcf2.sh
echo "export RHWL_LICENSE_SERVER=12.12.12.101:9000" >> run_gtx_gvcf2vcf2.sh
echo -e "/home/apps/gtx/bin/gtx joint -r ./kronos.final.split.fa" >> run_gtx_gvcf2vcf2.sh
cat sample_gvcf.list | while read f
do
echo -e "	-v ${f}_variant.g.vcf.gz"" \\" >> run_gtx_gvcf2vcf2.sh
done
echo " -o output.vcf.gz -t 128" >> run_gtx_gvcf2vcf2.sh
```
**Filtering for high-quality SNP and InDel variants**
```
########## extract SNP and InDel
genome=kronos.final.split.fa
cat ./chrom.id | while read f
do
mkdir -p ${f}
echo -e "#!/bin/bash
#SBATCH -o job.vcffilter.${f}.%j.out
#SBATCH -p fat1t128c
#SBATCH -J ${f}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128" > ${f}/01.vcffiltering_${f}.sh
vcf=TW.sample742
gatk=~/software/gatk/bin/gatk
echo -e "$gatk  --java-options \"-Xmx80g -Djava.io.tmpdir=./tmp\"  SelectVariants -V ../${vcf}.vcf.gz --select-type SNP -L ${f} -O ${vcf}.${f}.raw.snp.vcf.gz" >> ${f}/01.vcffiltering_${f}.sh
echo -e "$gatk  --java-options \"-Xmx80g -Djava.io.tmpdir=./tmp\"  SelectVariants -V ../${vcf}.vcf.gz --select-type INDEL -L ${f} -O ${vcf}.${f}.raw.indel.vcf.gz" >> ${f}/01.vcffiltering_${f}.sh
echo -e "$gatk --java-options "-XX:ParallelGCThreads=128" VariantFiltration \\
-V ${vcf}.${f}.raw.snp.vcf.gz \\
--filter-expression \"QD < 2.0\"    \\
--filter-name \"low_QD\"            \\
--filter-expression \"SOR > 3.0\"   \\
--filter-name \"high_SOR\"          \\
--filter-expression \"MQ < 40.0\"   \\
--filter-name \"low_MQ\"            \\
--filter-expression \"FS > 60.0\"   \\
--filter-name \"high_FS\"           \\
--filter-expression \"MQRankSum < -12.5\" \\
--filter-name \"low_MQRankSum\"     \\
--filter-expression \"ReadPosRankSum < -8.0\" \\
--filter-name \"low_ReadPosRankSum\" \\
-O ${vcf}.${f}.filter.snp.vcf.gz" >> ${f}/01.vcffiltering_${f}.sh
echo -e "$gatk --java-options "-XX:ParallelGCThreads=128" VariantFiltration \\
-V ${vcf}.${f}.raw.indel.vcf.gz \\
--filter-expression \"QD < 2.0\"    \\
--filter-name \"low_QD\"            \\
--filter-expression \"SOR > 10.0\"   \\
--filter-name \"high_SOR\"          \\
--filter-expression \"FS > 200.0\"   \\
--filter-name \"high_FS\"           \\
--filter-expression \"MQRankSum < -12.5\" \\
--filter-name \"low_MQRankSum\"     \\
--filter-expression \"ReadPosRankSum < -8.0\" \\
--filter-name \"low_ReadPosRankSum\" \\
-O ${vcf}.${f}.filter.indel.vcf.gz" >> ${f}/01.vcffiltering_${f}.1.sh
echo -e "$gatk  --java-options \"-Xmx80g -Djava.io.tmpdir=./tmp\" SelectVariants -V ${vcf}.${f}.filter.snp.vcf.gz --exclude-filtered -O ${vcf}.${f}.filtered.snp.vcf.gz" >> ${f}/01.vcffiltering_${f}.sh
echo -e "$gatk  --java-options \"-Xmx80g -Djava.io.tmpdir=./tmp\" SelectVariants -V ${vcf}.${f}.filter.indel.vcf.gz --exclude-filtered -O ${vcf}.${f}.filtered.indel.vcf.gz" >> ${f}/01.vcffiltering_${f}.sh
done
```
**Filtering for missing rate and minor allele frequency of SNP and InDel variants**
```
cat ./chrom.id | while read f
do
echo -e "#!/bin/bash
#SBATCH -o job.softfilter.${f}.%j.out
#SBATCH -p m256Partition
#SBATCH -J ${f}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10" > ./${f}/02.softfilter_${f}.sh
echo -e "python3 re_chrom.py genome_split.bed TW.sample742.${f}.filtered.snp.vcf.gz | bgzip -@ 128 -c > TW.sample742.${f}.filtered.snp.chr.vcf.gz" >> ./${f}/02.softfilter_${f}.sh
echo -e "python3 re_chrom.py genome_split.bed TW.sample742.${f}.filtered.indel.vcf.gz | bgzip -@ 128 -c > TW.sample742.${f}.filtered.indel.chr.vcf.gz" >> ./${f}/02.softfilter_${f}.sh
echo -e "tabix -C -p vcf TW.sample742.${f}.filtered.snp.chr.vcf.gz" >> ./${f}/02.softfilter_${f}.sh
echo -e "tabix -C -p vcf TW.sample742.${f}.filtered.indel.chr.vcf.gz" >> ./${f}/02.softfilter_${f}.sh
echo -e "vcftools --gzvcf TW.sample742.${f}.filtered.snp.chr.vcf.gz --max-missing 0.8 --min-alleles 2 --max-alleles 2 --maf 0.05 --recode --recode-INFO-all --stdout | bgzip -@ 50 -c > TW.sample742.${f}.filtered.snp.chr.miss0.8.maf0.05.vcf.gz" >> ./${f}/02.softfilter_${f}.sh
echo -e "vcftools --gzvcf TW.sample742.${f}.filtered.indel.chr.vcf.gz --max-missing 0.8 --min-alleles 2 --max-alleles 2 --maf 0.05 --recode --recode-INFO-all --stdout | bgzip -@ 50 -c > TW.sample742.${f}.filtered.indel.chr.miss0.8.maf0.05.vcf.gz" >> ./${f}/02.softfilter_${f}.sh
echo -e "tabix -C -p vcf TW.sample742.${f}.filtered.snp.chr.miss0.8.maf0.05.vcf.gz" >> ./${f}/02.softfilter_${f}.sh
echo -e "tabix -C -p vcf TW.sample742.${f}.filtered.indel.chr.miss0.8.maf0.05.vcf.gz" >> ./${f}/02.softfilter_${f}.sh
echo -e "bcftools annotate --rename-chrs change_chrom.txt TW.sample742.${f}.filtered.snp.chr.miss0.8.maf0.05.vcf.gz --threads 128 -Oz > TW.sample742.${f}.filtered.snp.miss0.8.maf0.05.vcf.gz" >> ./${f}/02.softfilter_${f}.sh
echo -e "bcftools annotate --rename-chrs change_chrom.txt TW.sample742.${f}.filtered.indel.chr.miss0.8.maf0.05.vcf.gz --threads 128 -Oz > TW.sample742.${f}.filtered.indel.miss0.8.maf0.05.vcf.gz" >> ./${f}/02.softfilter_${f}.sh
echo -e "tabix -C -p vcf TW.sample742.${f}.filtered.snp.miss0.8.maf0.05.vcf.gz" >> ./${f}/02.softfilter_${f}.sh
echo -e "tabix -C -p vcf TW.sample742.${f}.filtered.indel.miss0.8.maf0.05.vcf.gz" >> ./${f}/02.softfilter_${f}.sh
done
```
