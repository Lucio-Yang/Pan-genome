## Population sample SV calling

To identify SVs in the 742 accessions, Illumina short reads were mapped to the tetraploid wheat graph pangenome using the ‘vg giraffe’ pipeline with default parameters. SVs were genotyped using Paragraph with default parameters.

**Short reads mapping**
```
#!/usr/bin/bash

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

**Genotyping**
```
#!/usr/bin/bash

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
    multigrmpy.py -m ${f}_for_paragraph.txt -i TWGG.vcf -r Kronos_split.fasta -o paragraph_${f}_output -t 80 -M 15
done
```