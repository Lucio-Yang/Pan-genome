## Identification for centromere of each accession
**ChIP-seq data of each accession were collected**

**Step1 Filtering adapter sequences**
```
cat sample.list | while read f
do
trimmomatic PE ${f}_1.fq.gz ${f}_2.fq.gz \
	${f}_1_paired.fq.gz ${f}_1_unpaired.fq.gz \
	${f}_2_paired.fq.gz ${f}_2_unpaired.fq.gz \
	ILLUMINACLIP:/home/YG/software/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
	SLIDINGWINDOW:4:20 MINLEN:36
fastqc ${f} -o ./
done
```

**Step2 Mapping reads to each accession**
```
cat sample.list | parallel -j 12 "bowtie2 --very-sensitive --no-mixed --no-discordant -k 10 -X 1000 -p 128 -x {} -1 {}_CENH3_1_paired.fq.gz -2 {}_CENH3_2_paired.fq.gz | samtools sort -O bam -@ 128 -o - > {}.CENH3.bam"
cat sample.list | parallel -j 12 "bowtie2 --very-sensitive --no-mixed --no-discordant -k 10 -X 1000 -p 128 -x {} -1 {}_input_1_paired.fq.gz -2 {}_input_2_paired.fq.gz | samtools sort -O bam -@ 128 -o - > {}.input.bam"

cat sample.list | while read f
do
samtools view -F 1804 -f 2 -q 30 -bS -@ 52 ${f}.input.bam > ${f}.input.filtered.bam
samtools view -F 1804 -f 2 -q 30 -bS -@ 52 ${f}.CENH3.bam > ${f}.CENH3.filtered.bam
samtools index -@ 52 -c ${f}.input.filtered.bam
samtools index -@ 52 -c ${f}.CENH3.filtered.bam
```

**Step3 Calling peak**

```
cat sample.list2 | while read f
do
genome=$(echo ${f} | awk '{print$1}')
genome_length=${echo ${f} | awk '{print$2}'}
macs3 callpeak -f BAMPE -c ${genome}_input.filtered.bam -t ${genome}.CENH3.filtered.bam -n ${genome} -g ${genome_length} --outdir ./ --broad-cutoff 0.05 -q 0.05 --min-length 1000 --seed 123
done
```

**Step4 Merging and identification potential centromere region**<br>
We merged all peaks of each accession based on the 1 Mb distance. The genome region with the highest number of peaks was considered a potential centromere region.

```
cat sample.list | while read f
do
	for i in 1A 1B 2A 2B 3A 3B 4A 4B 5A 5B 6A 6B 7A 7B
	do
		awk '$1 == "'${i}'" {print$0}' ${f}_peaks.xls | awk -v OFS='\t' '{print$1,$2,$3,"1"}' >> ${f}_peaks.bed
	done
	bedtools merge -i ${f}_peaks.bed -d 1000000 -o sum -c 4 > ${f}_peaks.centromere
done

```

**Step5 Extracting sequence of each centromere region**<br>
The file 'all_centromere_region.txt' is obtained by merging the bed files of all genomic centromeres.
```
sed -n '2,$p' all_centromere_region.txt | while read f
do
chrom=$(echo ${f} | awk '{print$1}')
start=$(echo ${f} | awk '{print$2}')
end=$(echo ${f} | awk '{print$3}')
genome=$(echo ${f} | awk '{print$8}')
samtools faidx ${genome}.fa ${chrom}:${start}-${end} > ${genome}_${chrom}.centromere.fasta
done
```

**Step6 Sequence consistency analysis of the centromeric region**
```
cat sample.list | while read f
do
	for i in 1A 1B 2A 2B 3A 3B 4A 4B 5A 5B 6A 6B 7A 7B
	do
	snakemake -s ./StainedGlass/workflow/Snakefile --configfile=./StainedGlass/config/config.yaml --config sample=${f} fasta=${f}.centromere.fasta window=10000 --cores 12 make_figures -p
	done
done
```



