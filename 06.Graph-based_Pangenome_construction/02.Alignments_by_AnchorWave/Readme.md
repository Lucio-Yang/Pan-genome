## Alignments by [AnchorWave](https://github.com/baoxingsong/AnchorWave)
 ```
 #!/bin/bash
 
 chrom_number=14
 ref=Kronos
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