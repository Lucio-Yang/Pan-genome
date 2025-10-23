## genome survey:Estimate genome size using Illumina short reads and HIFI reads

`jellyfish count -C -m 51 -s 40G -t 20 -g raw_generators -G 3 -o IG77365_K51_raw.jf`

`jellyfish histo -h 100000 -t 128 -v -o IG77365_K51_raw.histo IG77365_K51_raw.jf`
`jellyfish stats -v -o IG77365_K51_raw.stat IG77365_K51_raw.jf  
genomescope2 -i IG77365_K51_raw.histo -k 51 -o t5.out`

`library("findGSE")`
`findGSE(histo="./IG77365_K51_raw.histo",sizek=51,outdir="IG77365")`





######
#genome survey:Estimate genome size using Illumina short reads and HIFI reads

jellyfish count -C -m 51 -s 40G -t 20 -g raw_generators -G 3 -o IG77365_K51_raw.jf  
jellyfish histo -h 100000 -t 128 -v -o IG77365_K51_raw.histo IG77365_K51_raw.jf  
jellyfish stats -v -o IG77365_K51_raw.stat IG77365_K51_raw.jf  
genomescope2 -i IG77365_K51_raw.histo -k 51 -o t5.out  

library("findGSE")
findGSE(histo="./IG77365_K51_raw.histo",sizek=51,outdir="IG77365")

######  
#assembly using hifiasm  
hifiasm -o IG77365_l0 -t 128 -l 0 --primary IG77365.all.hifi.fq.gz  > asm.all.out 2> asm.all.err  

#filter organelle contig  
minimap2 -t 52 --secondary=no -cx asm20 organ.fa pctg.fa > IG77365_to_organ.paf 2> map.log  
paf_qur_stat.pl IG77365_to_organ.paf > IG77365_to_organ.paf.stat  
awk '$3>0.95 && $4<0.05' IG77365_to_organ.paf.stat > organ.stat  
awk '{print $2}' organ.stat | numberStat.pl > organ.stat.stat   
contig_select.pl -r organ.stat pctg.fa > IG77365.l0.nucleus.fa  
fastaDeal.pl -attr id:len IG77365.l0.nucleus.fa > IG77365.l0.nucleus.fa.len  
awk '{print $2}' IG77365.l0.nucleus.fa.len | numberStat.pl > IG77365.l0.nucleus.fa.len.stat  

######  
#scaffolding  
#HIC-Pro
HiC-Pro -c config-hicpro.txt -o Analysis -i HIC_read  

#YAHS
bwa index  pctg.fa 2> index.log  
bwa mem -SP5M -t 128 IG77365.l0.nucleus.fa  IG77365_good_1.fq.gz IG77365_good_2.fq.gz  | samtools view -bo IG77365.hic.bam -  
samtools faidx IG77365.l0.nucleus.fa  
yahs -e GATC --no-contig-ec  --no-scaffold-ec IG77365.l0.nucleus.fa IG77365.hic.bam 2> yahs.log  
juicer pre -a -o out yahs.out.bin yahs.out_scaffolds_final.agp IG77365.l0.nucleus.fa.fai > out_JBAT.log 2> pre.log  
java -jar -Xmx200G juicer_tools.1.9.9_jcuda.0.8.jar pre out.txt out.hic <(cat pre.log  | grep PRE_C_SIZE | awk '{print $2" "$3}')  

juicebox manual correction

juicer post -o out.IG77365 out.review.assembly out.liftover.agp IG77365.l0.nucleus.fa 2> post.log  

matrix2heatmap.py hic_100000_abs.bed hic_100000.matrix 10  

######  
#near T2T genome assembly of Kronos  
#Nanopore ultral long reads basecalling and Quality control  
dorado basecaller ./dna_r10.4.1_e8.2_400bps_sup@v5.0.0  ./pod5/  --recursive -x 'cuda:all' > ./bam/ont.bam  
samtools fastq -@ 128 ont.bam > ont.fq  
NanoFilt -q 10 --length 100000 ./ont.fq | gzip > ont.100.fq.gz  
NanoStat -t 128 --fastq ./ont.100.fq.gz    -n 100.report  --outdir report  

#hifiasm HIFI data  
hifiasm -o kronos_l0 -t 96 -l 0 --primary Kronos.all.hifi.fq.gz > asm.all.out 2> asm.all.err  

#hifiasm: HIFI and ONT data assembly  
hifiasm -o kronos_l0 -t 128 -l 0 --telo-m CCCTAAA --primary Kronos.all.hifi.fq.gz  --ul Kronos.all.100.ont.fq.gz > asm.all.out 2> asm.all.err  

#verkko: HIFI and ONT data assembly  
verkko -d ./verkko_out --hifi Kronos.all.hifi.fq.gz  --no-correction --nano  Kronos.all.100.ont.fq.gz --threads 128 --local --local-memory 900 --local-cpus 1
28 2> verkko.log  

#hifiasm: ONT data assembly  
hifiasm -t 128 -k 63 --ont -o ONT.asm Kronos.all.100.ont.fq.gz > asm.all.out 2> asm.all.err  

#nextdenovo: ONT data assembly  
#Given the substantial size of the genome, we employed a Assembly by chromosome  strategy, using chromosome 1A as an example.  
minimap2 -t 52  -ax map-ont kronos.fa  Kronos.all.100.ont.fq.gz -o Map.result.sam  
awk '$1~/^@/'  Map.result.sam > header  
awk '$1~!/^@/' Map.result.sam > noMap.result.sam  
sort -t $'\t' -Vk3 --parallel=128 -o noMap.result.sam.sort noMap.result.sam  
cat header noMap.result.sam.sort > Map.result.sam.sort  
samtools faidx kronos.fa  
samtools view -@ 128 -bt kronos.fa.fai  Map.result.sam.sort > Map.result.sam.sort.bam  

samtools index -c Map.result.sam.sort.bam  

samtools view -@ 128 -h Map.result.sam.sort.bam 1A   > 1a.sam  
samtools view -@ 128 -bS 1a.sam | samtools sort - -o 1a.bam  
samtools bam2fq 1a.bam > 1a_reads.fq  

NextDenovo/nextDenovo ./run.cfg.1A  

#local assembly  
#use the ONT ultra-long reads unmapped to 14 chromosomes  
nextDenovo ./run.cfg.unchr

#merge different assemblys' results  

#Polishing  
minimap2 -ax map-hifi  --split-prefix k.ont.tmp --secondary=no  -t 128 006_27.patch.hifiasm1.fa Kronos.all.hifi.fq.gz  | samtools sort - -m 1500M --threads 128 -o lgs.sort.bam  
samtools index -c  lgs.sort.bam  

python /nextpolish2.py -g 006_27.patch.hifiasm1.fa -l lgs.sort.bam.fofn -r hifi -p 128 -sp -o polish1.fa  \

#assembly quality evaluation  
1. BUSCO  
busco -m genome -i kronos_l0.p_ctg.fa -o pctg_poa_out -l poales_odb10 -c 20 --offline  
busco -m prot -i longest_isoform.fasta  -l poales_odb10 -o busco -m  -c 20  

2. N50 length  
fastaDeal.pl -attr len kronos_l0.p_ctg.fa | numberStat.pl  

3. QV
meryl k=22 count output Kronos.all.hifi.fq.gz.meryl  Kronos.all.hifi.fq.gz 
meryl union-sum output Kronos.hifi.meryl Kronos.all.hifi.fq.gz.meryl  

4. LAI 
gt suffixerator -db genome.fa -indexname genome.fa -tis -suf -lcp -des -ssp -sds -dna  
gt ltrharvest -index genome.fa -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 85 -vic 10 -seed 20 -seqids yes > genome.fa.harvest.scn  
LTR_FINDER_parallel -seq genome.fa -threads 128 -harvest_out -size 1000000 -time 300  
cat genome.fa.harvest.scn genome.fa.finder.combine.scn > genome.fa.rawLTR.scn  
LTR_retriever -genome genome.fa -inharvest genome.fa.rawLTR.scn -threads 128  
