## Construction of tetraploid wheat graph genome (TWGG)

vg version = 1.56.0  

**VG indexing**

`vg autoindex --workflow giraffe -t 128 -r Kronos_split.fasta -v TWGG.vcf.gz -p TWGG`

`vg snarls --threads 52 -T TWGG.giraffe.gbz > TWGG.snarls`