######
#Annoate tandem repeats, use IG77365 as an example
trf  ./IG77365.pctg.fa 2 7 7 80 10 50 500 -d -l 10 -h

#Annoate TEs with EDTA for each accession;Annoate TEs, use IG77365 as an example
EDTA_raw.pl --type ltr --genome IG77365.l0.nucleus.fa  --threads 128
EDTA_raw.pl --type tir --genome IG77365.l0.nucleus.fa  --threads 128
EDTA_raw.pl --type helitron --genome IG77365.l0.nucleus.fa  --threads 128
EDTA.pl  --genome IG77365.l0.nucleus.fa  --overwrite 0 --sensitive 1 --threads 128

remask.pl IG77365.l0.nucleus.fa IG77365.l0.nucleus.fa.mod.EDTA.intact.gff3 > IG77365.fa.intactTE.masked.fa 

seqkit grep -nr -i -p Unknown IG77365.l0.nucleus.fa.mod.EDTA.TElib.fa > un.fa 
seqkit grep -nr -v -i -p Unknown IG77365.l0.nucleus.fa.mod.EDTA.TElib.fa > k.fa
terl_test.py -m DS3 -f un.fa
TE_rename.pl un.fa TERL.IG77365.l0.nucleus.fa.mod.EDTA.TElib.fa > un.terl_ds3.fa
cat k.fa un.terl_ds3.fa > all.terl_ds3.fa 

mkdir 03.EDTA.TElib.mask
RepeatMasker -nolow -no_is -norna -engine ncbi -parallel 32 -gff -dir ./edtaTE_out -lib all.terl_ds3.fa IG77365.fa.intactTE.masked.fa

repeat_to_gff.pl --prefix EDTA IG77365.fa.intactTE.masked.fa.out 
repeat_conjoin_classify.pl IG77365.fa.intactTE.masked.fa.out.gff > IG77365.fa.intactTE.masked.fa.out.nr.gff 
less IG77365.l0.nucleus.fa.mod.EDTA.intact.gff3 | awk '$3!="repeat_region" && $3!="target_site_duplication" && $3!="long_terminal_repeat"' >  IG77365.l0.nucleus.fa.mod.EDTA.intact.gff3.TE.gff
repeat_conjoin_classify.pl IG77365.l0.nucleus.fa.mod.EDTA.intact.gff3.TE.gff > IG77365.l0.nucleus.fa.mod.EDTA.intact.gff3.TE.gff.nr.gff 
cat IG77365.l0.nucleus.fa.mod.EDTA.intact.gff3.TE.gff.nr.gff IG77365.fa.intactTE.masked.fa.out.nr.gff > IG77365.l0.nucleus.fa.all.TE.gff
repeat_conjoin_classify.pl IG77365.l0.nucleus.fa.all.TE.gff > IG77365.l0.nucleus.fa.all.TE.nr.gff

perl filter_small.pl IG77365.l0.nucleus.fa.all.TE.nr.gff 100 > IG77365.all.TE.100.gff
remask.pl IG77365.l0.nucleus.fa  IG77365.all.TE.100.gff soft > IG77365.l0.nucleus.fa.all.TE.soft-masked.fa.100 

stat_TE.pl -gff IG77365.all.TE.100.gff -rank all > te.stat.100
stat_TE.pl -gff IG77365.all.TE.100.gff -rank order >> te.stat.100
stat_TE.pl -gff IG77365.all.TE.100.gff -rank superfamily >> te.stat.100

gff_ctg2chr_yahs.pl IG77365.all.TE.100.gff z.IG77365.FINAL.agp  > IG77365.chr.TE.100.gff
