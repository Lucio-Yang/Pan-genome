# 3. Pan-genome analysis

## Pan-genome analysis overview:

1. Pan-genome analysis
   1. pan-genome construct
   2. Ka/Ks analysis
   3. TE density analysis
2. Evolutionary history analysis
   1. Species tree construction
   2. Divergence time estimation
   3. Gene family analysis
   4. GO enrichment
3. NLR-pan-genome analysis
4. A and B subgenomes compare analysis

## Pan-genome analysis

### Pan-genome construct

*dagchainer    ——**yangguang** shige*

```text
commonds
```

#### Ka/Ks_analysis

To calculate the Ka/Ks values for analyzing the select stress analysis among the gene pairs.

```text
python run_ks.py
```

#### TE density analysis

To statistic the density of transposable elements (TEs) in the gene ontology and its upstream and downstream 2,000 bp regions.

```text
perl genebody.and.down.up.stream.TEs.count.pl
```

### Evolutionary history analysis

Analyzing the evolutionary history of tetraploid wheats. 

#### Species tree construction

Based on the single-copy syntenic groups, constructing the species tree. 

```text
#extracte the single-copy syntenic groups (SGs)
awk -F'\t' '$2 == 13 && $3 == 13 {print $1}' sample.list.SG.pan > SingleCopySGs.txt
awk -F'\t' 'BEGIN{OFS="\t"} {$3=$4=$5=""; $1=$1; gsub(/^\t+|\t+$/, ""); gsub(/\t+/, "\t"); print}' sample.list.SG.pan > sample.list.SG.pan.modified

#consruct the maximum likelihood gene trees of all single-copy SGs
python 1extract_cds.py 
python 2run_iqtree.py
python 3change_name.py

#construct the specise tree
java -jar astral.5.7.8.jar --outgroup Hv -i all.reroot.tre -o species_tree.tree -t 3
```

#### Divergence time estimation 

Estimating the divergence time

```text
bash SingleCopySGs_alignment.sh
mcmctree mcmctree.ctl
```

#### Gene family analysis

Analyzing the expansion and construction of gene families among each node of the evolutionary history.

```text
orthofinder -f input_last -M msa -S diamond -t 96 -a 8 -I 1.5 -T iqtree -A mafft
cafe5 -i gene_families_filter.txt -t tree.txt -p -k 2 -o k2p
```

#### Go enrichment

GO_enrichment based on the gene annotation by interproscan.

```text
python3 interproscan_GO_Ath.py
```

### NLR-pan-genome analysis

Identificating and classifying the NLR genes, and extract the NLR-pangenome.

```text
#1. Identificating and classifying the NLR genes.
cat sample.list | while read f
do
java -jar ./NLR-Annotator/NLR-Annotator-v2.1b.jar -i ./genome/${f}.fa -x ./NLR-Annotator/src/mot.txt -y ./NLR-Annotator/src/store.txt -o ${f}.NLR.txt -g ${f}.NLR.gff -b ${f}.NLR.gff -t 20
done
#2. corresponding to the genes based on the GFF3 and ${f}.NLR.bed files.
python NLR.geneid.py
#3. extract the NLR-pangenome
bash NLR-pangenome.sh
```

### A and B subgenomes compare analysis

*RBH    ——**yangguang** shige*

```text
commonds
```
