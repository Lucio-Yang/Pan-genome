#### Pan-genome construct<br>

SG were identified using the Dagchainer and SynPan tool. Input files is the gene bed file and fasta file (amino acid sequence) of each accession. The file 'sample.list' is the prefix of bed and fasta of each accession.

```text
./synpan_build.pl sample.list
```
Saturation analysis of core genes and pan-genes
```
for f in `seq 2 12`
do
	for i in `seq 1 100`
	do
	python SG_step1.py sample.list sample.list.SG.pan ${f} > temp${f}
	python SG_step2.py temp${f} ${f} >> sample_${f}.stat
	done
done
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

### 
