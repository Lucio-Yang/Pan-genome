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
