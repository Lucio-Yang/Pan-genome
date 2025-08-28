## Identifying homologoues between A and B genomes

**Running GeneTribe software for each sample**

```
cat sample.list | while read f
do
	mkdir -p ${f}; cd ${f}
	genetribe core -l ${f}_A -f ${f}_B -n 128
	cd ../
done
```