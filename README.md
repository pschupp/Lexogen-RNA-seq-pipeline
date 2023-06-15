# How to use this repo
Each ```snakemake``` script can be run the usual way with Snakemake in its own virtual environment via:
```{bash}
# for dry run
#   -n : dry run
#   -r : reason for each excecuted rule
#   -s : snakemake file
snakemake -n -r -s [name of file.smk] 
# create DAG diagram
snakemake --dag -s [name of file.smk]| dot -Tsvg > dag.svg
# create DAG diagram
snakemake --rulegraph -s [name of file.smk]| dot -Tsvg > rules.svg
# for real run
snakemake --cores 15 -s [name of file.smk]
```
