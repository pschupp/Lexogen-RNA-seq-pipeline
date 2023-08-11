# Usage

Each `snakemake` script can be run the usual way with Snakemake in its own virtual environment via:

```{bash}
# for dry run
#   -n : dry run
#   -r : reason for each excecuted rule
#   -s : snakemake file
snakemake -n -r -s rna-seq-lexogen-3p.smk --configfile human-config.yaml
# create DAG diagram
snakemake --dag -s rna-seq-lexogen-3p.smk --configfile human-config.yaml | dot -Tsvg > dag.svg
# create DAG diagram
snakemake --rulegraph -s rna-seq-lexogen-3p.smk --configfile human-config.yaml | dot -Tsvg > rules.svg
# for real run
snakemake --cores 15 -s rna-seq-lexogen-3p.smk --configfile human-config.yaml
```

# Dependencies

Python packages:
- `os`
- `numpy`
- `pandas`
- `re`
- `math`

