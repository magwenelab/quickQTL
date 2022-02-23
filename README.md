# quickQTL

Dead simple QTL mapping for haploids based on per-site t-tests. Returns a plot and a CSV file of p-values.s

## Dependencies

- Numpy
- Scipy
- Matplotlib
- Pandas
- Click
- pandarallel

All of these dependencies can be installed via Conda.


## Input file specifications


### Phenotype file

- CSV file with at least 2 columns
- Column 1: Sample Names (header `Sample_Name`)
- Columns 2... : Phenotypes (headers `Phenotype_XXXX`, `Phenotype_YYYY`, etc).

### Genotype file

- CSV file with at least 3 columns
- Column 1: Chromosome Name (header `Chromosome`)
- Column 2: Coordinate (header `Coordinate`)
- Column 3:... : Genotypes for `Sample_N` where each Sample_N header matches a Sample_N in the phenotype files

### Chromosome file
- CSV file with exactly 2 columns
- Column 2: Chromosome Name (header `Chromosome`)
- Column 3: Chromosome Length (header `Length`)

