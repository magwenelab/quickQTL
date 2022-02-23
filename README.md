# quickQTL

Dead simple QTL mapping for haploids based on per-site test of mean differences by genotype (defaults to ANOVA). Returns a plot and a CSV file of p-values.

## Example usage

Generate a file giving F-statistics and p-values on a per site basis:

```bash
./quickQTL.py stats example-geno-table.csv example-pheno-table.csv example.out
```

Generate an Matplotlib interactive plot from the output of the command above.

```bash
./quickQTL.py plot example.out example-chrom.csv
```

Generate an PNG plot without interactive plot (e.g. appropriate from a non-GUI session):

```bash
./quickQTL.py plot example.out example-chrom.csv --output example.png
```

Generate an interactive plot via Bokeh:

```bash
./bokehQTL.py example.out example-chrom.csv
```


## Dependencies

- Numpy
- Scipy
- Matplotlib
- Pandas
- Click
- pandarallel
- Bokeh

All of these dependencies can be installed via Conda.


## Input file specifications

See the include `example-*.csv` files for exemplar versions of these file formats. 

### Phenotype file

- CSV file with at least 2 columns
- Column 1: Sample Names (header `Sample_Name`)
- Columns 2... : Phenotypes (headers `Phenotype_XXXX`, `Phenotype_YYYY`, etc).

### Genotype file

- CSV file with at least 3 columns
- Column 1: Chromosome Name (header `Chromosome`)
- Column 2: Coordinate (header `Coordinate`)
- Column 3:... : Genotypes for `Sample_N` where each Sample_N header matches a Sample_N in the phenotype files. 

**NOTE:** the order of sample headers in the genotype file does not need to match the order of sample names in the phenotype file.  There can also be samples in the phenotype file that are not in the genotype file and vice versa; those cases are ignored in the analysis.

### Chromosome file
- CSV file with exactly 2 columns
- Column 2: Chromosome Name (header `Chromosome`)
- Column 3: Chromosome Length (header `Length`)

