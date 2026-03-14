# quickQTL

Simple marker-regression style QTL mapping. 

## Example usage

Generate a file giving F-statistics and p-values on a per site basis:

```bash
./quickQTL.py stats example-geno-table.csv example-pheno-table.csv example.out
```


Generate an Matplotlib interactive plot from the output of the command above.

```bash
./quickQTL.py plot example.out example-chrom.csv
```

Generate an PNG plot non-interactively (e.g. appropriate from a non-GUI session):

```bash
./quickQTL.py plot example.out example-chrom.csv --output example.png
```

Generate an interactive genome-wide QTL plot via Bokeh:

```bash
./bokehQTL.py genomewide example.out example-chrom.csv
```

Generate a plot for a specific QTL in Bokeh, including markers indicating annotated genes:

```bash
./bokehQTL.py chromosome -n 8 example.out example-chrom.csv example.gff
```

## Phenotype plots


Calculate mean phenotypes, conditioned on allelic state, at each variable site across the genome:

```bash
./quickQTL.py phenotype -p Phenotype_01 example-geno-table.csv example-pheno-table.csv example-phenotype-means.out
```

Plot mean phenotypes, conditioned at allelic state, across the genome:

```bash
./quickQTL.py plot-phenotype example-phenotype-means.out example-chrom.csv
```

Plot difference in mean phenotypes, conditioned at allelic state, across the genome:

```bash
./quickQTL.py plot-phenotype --diff example-phenotype-means.out example-chrom.csv
```

## Details about the genotype and phenotype data sets

Report the number of samples in the genotype and phenotype files, the size of the sample intersection set,  the number of variable sites, and the names of the phenotypes.

```bash
/quickQTL.py report example-geno-table.csv example-pheno-table.csv
```

## Customizing plot output

The `plot` sub-command include a number of options to customize the Manhattan plot output.

These include:

* `--ylim FLOAT`: Specifies maximum -log10 pvalue of the y-axis for plotting consistency.

* `--threshold FLOAT`: specify the -log10 pvalue at which to draw a horizontal line indicating the threshold for significance. By default, no threshold line is drawn. See description of `qtl_permutations.py` for estimating genome-wide significance thresholds.

* `--rotlabels`: whether chromosome names should be rotated 45 degrees in the plotting

* `--dim FLOAT (width) FLOAT (height)`: Two floating point values, specifying dimensions of the output figure in inches. Defaults to 12 inches wide and 4 inches tall.

* `-c, --colors TEXT`: Specify colors for each chromosome. Defaults to a cycling color scheme. Can be repeated multiple times with chromosomes colored in order given.  Will cycle if colors less than number of chromosomes

* `--colormap`: Specify color for each chromosome using a pre-defined [matplotlib colormap](https://matplotlib.org/stable/gallery/color/colormap_reference.html).  `--colors` argument takes precedent if both specified.

* `--marksize` and `--markeralpha`: size and alpha transparency of plotted points.

## Permutations

Genome-wide significance thresholds can be estimated by performing permutations of the phenotype data.  This is implemented in the `qtl_permutations.py` script.

Run 100 permutations of phenotype data. For each permutation keep the -log10 p-values above the specified quantile threshold (e.g. 99th percentile) and save those to a file.

```bash
./qtl_permutations.py -n 100 -q 0.99 example-geno-table.csv example-pheno-table.csv example.permutations.out
```

Estimate the genome-wide significance threshold from the permutation output file by taking the specified quantile (e.g. 99th percentile) of the maximum -log10 p-value across all permutations.

```bash
./qtl_permutations.py threshold -q 0.99 example.permutations.out example.threshold.out
```

The output of threshold is a 1-row CSV file with the following columns: `threshold_quantile`, `empirical_threshold`, `theoretical_threshold`


Plot the distribution of maximum -log10 p-values across all permutations, along with the theoretical extreme value distribution fit to those values.

```bash
./qtl_permutations.py plot example.permutations.out example.threshold.out
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

