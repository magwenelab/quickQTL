#!/usr/bin/env python


from itertools import cycle

import numpy as np
from scipy import stats
import pandas as pd
from matplotlib import pyplot as plt

import dask.dataframe as dd
import click


def calculate_association_stats(genotypes, phenotypes, focal_phenotype=None):

    phenotypes = phenotypes.set_index("Sample_Name")
    if focal_phenotype is None:
        focal_phenotype = phenotypes.columns[0]

    qtl_stats = []

    for i, site in genotypes.iterrows():
        chrom_coord_pval = site[:2].values.tolist()
        lgeno = site[2:]
        has_geno0 = lgeno.index[lgeno == 0]
        has_geno1 = lgeno.index[lgeno == 1]
        if (len(has_geno0) < 2) or (len(has_geno1) < 2):
            chrom_coord_pval.append(np.NaN)
            qtl_stats.append(chrom_coord_pval)
            continue
        pheno0 = phenotypes.loc[has_geno0, focal_phenotype]
        pheno1 = phenotypes.loc[has_geno1, focal_phenotype]
        tstat, pvalue = stats.ttest_ind(pheno0, pheno1)
        chrom_coord_pval.append(pvalue)
        qtl_stats.append(chrom_coord_pval)

    return pd.DataFrame(qtl_stats, columns=("Chromosome", "Coordinate", "p_value"))


def plot_log10pval(association_stats, chroms):
    chroms = chroms.set_index("Chromosome")
    chroms["Offset"] = [0] + chroms.Length[:-1].cumsum().values.tolist()
    chroms["Midpoint"] = (chroms.Length + 2 * chroms.Offset) / 2
    chroms["Endpoint"] = chroms.Length + chroms.Offset

    grouped_stats = association_stats.groupby("Chromosome")

    colorcycle = cycle(plt.cm.Dark2.colors)
    fig, ax = plt.subplots(1, 1)

    for chrom in chroms.index:
        grp = grouped_stats.get_group(chrom)
        ax.plot(
            grp.Coordinate + chroms.loc[chrom, "Offset"],
            -np.log10(grp.p_value),
            markersize=2,
            marker="o",
            linestyle="None",
            alpha=0.75,
            color=next(colorcycle),
        )
    ax.set_xticks(chroms.Midpoint, labels=chroms.index)
    ax.set_xticks(chroms.Endpoint, minor=True)

    return fig, ax


@click.command()
@click.argument("genotypes", type=click.File("r"))
@click.argument("phenotypes", type=click.File("r"))
@click.argument("chromosomes", type=click.File("r"))
@click.argument("outfile", type=click.File("w"))
def map_qtls(genotypes, phenotypes, chromosomes, outfile):
    """Carry out simple QTL mapping for haploid individuals.

    GENOTYPES, PHENOTYPES, and CHROMOSOMES are CSV files. See README for formatting.

    OUTFILE is the name of the CSV file to write that contains p-values.

    A dash ('-') can be substituted OUTFILE to write to stdout.
    """

    genos = dd.read_csv(genotypes)
    phenos = dd.read_csv(phenotypes)
    chroms = dd.read_csv(chromosomes)

    genos.drop(set(genos.columns[2:]) - set(phenos.Sample_Name), axis=1, inplace=True)

    # phenos.drop(set(phenos.Sample_Name) - set(genos.columns[2:]), axis=0, inplace=True)

    qtl_stats = calculate_association_stats(genos, phenos, None)
    qtl_stats.to_csv(outfile, index=False)

    fig, ax = plot_log10pval(qtl_stats, chroms)
    plt.show()


if __name__ == "__main__":
    map_qtls()
