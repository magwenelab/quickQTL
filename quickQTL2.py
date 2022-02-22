#!/usr/bin/env python

from itertools import cycle

import numpy as np
from scipy import stats
import pandas as pd
from matplotlib import pyplot as plt

import click

from pandarallel import pandarallel

pandarallel.initialize(verbose=0)


@click.group()
def cli():
    pass


def load_files(genotypefile, phenotypefile):
    with open(genotypefile, "r") as gfile, open(phenotypefile, "r") as pfile:

        genotypes_df = pd.read_csv(gfile, index_col=("Chromosome", "Coordinate"))
        phenotypes_df = pd.read_csv(pfile, index_col="Sample_Name")

        # reindexing makes phenotype sample order conform to genotype sample order
        # sample names in genotype header that don't exist get a NaN in the
        # phenotype columns; sample names in phenotypes that don't exist in
        # genotype header are dropped
        phenotypes_df = phenotypes_df.reindex(genotypes_df.columns)

        return genotypes_df, phenotypes_df


assoc_test_dict = {
    "anova": stats.f_oneway,
    "alexander": stats.alexandergovern,
    "kruskal": stats.kruskal,
}


def assoc_by_genotype(genorow, phenodf, focalpheno, genos=(0, 1), test="anova"):
    groups = [
        phenodf[genorow == i][focalpheno].dropna().values.flatten() for i in genos
    ]
    statistic, pval = assoc_test_dict[test](*groups)
    return statistic, pval


def assoc_test_pandarallel(genodf, phenodf, focalpheno, genos=(0, 1), test="anova"):
    assoc = genodf.parallel_apply(
        assoc_by_genotype,
        args=(phenodf, focalpheno, genos, test),
        axis=1,
        result_type="expand",
    )
    assoc.columns = ("statistic", "Pvalue")
    return assoc


def assoc_test_pandas(genodf, phenodf, focalpheno, genos=(0, 1), test="anova"):
    assoc = genodf.apply(
        assoc_by_genotype,
        args=(phenodf, focalpheno, genos, test),
        axis=1,
        result_type="expand",
    )
    assoc.columns = ("F", "Pvalue")
    return assoc


@cli.command()
@click.option("--show", type=bool, default=True)
@click.argument("assocfile", type=click.File("w"))
@click.argument("chromfile", type=click.File("w"))
def plot(statfile, chromfile, show):
    association_stats = pd.read_csv(statfile)
    chroms = pd.read_csv(chromfile)

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
            -np.log10(grp.Pvalue),
            markersize=2,
            marker="o",
            linestyle="None",
            alpha=0.75,
            color=next(colorcycle),
        )
    ax.set_xticks(chroms.Midpoint, labels=chroms.index)
    ax.set_xticks(chroms.Endpoint, minor=True)

    if show:
        plt.show()

    return fig, ax


@cli.command()
@click.option("--phenoname", "-p", type=str, default="")
@click.option("--test", type=click.Choice(assoc_test_dict.keys()), default="anova")
@click.argument("genotypefile", type=click.Path(exists=True))
@click.argument("phenotypefile", type=click.Path(exists=True))
@click.argument("outfile", type=click.File("w"))
def stats(genotypefile, phenotypefile, outfile, phenoname, test):
    """
    Calculate association statistics between genotype and phenotype.

    GENOTYPES, PHENOTYPES, and CHROMOSOMES are CSV files. See README for formatting.

    OUTFILE is the name of the CSV file to write that contains p-values.
    A dash ('-') can be substituted OUTFILE to write to stdout.
    """
    genodf, phenodf = load_files(genotypefile, phenotypefile)
    if phenoname == "":
        phenoname = phenodf.columns[0]
    assoc = assoc_test_pandarallel(genodf, phenodf, phenoname, genos=(0, 1), test=test)
    assoc.to_csv(outfile)


if __name__ == "__main__":
    cli()
