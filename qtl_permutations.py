#!/usr/bin/env python

import numpy as np
import scipy.stats as stats
import pandas as pd
from matplotlib import pyplot as plt

import click


import quickQTL as qq


@click.group()
def cli():
    pass


def permuted_assoc(outfile, genodf, phenodf, focalpheno, genos=(0, 1), test="anova", n = 10, q = 0.99):
    """
    Generate n random permutations of the focal phenotype, and for each permutation, calculate the association test statistics for all sites, determine the log10 Pvalue threshold corresponding to the q quantile of the distribution of log10 Pvalues across all sites, and the extreme log10 Pvalues above the threshold to the output file.
    """
    for i in range(n):
        phenodf[focalpheno] = np.random.permutation(phenodf[focalpheno])
        tstats = qq.assoc_test_pandarallel(genodf, phenodf, focalpheno, genos, test)
        qthresh = float(tstats.log10Pvalue.quantile(q))
        above_threshold = tstats.log10Pvalue[tstats.log10Pvalue >= qthresh]
        df = pd.DataFrame({
            "extreme_log10Pval":above_threshold})   
        df.to_csv(outfile, index=False, mode="a", header=False)
    


@cli.command()
@click.option(
    "--phenoname",
    "-p",
    type=str,
    default="",
    help=(
        "Name of the phenotype column to analyze. "
        "Defaults to first column other than Sample_Name."
    ),
)
@click.option(
    "--genos",
    type=click.Choice(["01", "012", "-11"]),
    default="01",
    help="Genotype states that are used for analysis.",
)
@click.option(
    "--test",
    type=click.Choice(qq.assoc_test_dict.keys()),
    default="anova",
    help="Which statistical test of mean differences to apply. Defaults to ANOVA.",
)
@click.option(
    '-n',
    type=int,
    default=100,
    help="Number of permutations to run for estimating EVD of log10 Pvals. Default=100"
)
@click.option(
    '-q',
    type=click.FloatRange(0,1),
    default=0.99,
    help="Quantile threshold for EVD of log10 Pvals. Default=0.99"
)
@click.argument("genotypefile", type=click.Path(exists=True))
@click.argument("phenotypefile", type=click.Path(exists=True))
@click.argument("outfile", type=click.File("w"))
def generate(genotypefile, phenotypefile, outfile, phenoname, genos, test, n, q):
    """
    Generate n random permutations of the focal phenotype, and for each permutation, calculate the association test statistics for all sites, determine the log10 Pvalue threshold corresponding to the q quantile of the distribution of log10 Pvalues across all sites, and the extreme log10 Pvalues above the threshold to the output file.
    """
    genodf, phenodf = qq.load_files(genotypefile, phenotypefile)
    if phenoname == "":
        phenoname = phenodf.columns[0]

    if genos == "012":
        genotup = (0, 1, 2)
    elif genos == "-11":
        genotup = (-1, 1)
    else:
        genotup = (0, 1)

    permuted_assoc(outfile, genodf, phenodf, phenoname, genos=genotup, test=test, n=n, q=q)




@cli.command()
@click.option(
    "--threshold",
    "-t",
    type=click.FloatRange(0,1),
    default=0.99,
    help="Quantile of the empirical/GEV distribution to calculate as threshold.",
)
@click.argument("infile", type=click.File("r"))
@click.argument("outfile", type=click.File("w"))
def threshold(infile, outfile, threshold):
    """
    Estimate log10 Pvalue QTL threshold from extreme value distribution of log10 Pvals of permuted phenotypes.

    Similar to:

    Knijnenburg TA, Wessels LF, Reinders MJ, Shmulevich I. Fewer permutations, more accurate P-values. Bioinformatics. 2009 Jun 15;25(12):i161-8. doi: 10.1093/bioinformatics/btp211. PMID: 19477983; PMCID: PMC2687965.

    """
    pvals = pd.read_csv(infile, names=["extreme_log10Pval"])

    empirical_threshold = pvals.extreme_log10Pval.quantile(threshold)

    shape, loc, scale = stats.genextreme.fit(pvals.extreme_log10Pval)
    gev_threshold= stats.genextreme.isf(1.0 - threshold, shape, loc, scale)


    df = pd.DataFrame({
        "threshold_quantile": [threshold],
        "empirical_threshold": [empirical_threshold],
        "theoretical_threshold": [gev_threshold],       
    })

    df.to_csv(outfile, index=False)
     

@cli.command()
@click.argument("infile", type=click.File("r"))
def plot(infile):
    """
    Visualize the distribution of extremal log10 Pvalues from permuted phenotypes, and the corresponding PDF estimated from the GEV distribution.
    """
    pvals = pd.read_csv(infile, names=["extreme_log10Pval"])
    plt.hist(pvals.extreme_log10Pval, bins=50, density=True, 
             label="Empirical distribution")
    plt.xlabel("Extreme log10 Pvalue")
    plt.ylabel("Density")
    plt.title("Distribution of extreme log10 Pvalues from permuted phenotypes")

    # draw theoretical GEV distribution with same shape, loc, and scale parameters as the empirical distribution
    shape, loc, scale = stats.genextreme.fit(pvals.extreme_log10Pval)
    x = np.linspace(pvals.extreme_log10Pval.min(), pvals.extreme_log10Pval.max(), 100)
    pdf = stats.genextreme.pdf(x, shape, loc, scale)
    plt.plot(x, pdf, label="Theoretical GEV distribution", color="red")

    # add legend
    plt.legend()

    plt.show()



if __name__ == "__main__":
    cli()

