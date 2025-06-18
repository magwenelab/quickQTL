#!/usr/bin/env python

import numpy as np
import scipy.stats as stats

import pandas as pd
import click


@click.command()
@click.option(
    "--threshold",
    type=click.FloatRange(0,1),
    default=0.95,
    help="Quantile of the empirical/GEV distribution to calculate as threshold.",
)
@click.argument("infile", type=click.File("r"))
@click.argument("outfile", type=click.File("w"))
def permute(infile, outfile, threshold):
    """
    Estimate log10 Pvalue QTL threshold from permuted max log10 Pvals
    """
    pvals = pd.read_csv(infile)

    empirical_threshold = pvals.max_log10Pval.quantile(threshold)
    shape, loc, scale = stats.genextreme.fit(pvals["max_log10Pval"])

    gev_threshold = stats.genextreme.isf(1.0 - threshold, shape, loc, scale)

    df = pd.DataFrame({
        "threshold_quantile": [threshold],
        "empirical_cutoff": [empirical_threshold],
        "theoretical_cutoff": [gev_threshold],
    })

    df.to_csv(outfile, index=False)
    

if __name__ == "__main__":
    permute()

