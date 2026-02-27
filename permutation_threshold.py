#!/usr/bin/env python

import numpy as np
import scipy.stats as stats

import pandas as pd
import click


import quickQTL





# @click.command()
# @click.option(
#     "--threshold",
#     type=click.FloatRange(0,1),
#     default=0.95,
#     help="Quantile of the empirical/GEV distribution to calculate as threshold.",
# )
# @click.argument("infile", type=click.File("r"))
# @click.argument("outfile", type=click.File("w"))
# def permute(infile, outfile, threshold):
#     """
#     Estimate log10 Pvalue QTL threshold from extreme permuted  log10 Pvals
#     """
#     pvals = pd.read_csv(infile)

#     empirical_threshold_max = pvals.max_log10Pval.quantile(threshold)
#     shape, loc, scale = stats.genextreme.fit(pvals["max_log10Pval"])
#     gev_threshold_max= stats.genextreme.isf(1.0 - threshold, shape, loc, scale)


#     empirical_threshold_quant = pvals.quant_log10Pval.quantile(threshold)
#     shapeq, locq, scaleq = stats.genextreme.fit(pvals["quant_log10Pval"])
#     gev_threshold_quant= stats.genextreme.isf(1.0 - threshold, shapeq, locq, scaleq)    

#     df = pd.DataFrame({
#         "threshold_quantile": [threshold],
#         "empirical_cutoff_max": [empirical_threshold_max],
#         "theoretical_cutoff_max": [gev_threshold_max],
#         "empirical_cutoff_quant": [empirical_threshold_quant],
#         "theoretical_cutoff_quant": [gev_threshold_quant],        
#     })

#     df.to_csv(outfile, index=False)
    

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
    Estimate log10 Pvalue QTL threshold from extreme permuted  log10 Pvals
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
     

if __name__ == "__main__":
    permute()

