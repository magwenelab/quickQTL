#!/usr/bin/env python

import pandas as pd
import numpy as np

import bokeh.plotting as bplt
import bokeh.palettes as palletes
from bokeh.transform import factor_cmap
from bokeh.models import ColumnDataSource

import click


def create_data(statfile, chromfile):
    stats = pd.read_csv(statfile)
    chroms = pd.read_csv(chromfile)  # , index_col="Chromosome")

    chroms["PlotOffset"] = [0] + chroms.Length[:-1].cumsum().values.tolist()
    chroms["PlotMidpoint"] = (chroms.Length + 2 * chroms.PlotOffset) / 2
    chroms["PlotEndpoint"] = chroms.Length + chroms.PlotOffset

    stats = stats.merge(chroms[["Chromosome", "PlotOffset"]], on="Chromosome")
    stats["PlotCoordinate"] = stats.Coordinate + stats.PlotOffset
    stats["Log10Pvalue"] = -np.log10(stats.Pvalue)

    return stats, chroms, ColumnDataSource(stats)


def manhattan_plot(statsdf, chromsdf, src):

    TOOLTIPS = [
        ("Chrom", "@Chromosome"),
        ("Coord", "@Coordinate"),
    ]

    p = bplt.figure(
        x_axis_label="Coordinate",
        y_axis_label="-log10(p-value)",
        sizing_mode="stretch_width",
        tooltips=TOOLTIPS,
    )

    p.xaxis.ticker = chromsdf.PlotMidpoint
    p.xaxis.major_label_overrides = dict(
        zip(chromsdf.PlotMidpoint, chromsdf.Chromosome)
    )

    cmapper = factor_cmap(
        field_name="Chromosome",
        palette=palletes.Dark2_6,
        factors=chromsdf.Chromosome.unique(),
    )

    p.circle(x="PlotCoordinate", y="Log10Pvalue", color=cmapper, source=src)
    bplt.show(p)


@click.command()
@click.argument("statfile", type=click.File("r"))
@click.argument("chromfile", type=click.File("r"))
def cli(statfile, chromfile):
    stats, chroms, src = create_data(statfile, chromfile)
    manhattan_plot(stats, chroms, src)


if __name__ == "__main__":
    cli()
