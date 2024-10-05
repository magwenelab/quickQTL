#!/usr/bin/env python

import itertools
from json import tool


import pandas as pd
import numpy as np

import bokeh.plotting as bplt
import bokeh.palettes as palletes
from bokeh.transform import factor_cmap

from bokeh.models import (
    ColumnDataSource,
    FixedTicker,
    Range1d,
    NumeralTickFormatter,
    WheelZoomTool,
    PanTool,
    HoverTool,
)

import click


@click.group()
def cli():
    pass


def gff2DF(fname):
    gff = pd.read_table(
        fname,
        comment="#",
        names=(
            "Seqid",
            "Source",
            "Type",
            "Start",
            "End",
            "Score",
            "Strand",
            "Phase",
            "Attributes",
        ),
    )
    return gff


def arrow_patch_coords(start, end, strand, width=5, level=0, rel_arrow_width=0.10):
    hw = width / 2.0
    x1, x2 = (start, end) if (strand == "+") else (end, start)
    delta = abs(end - start) * rel_arrow_width
    if strand == "+":
        head_base = max(x1, x2 - delta)
    elif strand == "-":
        head_base = min(x1, x2 + delta)
    else:
        head_base = x2
    return (
        [x1, x1, head_base, x2, head_base, x1],
        [e + level for e in [-hw, hw, hw, 0, -hw, -hw]],
    )


def chromosome_plot(statfile, chromfile, gfffile, nchrom):

    stats = pd.read_csv(statfile)
    chroms = pd.read_csv(chromfile)
    focalchrom = chroms.Chromosome[nchrom - 1]
    stats = stats[stats.Chromosome == focalchrom]
    stats["Log10Pvalue"] = -np.log10(stats.Pvalue)
    maxPval = stats.Log10Pvalue.max()


    gff = gff2DF(gfffile)
    
    gff = gff[(gff.Seqid == focalchrom) & ((gff.Type == "gene") | (gff.Type == "protein_coding_gene"))]
    gff["FtrID"] = gff.Attributes.str.extract(r"ID=([^;]+);")
    gff["Description"] = gff.Attributes.str.extract(r"description=([^;]+);")



    fig = bplt.figure(
        x_axis_label="Coordinate",
        y_axis_label="",
        sizing_mode="stretch_width",
        tools="pan,box_zoom,reset,save",
    )

    xs, ys = [], []

    for ftr in gff.itertuples():
        level = maxPval + 1 if ftr.Strand == "+" else maxPval + 2
        x, y = arrow_patch_coords(
            ftr.Start, ftr.End, ftr.Strand, level=level, width=1, rel_arrow_width=0.10
        )
        xs.append(x)
        ys.append(y)

    ftrsource = ColumnDataSource(
        dict(
            xs=xs,
            ys=ys,
            FtrType=gff.Type,
            FtrID=gff.FtrID,
            FtrDesc = gff.Description,
            start=gff.Start,
            end=gff.End,
        )
    )

    TOOLTIPS_patches = [
        ("ID", "@FtrID"),
        ("Type", "@FtrType"),
        ("start:end", "@start:@end"),
        ("Desc", "@FtrDesc")
    ]
    patches = fig.patches("xs", "ys", alpha=0.5, source=ftrsource)
    fig.add_tools(HoverTool(renderers=[patches], tooltips=TOOLTIPS_patches, mode="mouse"))


    statsource = ColumnDataSource(
        dict(
            xs=stats.Coordinate,
            ys=stats.Log10Pvalue,
            Chromosome = stats.Chromosome,
            Coordinate = stats.Coordinate
        )
    )
    TOOLTIPS_points = [
        ("Chrom", "@Chromosome"),
        ("Coord", "@Coordinate"),
    ]
    points = fig.scatter(x="xs", y="ys", color="firebrick", alpha=0.5, source=statsource)
    fig.add_tools(HoverTool(renderers=[points], tooltips=TOOLTIPS_points, mode="mouse"))


    # p.x_range = Range1d(mid - 10000, mid + 10000, bounds=(0, stats.Coordinate.max()))
    # p.x_range = Range1d(mid - 10000, mid + 10000, bounds=(0, stats.Coordinate.max()))
    fig.x_range.bounds = (0, stats.Coordinate.max())

    # p.y_range = Range1d(-20, 50)
    fig.xaxis[0].formatter = NumeralTickFormatter(format="0,0")

    fig.add_tools(WheelZoomTool(dimensions="width"))

    pantool = fig.select_one({"type": PanTool})
    pantool.dimensions = "width"

    bplt.show(fig)


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
        x_axis_label="Chromosome",
        y_axis_label="-log10(p-value)",
        sizing_mode="stretch_width",
        tooltips=TOOLTIPS,
    )

    p.xaxis.formatter = NumeralTickFormatter(format="0,0")
    p.xaxis.ticker = FixedTicker(ticks=chromsdf.PlotMidpoint.values, num_minor_ticks=0)
    p.xaxis.major_label_overrides = dict(
        zip(chromsdf.PlotMidpoint.values, chromsdf.Chromosome)
    )

    cmapper = factor_cmap(
        field_name="Chromosome",
        palette=list(palletes.Dark2_8) * 10,  # palletes.Dark2_8,
        factors=chromsdf.Chromosome.unique(),
    )

    p.scatter(x="PlotCoordinate", y="Log10Pvalue", color = cmapper, source=src)



    bplt.show(p)


@cli.command()
@click.argument("statfile", type=click.File("r"))
@click.argument("chromfile", type=click.File("r"))
def genomewide(statfile, chromfile):
    stats, chroms, src = create_data(statfile, chromfile)
    manhattan_plot(stats, chroms, src)


@cli.command()
@click.option("--nchrom", "-n", type=int, default=1)
@click.argument("statfile", type=click.File("r"))
@click.argument("chromfile", type=click.File("r"))
@click.argument("gfffile", type=click.File("r"))
def chromosome(statfile, chromfile, gfffile, nchrom):
    chromosome_plot(statfile, chromfile, gfffile, nchrom)


if __name__ == "__main__":
    cli()
