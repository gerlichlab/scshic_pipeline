import pandas as pd
import numpy as np
import seaborn as sbn
import os
import matplotlib.pyplot as plt
import pairtools
import re
import click


def parseBarcodes(path):
    """Takes input directory containing pairsam files and parses barcodes
    of the form ADAPTOR_SX, with ADAPTOR containing 6 letters and X being a number
    with one or more digits."""
    barcodes = [re.findall(r"[A-Z]{6}\_S[0-9]+", i)[0] for i in os.listdir(path)]
    return barcodes


def getInStream(path):
    """Gets inputstream of pairsfile"""
    header, pairs_body = pairtools._headerops.get_header(
        pairtools._fileio.auto_open(path, nproc=4, mode="r")
    )
    cols = pairtools._headerops.extract_column_names(header)
    return header, cols, pairs_body


def countPairs(path):
    """Iterates through .pairs file and counts
    the amount of reads"""
    header, cols, pairs_body = getInStream(path)
    count = 0
    for line in pairs_body:
        count += 1
    return count


@click.command()
@click.option(
    "--inputdir",
    help="Directory with the pipline outputs containing working directory of pipeline.",
)
@click.option("--resultsdir", help="Directory for the Results")
def qc2(inputdir, resultsdir):
    # set wd
    os.chdir(inputdir)

    # parse barcodes

    barcodes = parseBarcodes("./s4t_pairsam/")

    # count all reads

    countAll = {}

    for barcode in barcodes:
        countAll[barcode] = countPairs(f"./s4t_pairsam/{barcode}.dedup.s4t.pairsam.gz")

    # count doubly labeled

    countDouble = {}

    for barcode in barcodes:
        cis = countPairs(f"./s4t_merged_pairsam/{barcode}.cis.pairs.gz")
        trans = countPairs(f"./s4t_merged_pairsam/{barcode}.trans.pairs.gz")
        countDouble[barcode] = cis + trans

    # calculate percentages

    percDouble = {}

    for barcode in barcodes:
        percDouble[barcode] = (
            countDouble[barcode] / (countDouble[barcode] + countAll[barcode])
        ) * 100

    plotFrame = pd.DataFrame(percDouble, index=[0])
    plotMolt = pd.melt(plotFrame)

    # plot percentages

    f, ax = plt.subplots()
    sbn.barplot(x="variable", y="value", data=plotMolt, ax=ax, color="red")
    ax.set_ylabel("Perc Doubly labeled reads")
    ax.set_xlabel("Conditions")
    sbn.despine()
    locs, labels = plt.xticks()
    plt.xticks(locs, labels, rotation=90)
    f.savefig(os.path.join(resultsdir, "Perc_double.png"), bbox_inches="tight")


if __name__ == "__main__":
    qc2()
