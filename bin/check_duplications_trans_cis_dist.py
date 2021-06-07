import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sbn
import click


# define functions


def getUnique(df):
    """Calcualtes percentage of unique reads"""
    return (
        df.loc[df[0] == "total_nodups", 1].values[0]
        / df.loc[df[0] == "total_mapped", 1].values[0]
    )


def getCisTrans(df):
    """Calcualtes percentage of cis and trans chromosomal reads"""
    cis = (
        df.loc[df[0] == "cis", 1].values[0]
        / df.loc[df[0] == "total_nodups", 1].values[0]
    )
    trans = (
        df.loc[df[0] == "trans", 1].values[0]
        / df.loc[df[0] == "total_nodups", 1].values[0]
    )
    return {"cis": cis, "trans": trans}


def getCisDist(df):
    """Calcualtes percentage of reads at certain distance"""
    result = {}
    for dist in [
        "cis_1kb+",
        "cis_2kb+",
        "cis_4kb+",
        "cis_10kb+",
        "cis_20kb+",
        "cis_40kb+",
    ]:
        result[dist] = (
            df.loc[df[0] == dist, 1].values[0]
            / df.loc[df[0] == "pair_types/UU", 1].values[0]
        )
    return result


@click.command()
@click.option("--inputdir", help="Directory with the stats files.")
@click.option("--resultsdir", help="Directory for the Results")
def qc1(inputdir, resultsdir):
    # set wd
    os.chdir(inputdir)

    # load in stats

    stats = {
        i.split("_")[0]: pd.read_csv(i, sep="\t", header=None)
        for i in os.listdir()
        if "stats" in i
    }

    # get unique reads

    unique = {barcode: getUnique(frame) for barcode, frame in stats.items()}
    uniqueF = pd.DataFrame(unique, index=[0])
    uniqueMolt = pd.melt(uniqueF)

    # plot unique reads
    f, ax = plt.subplots()
    sbn.barplot(x="variable", y="value", data=uniqueMolt, color="red")
    ax.set_xlabel("Barcodes")
    ax.set_ylabel("Fraction Unique")
    ax.set_title("Unique Reads")
    f.set_size_inches(10, 10)
    locs, labels = plt.xticks()
    plt.xticks(locs, labels, rotation=90)
    sbn.despine()
    f.savefig(os.path.join(resultsdir, "Unqiue_Reads.png"), bbox_inches="tight")

    # get cisTrans

    cisTrans = {barcode: getCisTrans(frame) for barcode, frame in stats.items()}
    cisTransF = pd.DataFrame(cisTrans).reset_index()
    cisTransMolt = pd.melt(cisTransF, id_vars="index")

    # plot cisTrans

    f, ax = plt.subplots()
    sbn.barplot(x="variable", y="value", hue="index", data=cisTransMolt)
    ax.set_xlabel("Barcodes")
    ax.set_ylabel("Fraction Reads")
    ax.set_title("Cis Trans distribution")
    f.set_size_inches(10, 10)
    locs, labels = plt.xticks()
    plt.xticks(locs, labels, rotation=90)
    sbn.despine()
    # Shrink current axis by 20 %
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # Put a legend to the right of the current axis
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    f.savefig(os.path.join(resultsdir, "CisTrans_Dist.png"), bbox_inches="tight")

    # get cisDistance

    cisDist = {barcode: getCisDist(frame) for barcode, frame in stats.items()}
    cisDistF = pd.DataFrame(cisDist).reset_index()
    cisDistMolt = pd.melt(cisDistF, id_vars="index")

    # plot cisDistance

    f, ax = plt.subplots()
    sbn.barplot(
        x="variable",
        y="value",
        hue="index",
        data=cisDistMolt,
        hue_order=[
            "cis_1kb+",
            "cis_2kb+",
            "cis_4kb+",
            "cis_10kb+",
            "cis_20kb+",
            "cis_40kb+",
        ],
    )
    ax.set_xlabel("Barcodes")
    ax.set_ylabel("Fraction Reads")
    ax.set_title("Cis Distance distribution")
    f.set_size_inches(10, 10)
    locs, labels = plt.xticks()
    plt.xticks(locs, labels, rotation=90)
    sbn.despine()
    # Shrink current axis by 20 %
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # Put a legend to the right of the current axis
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    f.savefig(os.path.join(resultsdir, "CisTrans_Distance.png"), bbox_inches="tight")


if __name__ == "__main__":
    qc1()
