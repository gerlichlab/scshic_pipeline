# scshic_pipeline
## scsHi-C Preprocessing Nextflow Pipeline
IMBA - Gerlich Groupe <br>
christoph.langer@imba.oeaw.ac.at <br>
C.C.H. Langer, M. Mitter, A. Goloborodko

## A modular Hi-C mapping pipeline for reproducible data analysis.
[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.01.0-brightgreen.svg)]((http://nextflow.io))
![build](https://github.com/gerlichlab/scshic_pipeline/workflows/CI/badge.svg)
The `scshic` pipeline aims to provide the following preprocessing functionality:

- Align the sequences of Hi-C molecules to the reference genome
- Parse .sam alignment and create files with Hi-C pairs
- Filter PCR duplicates
- Annotate S4T mutation
- Filter cis and trans contacts
- Aggregate pairs into binned matrices of Hi-C interactions

## Installation

### Requirements

- java 8
- [nextflow](https://www.nextflow.io/)
- docker/singularity (should be able to run w/o root privileges, [tutorial](https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-16-04))
- an [optional container](https://hub.docker.com/repository/docker/gerlichlab/bcl2fastq) if you start with raw bcl files instead of demultiplexed fastq files
- the [main container](https://hub.docker.com/repository/docker/gerlichlab/scshic_docker) including all of the software needed to run the rest of the pipeline  
  
Here an example of how to prime our cluster for the pipeline:
```
module load singularity/3.7.3 # or whatever the newest version currently is (find via the command: module spider singularity)
module load nextflow/20.01.0 # or any newer version

singularity pull docker://docker.io/gerlichlab/bcl2fastq:latest
singularity pull docker://docker.io/gerlichlab/scshic_docker:release-1.5
``` 

### Pulling the pipeline
Now pull the pipeline directly with nextflow:

```
nextflow clone gerlichlab/scshic_pipeline ./
```

or if you have git installed, simply type:

```
git clone git@github.com:gerlichlab/scshic_pipeline ./
```

This will download the scsHi-C pipeline and the configuration files into the current working directory.

### Running the pipeline
**Attention: Running the pipeline is only recommended on a scientific cluster. For a NovaSeq flowcell, the working directory of nextflow will be >15TB large.**

Here is the example set up for a [SLURM](https://slurm.schedmd.com/documentation.html) batch system.

To run on the pipeline on a cluster, you need to modify not only `run_cbe.sh` but also the `/conf/cbe.conf` to adapt it to your cluster.
We have provided our run script and config file for our SLURM cluster as an example.

To launch the pipeline, you would run: 

```
bash run_cbe.sh
```
### Sample sheet for demultiplexing
For the demultiplexing with [bcl2fastq2](https://emea.support.illumina.com/downloads/bcl2fastq-conversion-software-v2-20.html) you need to provide it with a [sample sheet](https://emea.support.illumina.com/downloads/sample-sheet-v2-template.html)/
A minimal example is provided with `.\examples\Samplesheet.csv`; in it, the samples are named (Sample_ID) like the barcodes, and the sheet is for a NavaSeq flowcell with two lanes. 

**Attention: Only use alphanumeric characters for the Sample_ID.**

If you misformat the file, you will only get two files starting Undetermined*, and the pipeline will stop.

### Skipping demultiplexing:

If you provide the pipeline with a folder of demultiplexed .fastq files,
add the following parameter to the nextflow command in `run_cbe.sh`

**Attention: Make sure your filenames of your .fastq files are alphanumeric only**

```
--skipDemultiplexing "true"
```




## Detailed steps of the Pipeline:

1. Demultiplexing with [bcl2fastq2](https://emea.support.illumina.com/downloads/bcl2fastq-conversion-software-v2-20.html) 
If the Sequencing facility does
not provide you with demultiplexed.fastq this step is necessary to distinguish between the individual samples (often
different experimental conditions) and pooled samples.
2. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (on each sample) for raw sequencing reads quality control: FastQC is a tool that provides series of quality control checks on
raw sequence data. It provides a report including visualizations to give a
quick impression if the data has any heavier problems one should be aware of or consider before doing any further analysis.
3. Merge the fq-files of two lanes: Since a full Illumina NovaSeq run
includes two flowcells, a simple operation that concatenates two gzipped
fastq files performed with the unix ’cat’ command.
4. Alignment using [bwa-mem](https://doi.org/10.1109/IPDPS.2019.00041): The reads were aligned to a hg19 human genome with a correction for aberrant genome situations (SNPs and chromosomal translocations) of our HeLa cervical cancer cell lines used in our experiments. 
5. Parse: find ligation junctions in .sam, make .pairsam
6. Sort: Sorting the .pairsam files
7. Dedup: Find and remove PCR/optical duplicates. We used [pairtools](https://github.com/open2c/pairtools)
to detect ligation junctions of HiC pairs in aligned paired-end sequences
of Hi-C DNA molecules (done by parse). These .pairs files were sorted for
downstream analyses and PCR/optical duplicates were detected, tagged
and removed.
8. Annotate S4T mutations: We labeled a subset of DNA using synthetic
nucleotides, in our case, the newly synthesized sister chromatid and developed a novel method coined sister-chromatid sensitive Hi-C (scsHi-C) that allows the elucidation of sister chromatid structure at unprecedented resolution. During this step we were looking for the characteristic T-to-C and A-to-G mutations and decided afterward, according to heuristic, if the mutations are due to the labeling.
9. Filter cis and trans contacts: Filter for contacts within or between
sister chromatids.
10. Merge ref and comp .pairs file: Merge reference and complementary
strand reads using pairtools.
11. Generate cools: Create a cooler (i.e. .cool files) from genomic pairs
and bins, which is an efficient storage format for high-resolution genomic
interaction matrices using [cooler](https://github.com/open2c/cooler).
12. Zoomify and balance: generates a multi-resolution cooler file by coarsening and out-of-core matrix balancing. It is the final output and is loaded by the [HiGlass](https://docs.higlass.io/index.html) visualization system.
13. scsHi-C specific QC:
    - Fraction of labelled contacts for all samples.
    - Fraction of unique reads.
    - Trans-chromosomal/Cis-chromosomal contacts. 
    - Size-distribution of contacts.

### Citing
This tool was developed for the following paper:
[M. Mitter et al.](https://doi.org/10.1038/s41586-020-2744-4)

The protocol paper: M. Mitter, Z. Takacs. et al. (submitted)

### Versions
#### Versions 1.1.0
Is used for the papers.
#### Version 1.2.0
Is built to automatically start from a [FileMaker17](https://fmhelp.filemaker.com/help/17/fmp/en/) database.
The location of all output files are entered to the database after completion.


### License

MIT License
Copyright (c) 2019 Christoph C. H. Langer

### Components
The scsHi-C pipeline uses the following software components and tools:
- [bcl2fastq2](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html) (2.20.0) 
- [bwa mem](http://bio-bwa.sourceforge.net/) (0.7.17) as in other [Hi-C processing pipelines](https://data.4dnucleome.org/help/analysis-and-visualization/hi_c-processing-pipeline) we are using the -SP5M flag
    - -SP option is used to ensure the results are equivalent to that obtained by running bwa mem on each mate separately,
while retaining the right formatting for paired-end reads. This option skips a step in bwa mem that forces 
alignment of a poorly aligned read given an alignment of its mate with the assumption 
that the two mates are part of a single genomic segment.
    - -5 option is used to report the 5' portion of chimeric alignments as the primary alignment. 
In Hi-C experiments, when a mate has chimeric alignments, typically, the 5' portion is the position of interest,
while the 3' portion represents the same fragment as the mate. 
For chimeric alignments, bwa mem reports two alignments:
one of them is annotated as primary and soft-clipped, retaining the full-length of the original sequence. 
The other end is annotated as hard-clipped and marked as either 'supplementary' or 'secondary'. The -5 option forces the 5'end to be always annotated as primary.
    - -M option is used to annotate the secondary/supplementary clipped reads as secondary rather than supplementary, 
for compatibility with some public software tools such as picard MarkDuplicates.
- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (0.11.8)
- [pairtools](https://github.com/mirnylab/pairtools) (0.3.0)
- a custom s4T detection python/cython script included in this repo
- [cooler](https://github.com/mirnylab/cooler) (0.8.6)

### Containers
All components have been packaged into containers:
- bcl2fastq container can be found on [dockerhub](https://hub.docker.com/r/gerlichlab/bcl2fastq) and [git](https://github.com/cchlanger/bcl2fastq_docker)
- scshic container can be found on [dockerhub](https://hub.docker.com/r/gerlichlab/scshic_docker) and [git](https://github.com/gerlichlab/scshic_docker)

### Related Projects
- [Upstream analysis](https://github.com/gerlichlab/scshic_analysis) with the Jupyter notebooks that create all Figures of the scsHi-C paper.
- Use Michael's upstream [HiCTools](https://github.com/gerlichlab/ngs); they facilitate analysis of HiC data based on the cooler and cooltools interfaces.
- A python [wrapper](https://github.com/cchlanger/cooler_ontad) for [OnTAD](https://github.com/anlin00007/OnTAD) to use with resulting mcoolers. 
- Downstream analysis with [cooltools](https://github.com/mirnylab/cooltools)!
- Visualize your cooler data with [HiGlass](http://higlass.io)!
