# mmhic-nf
## mmHiC Preprocessing Nextflow Pipline
IMBA - Gerlich Groupe <br>
christoph.langer@imba.oeaw.ac.at <br>
C.C.H. Langer, M. Mitter, A. Goloborodko

## A modular Hi-C mapping pipeline for reproducible data analysis.

The `scshic` pipeline aims to provide the following preprocessing functionality:

- Align the sequences of Hi-C molecules to the reference genome
- Parse .sam alignment and form files with Hi-C pairs
- Filter PCR duplicates
- Annotate S4T mutation
- Filter cis and trans contacts
- Aggregate pairs into binned matrices of Hi-C interactions

### Installation

Requirements:

- java 8
- [nextflow](https://www.nextflow.io/)
- docker/singularity (should be able to run w/o root privileges, 
[tutorial](https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-16-04))

As long as the repository is privat ad this to your scm file (or create it) @ ~/.nextflow/scm
```
providers {
    github {
        user = 'me'
        password = 'my-secret-password'
    }
}
```
TODO: make static location of lib possible
To setup a new project, execute the following line in the project folder:

```
nextflow clone gerlichlab/mmhic-nf ./
```

or if you have git installed simply type:

```
git clone git@github.com:cchlanger/mmhic_nf.git ./
```

This will download the mmhic pipeline and the configuration files.

Then:

- modify set/change the parameters in run_cbe.sh
    An example for the 5 mandatory parameters for the Novaseq and Iseq parameters are included
- the default hardware configuration is cbe
- in the future (TODO) use provided hardware configurations using `cbe-novaseq` or `cbe-iseq` profiles (or `cbe-slow` in case of fast memory problems on CBE), or provide your own using `custom` profile

Launch mmhic depending on your usage scenario:

```
bash run_cbe.sh
```

## Steps of the Pipline:
1. Demultiplexing with bcl2fastq2 [17] 
If the Sequencing facility does
not provide you with demultiplexed.fastq or unalligned .sam/.bam files
this step is necessary to distinguish between the individual samples (often
different experimental conditions) and pooled samples.
2. FastQC (on each sample) for raw sequencing reads quality con-
trol: FastQC is a tool that provides series of quality control checks on
23raw sequence data. It provides a report including visualizations, to give a
quick impression if the data has any heavier problems one should be aware
of or consider,before doing any further analysis [39].
3. Merge the fq-files of two lanes: Since a full Illumina Novaseq run
includes two flowcells, a simple operation that concatenates two gziped
fastq files performed with the unix ’cat’ command.
4. Alignment using bwa mem[27][25]: The reads where aligned to the
hg19 human genome with a correction for aberrant genome situation (he
SNPs and chromosomal translocations) of our HeLa cervical cancer cell
line used in our experiments. Unfortunately the bwa mem algorithm has
not yet been published.
5. Parse: find ligation junctions in .sam, make .pairsam
6. Sort: Sorting the .pairsam files
7. Dedup: Find and remove PCR/optical duplicates. We used pairtools [4]
to detect ligation junctions of HiC pairs in aligned paired-end sequences
of Hi-C DNA molecules (done by parse), These .pairs files were sorted for
downstream analyses and PCR/optical duplicates were detected, tagged
and removed.
8. Annotate S4T mutations: We labeled a subset of DNA using synthetic
nucleotides, in our case the newly synthesized sister chromatid and de-
veloped a novel method coined sister-chromatid sensitive Hi-C (scsHi-C)
that allows the elucidation of sister chromatid structure at unprecedented
resolution. During this step we were looking for the characteristic T-to-C
and A-to-G mutations and decided afterwards, according to heuristic, if
the mutations is due to the labeling.
9. Filter cis and trans contacts: Filter for contacts within or between
sister chromatids.
10. Merge ref and comp .pairs file: Merge reference and complementary
strand reads. Using pairtools.
11. Generate cools: Create a cooler (i.e. .cool files) from genomic pairs
and bins, which is an efficient storage format for high resolution genomic
interaction matrices using cooler [1].
12. Zoomify and balance: generates a multi-resolution cooler file by coarsen-
ing and out-of-core matrix balancing. It is the final output and is loaded
by the HiGlass visualization system [21]

### Test example

In a new project folder, execute:

```bash
$ nextflow clone gerlichlab/mmhic-nf  ./
$ bash run_cbe.sh 
```
### Citing
TODO: Add Michael's paper
Mitter, M., et al. (2019). mmHiC: ... doi: [xxx](xxx).

```bibtex
@article{mmhic2019,
    author = {Mitter, Michael et al.},
    title = "{mmhic: ...}",
    journal = {},
    year = {2019},
    month = {},
    doi = {},
    url = {},
}
```

### License

MIT License
Copyright (c) 2019 Christoph C. H. Langer

### Components
TODO add more Details! <br>
mmhic-nf uses the following software components and tools:
- [bcl2fastq2](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html) 2.20.0 
- [bwa mem](http://bio-bwa.sourceforge.net/) 0.7.17 as in other [Hi-C processing pipelines](https://data.4dnucleome.org/help/analysis-and-visualization/hi_c-processing-pipelin) we are using the -SP5M flag
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
- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 0.11.8
- [pairtools](https://github.com/mirnylab/pairtools) 0.3.0
- a custom s4T detection python/cython scripts included in this repo
- [cooler](https://github.com/mirnylab/cooler) 0.8.6

### Containers
All components have been packaged into containers:
- bcl2fastq container can be found on [dockerhub](https://hub.docker.com/r/gerlichlab/bcl2fastq) and [git](https://github.com/cchlanger/bcl2fastq_container)
- mmhic-nf container can be found on [dockerhub](https://hub.docker.com/r/gerlichlab/mmhic) and [git](https://github.com/cchlanger/mmHiC)

### Related Projects
- Use Michael's [HiCTools](https://github.com/gerlichlab/NGS)!
- Downstream analysis with [cooltools](https://github.com/mirnylab/cooltools)!
- Visualize your cooler data with [HiGlass](http://higlass.io)!
- A modular Hi-C mapping pipeline called [Distiller](https://github.com/mirnylab/distiller-nf)
