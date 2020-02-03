#!/usr/bin/env bash

DATE=$(date +"%Y%m%d%H%M%S")

#If you can not use containers below you will find the minimum installation and conda environment
#module load singularity/3.1.0
#module load anaconda3/2019.03
#module load gcc/7.3.0-2.30
#module load bwa/0.7.17-foss-2018b
#module load fastqc/0.11.8-java-1.8

#This is my very special environment with everthing installed
#conda create --name snowflake
#conda activate snowflake
#conda install -c conda-forge -c bioconda nextflow pbgzip pairtools cooler pysam numpy pandas click cython conda graphviz 
#conda activate snowflake

echo "========================================================
Log of the pipeline run is written to 'mmhic_$DATE.log' and can be observed in this terminal.
It will contain a message when finished. This terminal can be closed.
========================================================" > mmhic_$DATE.log


########################
#CHANGE FOLLOWING 5 PARAMETERS
########################
#Small example
OUT="./nf-out"
#Only use a number if you have multiple flowcells in the experiment, extend it like this: 046551, 046552, ... TODO: Fix this bahaviour.
EXP_ID="04655"
MACHINE="Iseq"
SAMPLE="/groups/gerlich/labinfo/scratch/samplesheet_4655.csv"
INPUT="/groups/gerlich/labinfo/scratch/20190730_FS10000507_18_BPC29604-1714/"

#Big example
#OUT="./nf-out"
#EXP_ID="04655"
#MACHINE="Novaseq"
#SAMPLE="/groups/gerlich/labinfo/scratch/samplesheet_4655_2.csv"
#INPUT="/groups/gerlich/labinfo/scratch/190802_A00700_0046_AHF2CFDRXX/"

########################
#ARE READY TO USE AT THE IMBA
#Needs to be done for consistancy because nextflow would eat first zero anyway
EXP_ID_CLEAN=$((10#$EXP_ID))
#Workdir on fast scratch
WORKDIR=/scratch-cbe/users/${USER}/nf-workdir/${EXP_ID_CLEAN}-${MACHINE}/
#Slow backup dir
#WORKDIR=./nf-workdir/${EXP_ID_CLEAN}-${MACHINE}/
BASEDIR=$(pwd)
#Sends Email to default institute adresse
NOTIFICATION_EMAIL=${USER}@imba.oeaw.ac.at

mkdir -p $WORKDIR
cd $WORKDIR

if [ $MACHINE == "Iseq" ]
then
  FLOWCELL="false"
fi

if [ $MACHINE == "Novaseq" ]
then
  FLOWCELL="true"
fi

#Otherwise it saves the singularity image into the workdir
export NXF_SINGULARITY_CACHEDIR=$HOME/.singularity

#Redirects nohub into a log file, redirect the stderr to the same place we are redirecting the stdout and then starts tail to keep displaying the changing file.
nohup nextflow ${BASEDIR}/main.nf -profile cbe \
    --inputfolder $INPUT \
    --sampleheet $SAMPLE \
    --machinetype $MACHINE \
    --experimentID $EXP_ID_CLEAN \
    --outdir $OUT \
    --doubleflowcell $FLOWCELL \
    -with-timeline ${OUT}/${EXP_ID_CLEAN}-${MACHINE}/nxf_log/timeline.html \
    -with-report ${OUT}/${EXP_ID_CLEAN}-${MACHINE}/nxf_log/report.html \
    -with-dag ${OUT}/${EXP_ID_CLEAN}-${MACHINE}/nxf_log/nf_DAG.svg \
    -with-trace ${OUT}/${EXP_ID_CLEAN}-${MACHINE}/nxf_log/nf_trace.txt \
    -N $NOTIFICATION_EMAIL \
    -resume >> ${BASEDIR}/mmhic_$DATE.log 2>&1 &
#TODO Test new logging


tail -n 1000 -f ${BASEDIR}/mmhic_$DATE.log
