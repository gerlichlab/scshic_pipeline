#!/usr/bin/env bash

DATE=$(date +"%Y-%m-%d-%H:%M:%S")

echo "========================================================
Log of the pipeline run is written to 'scshic_$DATE.log' and can be observed in this terminal.
It will contain a message when finished. This terminal can be closed.
========================================================" > scshic_$DATE.log


########################
#CHANGE FOLLOWING 6 PARAMETERS
########################
#Output location:
OUT=$(pwd)"/nf-out"
#Reference Genome Location, ours: hg19 fixed for Hela SNPs:
REFDIR="/groups/gerlich/members/MichaelMitter/Reference_genomes/Fasta/hg19_SNPs/"

#Iseq example
#Only use a number if you have multiple flowcells in the experiment, extend it like this: 046551, 046552 ...
EXP_ID="04655"
MACHINE="Iseq"
SAMPLE="/groups/gerlich/labinfo/scratch/samplesheet_4655_2.csv"
#Provided as BCL folder:
INPUT="/groups/gerlich/labinfo/scratch/20190730_FS10000507_18_BPC29604-1714/"

#NovaSeq example
#EXP_ID="04655"
#MACHINE="Novaseq"
#SAMPLE="/groups/gerlich/labinfo/scratch/samplesheet_4655_2.csv"
#INPUT="/groups/gerlich/labinfo/scratch/190802_A00700_0046_AHF2CFDRXX/"

########################
#READY TO USE AT THE IMBA - need to be modified at your institute
########################

#Needs to be done for consistency because nextflow would eat the first zero anyway
EXP_ID_CLEAN=$((10#$EXP_ID))
#Workdir on fast scratch - point this to your fastest filesystem
WORKDIR=/scratch-cbe/users/${USER}/nf-workdir/${EXP_ID_CLEAN}-${MACHINE}/
#Uncomment  the following line to use the current working directory for the temporary files
#WORKDIR=./nf-workdir/${EXP_ID_CLEAN}-${MACHINE}/
BASEDIR=$(pwd)
#Sends Email to default institute address
NOTIFICATION_EMAIL=${USER}@imba.oeaw.ac.at

mkdir -p $WORKDIR
cd $WORKDIR

#Used for the parameter --doubleflowcell, since NovaSeq machines have two lanes in a flowcell. 
if [ $MACHINE == "Iseq" ]
then
  FLOWCELL="false"
fi

if [ $MACHINE == "Novaseq" ]
then
  FLOWCELL="true"
fi

#Otherwise, it saves the singularity image into the current working directory
export NXF_SINGULARITY_CACHEDIR=$HOME/.singularity

#Redirects nohub into a log file, redirect the stderr to the same place we are redirecting the stdout and then starts tail to keep displaying the changing file.
nohup nextflow ${BASEDIR}/main.nf -profile cluster \
    --inputfolder $INPUT \
    --refDir $REFDIR \
    --sampleheet $SAMPLE \
    --machinetype $MACHINE \
    --experimentID $EXP_ID_CLEAN \
    --outdir $OUT \
    --stopBeforeCoolers true \
    --doubleflowcell $FLOWCELL \
    -with-timeline ${OUT}/${EXP_ID_CLEAN}-${MACHINE}/nxf_log/timeline.html \
    -with-report ${OUT}/${EXP_ID_CLEAN}-${MACHINE}/nxf_log/report.html \
    -with-dag ${OUT}/${EXP_ID_CLEAN}-${MACHINE}/nxf_log/nf_DAG.svg \
    -with-trace ${OUT}/${EXP_ID_CLEAN}-${MACHINE}/nxf_log/nf_trace.txt \
    -N $NOTIFICATION_EMAIL \
    -resume >> ${BASEDIR}/scshic_$DATE.log 2>&1 &


tail -n 1000 -f ${BASEDIR}/scshic_$DATE.log
