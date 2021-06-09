with open("run_generated.sh", "w+") as runfile:
    runfile.write(
        """
#!/usr/bin/env bash

DATE=$(date +"%Y-%m-%d-%H:%M:%S")

echo "========================================================
Log of the pipeline run is written to 'scshic_$DATE.log' and can be observed in this terminal.
It will contain a message when finished. This terminal can be closed.
========================================================" > scshic_$DATE.log

"""
    )

    # TODO: Change to final location
    OUT = '$(pwd)"/nf-out"'
    runfile.write(f"OUT={OUT}\n")
    # TODO: Change to final location
    REFDIR = "/groups/gerlich/members/MichaelMitter/Reference_genomes/Fasta/hg19_SNPs/"
    runfile.write(f"REFDIR={REFDIR}\n")
    # TODO: Import from Database
    EXP_ID = "04655"
    runfile.write(f"EXP_ID={EXP_ID}\n")
    # TODO: Import from Database
    MACHINE = "Iseq"
    runfile.write(f"MACHINE={MACHINE}\n")
    # TODO: Import from Database
    SAMPLE = "/groups/gerlich/labinfo/scratch/samplesheet_4655_2.csv"
    runfile.write(f"SAMPLE={SAMPLE}\n")
    # TODO: Import from Database
    INPUT = "/groups/gerlich/labinfo/scratch/20190730_FS10000507_18_BPC29604-1714/"
    runfile.write(f"INPUT={INPUT}\n")

    runfile.write(
        """
#Needs to be done for consistency because nextflow would eat the first zero anyway
EXP_ID_CLEAN=$((10#$EXP_ID))
#Workdir on fast scratch - point this to your fastest filesystem
WORKDIR=/scratch-cbe/users/${USER}/nf-workdir/${EXP_ID_CLEAN}-${MACHINE}/
#Uncomment  the following line to use the current working directory for the temporary files
#WORKDIR=./nf-workdir/${EXP_ID_CLEAN}-${MACHINE}/
BASEDIR=$(pwd)
"""
    )
    # Sends Email to default institute address
    # TODO get from the database
    NOTIFICATION_EMAIL = "${USER}@imba.oeaw.ac.at"
    runfile.write(f"NOTIFICATION_EMAIL={NOTIFICATION_EMAIL}\n")

    runfile.write(
        """
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
    --doubleflowcell $FLOWCELL \
    -with-timeline ${OUT}/${EXP_ID_CLEAN}-${MACHINE}/nxf_log/timeline.html \
    -with-report ${OUT}/${EXP_ID_CLEAN}-${MACHINE}/nxf_log/report.html \
    -with-dag ${OUT}/${EXP_ID_CLEAN}-${MACHINE}/nxf_log/nf_DAG.svg \
    -with-trace ${OUT}/${EXP_ID_CLEAN}-${MACHINE}/nxf_log/nf_trace.txt \
    -N $NOTIFICATION_EMAIL \
    -resume >> ${BASEDIR}/scshic_$DATE.log 2>&1 &


tail -n 1000 -f ${BASEDIR}/scshic_$DATE.log
"""
    )
