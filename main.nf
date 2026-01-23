#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================
            IMBA - Gerlich Groupe - scsHi-C PREPROCESSING PIPELINE
========================================================================

 christoph.langer@imba.oeaw.ac.at
 C.C.H. Langer, M. Mitter, A. Goloborodko

========================================================================================
----------------------------------------------------------------------------------------
Pipeline overview:
 - 1:   Demultiplexing with bcl2fastq
 - 2:   FastQC (on each sample) for raw sequencing reads quality control
 - 3:   Merge the fq-files of two lanes
 - 4:   Alignment using bwa mem
 - 5.1: Parse: find ligation junctions in .sam, make .pairsam
 - 5.2: Sort the .pairsam files
 - 5.3: Dedup: Find and remove PCR/optical duplicates.
 - 6:   Annotate S4T mutations
 - 7.1: Filter cis and trans contacts
 - 7.2: Merge ref and comp .pairs file
 - 8:   Generate cools
 - 9:   Zoomify and balance: generate a multi-resolution cooler file by coarsening and out-of-core matrix balancing
 - 10:  Run QC scripts and copy all output to final destination
 ----------------------------------------------------------------------------------------
*/

//TODO Naming into Camel Case

//Mandatory parameters need to be set in run script
params.inputfolder 
params.sampleheet 
params.machinetype 
params.experimentID 
params.doubleflowcell
//How the hg19 reference genome fixed for Hela SNPs is explained in the paper.
params.refDir
outputFolder = "$params.experimentID-$params.machinetype"

//Default parameters can be overwritten when calling the pipeline (i.e in the run_xxx.sh scripts):

params.skipFastQC = false
params.skipDemultiplexing = false
//TODO should this be a parameter? Or fetched automatically?
params.taskcpu = 38
params.chrSizes="${baseDir}/bin/hg19.chrom.sizes"
params.outdir="${baseDir}/nf-out"
// Stop before cooler generation (skip steps 8-10)
params.stopBeforeCoolers = false

//Parameters for s4t mutations
params.min_map_q = 30
params.min_right_muts = 2
params.max_wrong_muts = 0

//Performance Parameters
bwa_threads=38
memory_merge = "20G"
nproc_merge = 19
nproc_pairsam = 19

//Step 5.2 parameters
nproc_sort = 19
nproc_out_sort = 19
memory_sort = "20G"


def helpMessage() {
    log.info"""

========================================================================
            IMBA - Gerlich Groupe - scsHi-C PREPROCESSING PIPELINE
========================================================================

 christoph.langer@imba.oeaw.ac.at
 C.C.H. Langer, M. Mitter, A. Goloborodko

========================================================================

    Usage:
    The typical command for running the pipeline is using a run script,
    modified for your experiment, examples can be found on github:
    
    bash run_cbe.sh
    
    For more information and detailed documentation see: 
    www.github.com/gerlichlab/scshic_pipeline
========================================================================
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = 1.1

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

/*
 * STEP 1 - Demultiplexing with bcl2fastq
 */

process bcl2fastq{
    publishDir "$outputFolder/demultiplexed", mode: "symlink"
    output:
        file("*.fastq.gz") into CH_demult_raw
    when:
        !params.skipDemultiplexing
    script:
        """
        bcl2fastq -R $params.inputfolder --sample-sheet $params.sampleheet -o .
        """
}

/*
 * STEP 2 - FastQC (on each sample) for raw sequencing reads quality control 
 */

CH_demult_raw
    .flatMap() //Attention: Do not forget forget to use .flatMap() after this to detect multiple outputs, since bcl2fastq ran only once and put all files into one item of the channel. 
    .into { CH_demult_for_fastqc; CH_demult_temp }


process fastqc{
    publishDir "$outputFolder/fastqc", mode: 'symlink',
        saveAs: {filename -> filename.endsWith(".zip") ? "zips/$filename" : "$filename"}
    input:
        file(fastq) from CH_demult_for_fastqc
    output:
        file("*.{zip,html}") into CH_fastqc_reports
    when:
        !params.skipFastQC
    script:
        """
        fastqc -t $params.taskcpu -q $fastq
        """
}

CH_demult_temp
            .map{ item ->
                if (! "${item}".contains("Undetermined_")){
                    return item
                }
            }
            .set{ CH_demult_clean }

CH_demult_clean
    .into {CH_demult_iseq; CH_demult_novaseq}
CH_demult_novaseq
    .into {CH_demult_novaseq_temp1; CH_demult_novaseq_temp2}

//TODO replace with choice operator

CH_demult_novaseq_temp1
    .map { item -> 
            if ( "${item}".contains("_R1_001.fastq.gz")){
                return item
            }
    }
    .map { file ->
        def key = file.name.toString().tokenize('_').get(0)
        return tuple(key, file)
     }
    .groupTuple(size: 2, sort: true)
    .set{CH_lane_pairs_R1}

CH_demult_novaseq_temp2
    .map { item -> 
            if ( "${item}".contains("_R2_001.fastq.gz")){
            //TODO: replace with filter and "*_L00{1,2}_R2*.fastq.gz"
                return item
            }
    }
        .map { file ->
        def key = file.name.toString().tokenize('_').get(0)
        return tuple(key, file)
     }
    .groupTuple(size: 2, sort: true)
    .set{CH_lane_pairs_R2}

/*
 * STEP 3 - Merge the fq-files of two lanes
 */

//Iseq has only one lane, simply renames file
process rename_lanes{
    publishDir "$outputFolder/merged_fastq", mode: 'symlink'
    input:
        file(lane_1) from CH_demult_iseq.filter(~/^.*_001\.fastq\.gz/)
    output:
        file("*R{1,2}.fastq.gz") into CH_iseq_for_align
    when:
        params.machinetype == 'Iseq'
        params.doubleflowcell == false
    script:
        barcode = lane_1.getName().tokenize('_').get(0)
        sample_ID = lane_1.getName().tokenize('_').get(1)
        read_direction = lane_1.getName().tokenize('_').get(3)
        file_name =  "${barcode}_${sample_ID}_${read_direction}"
        """
        mv $lane_1 ${file_name}.fastq.gz
        """
}

// In a full Novaseq run two lanes need to be merged
process merge_lanes_R1{
    publishDir "$outputFolder/merged_fastq", mode: 'symlink'
    input:
        set pair_id, file(pair) from CH_lane_pairs_R1
    output:
        file("*.fastq.gz") into CH_merged_R1
    when:
        params.machinetype == 'Novaseq'
        params.doubleflowcell == true
    script:
        barcode = pair[0].getName().tokenize('_').get(0)
        assert pair_id == barcode
        sample_ID = pair[0].getName().tokenize('_').get(1)
        read_direction = pair[0].getName().tokenize('_').get(3)
        file_name = "${barcode}_${sample_ID}_${read_direction}"
        """
        cat ${pair} > $file_name".fastq.gz"
        """ 
}

process merge_lanes_R2{
    publishDir "$outputFolder/merged_fastq", mode: 'symlink'
    input:
        set pair_id, file(pair) from CH_lane_pairs_R2
    output:
        file("*.fastq.gz") into CH_merged_R2
    when:
        params.machinetype == 'Novaseq'
        params.doubleflowcell == true
    script:
        barcode = pair[0].getName().tokenize('_').get(0)
        assert pair_id == barcode
        sample_ID = pair[0].getName().tokenize('_').get(1)
        read_direction = pair[0].getName().tokenize('_').get(3)
        file_name = "${barcode}_${sample_ID}_${read_direction}"
        """
        cat ${pair} > $file_name".fastq.gz"
        """ 
}


CH_iseq_for_align.map { file ->
        def key = file.name.toString().tokenize('_').get(0)
        return tuple(key, file)
     }
    .groupTuple(size: 2, sort: true)
    .set{ CH_iseq_tuple_for_align}

CH_merged_R1.mix(CH_merged_R2)
    .map { file ->
        def key = file.name.toString().tokenize('_').get(0)
        return tuple(key, file)
     }
    .groupTuple(size: 2, sort: true)
    .set{CH_novaseq_tuple_for_align}

//Since both are exclusive this just merges two branches back again
CH_iseq_tuple_for_align.mix(CH_novaseq_tuple_for_align)
    .set{CH_tuple_for_align}

//Caution: groupTuple on ids creates: [id [file1, file2]], joining channels on id will produce: [id file1, file2]

/*
 * STEP 4 - Alignement using bwa mem
 */

// TODO sanity check fastq file names - parse till end

process align_w_bwa{
    publishDir "$outputFolder/bwa_results", mode: 'symlink'
    input:
        set pair_id, file(tuple) from CH_tuple_for_align
    output:
        file "*.sam" into CH_bwa_results
    script:
        barcode = tuple[0].getName().tokenize('_').get(0)
        assert pair_id == tuple[0].getName().tokenize('_').get(0)
        assert pair_id == tuple[1].getName().tokenize('_').get(0)
        sample_ID = tuple[0].getName().tokenize('_').get(1)
        file_name = "${barcode}_${sample_ID}"
        """
        bwa mem -SP5M -t $bwa_threads ${params.refDir}hg19_SNPs.fa ${tuple} > ${file_name}_bwa.sam
        """
}


/*
 * STEP 5 - Parse: find ligation junctions in .sam, make .pairsam
 */

process parse_pairs{
    publishDir "$outputFolder/parsed_pairsam", mode: 'symlink'
    input:
        file(bwa_sam) from CH_bwa_results
    output:
        file "*.pairsam" into CH_parsed_pairsam
        file "*.txt" into parsed_stats
    script:
        barcode = bwa_sam.getName().tokenize('_').get(0)
        sample_ID = bwa_sam.getName().tokenize('_').get(1)
        file_name = "${barcode}_${sample_ID}"
        """
        pairtools parse -c $params.chrSizes --add-columns mapq -o ${file_name}.pairsam --nproc-in $nproc_pairsam --output-stats ${file_name}_stats.txt $bwa_sam
        """
}


/*
 * STEP 5.2 - Sort the .pairsam files
 */

process sort_pairs{
    publishDir "$outputFolder/sorted_pairsam", mode: 'symlink'
    input:
        file(pairsam) from CH_parsed_pairsam
    output:
        file "*.pairsam.lz4" into CH_sorted_pairsam
    script:
        barcode = pairsam.getName().tokenize('_').get(0)
        sample_ID = pairsam.getName().tokenize('_').get(1).tokenize('.').get(0)
        file_name = "${barcode}_${sample_ID}"
        """
        pairtools sort --memory $memory_sort -o ${file_name}.sorted.pairsam.lz4 --nproc-out $nproc_out_sort  --nproc $nproc_sort $pairsam
        """
}


/*
 * STEP 5.3 - Dedup: Find and remove PCR/optical duplicates.
 */

process dedup_pairs{
    publishDir "$outputFolder/dedup_pairsam", mode: 'symlink',
        saveAs: {filename -> filename.endsWith(".stats") ? "stats/$filename" : "$filename"}
    input:
        file(pairsam) from CH_sorted_pairsam
    output:
        file "*.pairsam.gz" into CH_dedup_pairsam_main
        file "*.stats" into CH_dedup_stats
    script:
        barcode = pairsam.getName().tokenize('_').get(0)
        sample_ID = pairsam.getName().tokenize('_').get(1).tokenize('.').get(0)
        file_name = "${barcode}_${sample_ID}"
        """
        pairtools dedup -o ${file_name}.dedup.pairsam.gz --output-stats ${file_name}_dedup.stats $pairsam
        """
}

CH_dedup_pairsam_main.into { CH_dedup_pairsam; CH_dedup_for_cooler; CH_dedup_for_copy }
CH_dedup_for_copy.into { CH_dedup_for_copy_iseq; CH_dedup_for_copy_novaseq }

// Convert deduped pairsam (all contacts) to pairs for cooler generation and downstream copying
process pairsam_to_pairs_all{
    publishDir "$outputFolder/pairs_all", mode: 'symlink'
    input:
        file(pairsam) from CH_dedup_for_cooler
    output:
        file "*.pairs.gz" into CH_all_pairs_for_cooler
    script:
        barcode = pairsam.getName().tokenize('_').get(0)
        sample_ID = pairsam.getName().tokenize('_').get(1).tokenize('.').get(0)
        file_name = "${barcode}_${sample_ID}"
        """
        pairtools split --nproc ${nproc_pairsam} --output-pairs ${file_name}.all.pairs.gz --output-sam /dev/null ${pairsam}
        """
}

/*
 * STEP 6 - Annotate S4T mutations
 */

process detect_s4t{
    publishDir "$outputFolder/s4t_pairsam", mode: 'symlink'
    input:
        file(pairsam) from CH_dedup_pairsam
    output:
        file "*.pairsam.gz" into CH_s4T_pairsam
    script:
        barcode = pairsam.getName().tokenize('_').get(0)
        sample_ID = pairsam.getName().tokenize('_').get(1).tokenize('.').get(0)
        file_name = "${barcode}_${sample_ID}"
        //TODO make detect_s4t_mutations.py concurrent
        //Reminder: $baseDir is the location of main.nf
        """
        python ${baseDir}/bin/detect_s4t_mutations.py --chunksize 5000000 -o ${file_name}.dedup.s4t.pairsam.gz $pairsam
        """
}

/*
 * STEP 7.1 - Filter cis and trans contacts
 */

process filter_cis_trans{
    publishDir "$outputFolder/s4t_filtered_pairsam", mode: 'symlink'
    input:
        file(s4t_pairs) from CH_s4T_pairsam
    output:
        file "*right_muts*trans*pairs.gz" into CH_trans_pairs
        file "*right_muts*cis*pairs.gz" into CH_cis_pairs
    script:
        barcode = s4t_pairs.getName().tokenize('_').get(0)
        sample_ID = s4t_pairs.getName().tokenize('_').get(1).tokenize('.').get(0)
        sample_name = "${barcode}_${sample_ID}"
        """
        bash ${baseDir}/bin/filter_single_experiment.sh $sample_name . $s4t_pairs $params.min_map_q $params.min_right_muts $params.max_wrong_muts
        """
}

/*
 * STEP 7.2 - Merge ref and comp .pairs file
 */

CH_cis_pairs
    .flatMap()
    .map { file ->
        def key = file.name.toString().tokenize('_').get(0)
        return tuple(key, file)
     }
    .groupTuple(size: 2, sort: true)
    .set{ CH_cis_pairs_file }

process merge_cis_ref_comp{
    publishDir "$outputFolder/s4t_merged_pairsam", mode: 'symlink'
    input:
        set pair_id, file(set) from CH_cis_pairs_file
    output:
        file "*cis.pairs.gz" into CH_cis_merged_pairs
    script:
        barcode = set[0].getName().tokenize('_').get(0)
        sample_ID = set[0].getName().tokenize('_').get(1).tokenize('.').get(0)
        sample_name = "${barcode}_${sample_ID}"
        """
        pairtools merge --memory ${memory_merge} --nproc ${nproc_merge} -o ${sample_name}.cis.pairs.gz ${set}
        """
}

CH_trans_pairs
    .flatMap()
    .map { file ->
        def key = file.name.toString().tokenize('_').get(0)
        return tuple(key, file)
     }
    .groupTuple(size: 2, sort: true)
    .set{ CH_trans_pairs_file }


process merge_trans_ref_comp{
    publishDir "$outputFolder/s4t_merged_pairsam", mode: 'symlink'
    input:
        set pair_id2, file(set) from CH_trans_pairs_file
    output:
        file "*trans.pairs.gz" into CH_trans_merged_pairs
    script:
        barcode = set[0].getName().tokenize('_').get(0)
        sample_ID = set[0].getName().tokenize('_').get(1).tokenize('.').get(0)
        sample_name = "${barcode}_${sample_ID}"
        """
        pairtools merge --memory ${memory_merge} --nproc ${nproc_merge} -o ${sample_name}.trans.pairs.gz ${set}
        """
}

/*
 * STEP 8 - Generate cools
 */

process generate_cools{
    publishDir "$outputFolder/cooler", mode: 'symlink'
    input:
        file (all_pairs) from CH_all_pairs_for_cooler
        file (trans_pairs) from CH_trans_merged_pairs
        file (cis_pairs) from CH_cis_merged_pairs
    output:
        file "*cool" into CH_cools
    when:
        !params.stopBeforeCoolers
    script:
        trans_barcode = trans_pairs.getName().tokenize('_').get(0)
        trans_sample_ID = trans_pairs.getName().tokenize('_').get(1).tokenize('.').get(0)
        trans_sample_name = "${trans_barcode}_${trans_sample_ID}"
        
        all_barcode = all_pairs.getName().tokenize('_').get(0)
        all_sample_ID = all_pairs.getName().tokenize('_').get(1).tokenize('.').get(0)
        all_sample_name = "${all_barcode}_${all_sample_ID}"

        cis_barcode = cis_pairs.getName().tokenize('_').get(0)
        cis_sample_ID = cis_pairs.getName().tokenize('_').get(1).tokenize('.').get(0)
        cis_sample_name = "${cis_barcode}_${cis_sample_ID}"
        """
        # All contacts
        cooler cload pairs ${params.chrSizes}:1000 $all_pairs ${all_sample_name}.all.1000.cool --chrom1 2 --chrom2 4 --pos1 3 --pos2 5
        # Cis contacts
        cooler cload pairs ${params.chrSizes}:1000 $cis_pairs ${cis_sample_name}.cis.1000.cool --chrom1 2 --chrom2 4 --pos1 3 --pos2 5
        # Trans contacts
        cooler cload pairs ${params.chrSizes}:1000 $trans_pairs ${trans_sample_name}.trans.1000.cool --chrom1 2 --chrom2 4 --pos1 3 --pos2 5
        """
}

CH_cools.into { CH_cools_for_balancing; CH_cools_iseq }

/*
 * STEP 9 - Zoomify and balance: generate a multi-resolution cooler file by coarsening and out-of-core matrix balancing
 */

process zoomify_and_balance{
    publishDir "$outputFolder/balanced_cooler", mode: 'symlink'
    input:
        file (cooler_2) from CH_cools_for_balancing.flatMap()
    output:
        file "*mcool" into CH_mcools_novaseq
    when:
        params.machinetype=='Novaseq' && !params.stopBeforeCoolers
    script:
        """
        cooler zoomify ${cooler_2} -n 4 -r 1000,2000,4000,5000,6000,8000,10000,20000,50000,100000,200000,400000,500000,1000000,5000000 --balance --balance-args '--ignore-diags 0 --mad-max 5 --max-iters 500 --convergence-policy store_nan'
        """
}

/*
 * STEP 10 - Run QC scripts and copy all output to final destination
*/

process copy_to_output_iseq{
    //fake dependency for synchronization (barrier function)
    input:
        val (out) from CH_cools_iseq.collect()
    when:
        params.machinetype=='Iseq' && !params.stopBeforeCoolers
    script:
        //Creates the Folder is it does not exist yet
        //Copy to output and dereference symlinks
        """
        mkdir -p $params.outdir/$outputFolder/qc_and_stats/qc
        cp -rL ../../../$outputFolder/merged_fastq $params.outdir/$outputFolder/merged_fastq
        cp -rL ../../../$outputFolder/cooler $params.outdir/$outputFolder/unbalanced_cooler
        cp -rL ../../../$outputFolder/dedup_pairsam/stats $params.outdir/$outputFolder/qc_and_stats/stats
        cp -rL ../../../$outputFolder/fastqc $params.outdir/$outputFolder/qc_and_stats/fastqc
        cp -rL ../../../$outputFolder/pairs_all $params.outdir/$outputFolder/all_pairs
        cp -rL ../../../$outputFolder/s4t_merged_pairsam $params.outdir/$outputFolder/cis_trans_pairs
        python ${baseDir}/bin/check_duplications_trans_cis_dist.py --inputdir $params.outdir/$outputFolder/qc_and_stats/stats --resultsdir $params.outdir/$outputFolder/qc_and_stats/qc
        python ${baseDir}/bin/count_double_labeled.py --inputdir ../../../$outputFolder --resultsdir $params.outdir/$outputFolder/qc_and_stats/qc
        """
}

process copy_to_output_novaseq{
    //fake dependency for synchronization (barrier function)
    input:
        val (out) from CH_mcools_novaseq.collect()
    when:
        params.machinetype=='Novaseq' && !params.stopBeforeCoolers
    script:
        //Creates the Folder is it does not exist yet
        //Copy to output and dereference symlinks
        """
        mkdir -p $params.outdir/$outputFolder/qc_and_stats/qc
        cp -rL ../../../$outputFolder/merged_fastq $params.outdir/$outputFolder/merged_fastq
        cp -rL ../../../$outputFolder/balanced_cooler $params.outdir/$outputFolder/mcooler
        rm $params.outdir/$outputFolder/mcooler/*cis.1000.mcool
        rm $params.outdir/$outputFolder/mcooler/*trans.1000.mcool
        cp -rL ../../../$outputFolder/dedup_pairsam/stats $params.outdir/$outputFolder/qc_and_stats/stats
        cp -rL ../../../$outputFolder/cooler $params.outdir/$outputFolder/unbalanced_cooler
        cp -rL ../../../$outputFolder/fastqc $params.outdir/$outputFolder/qc_and_stats/fastqc
        cp -rL ../../../$outputFolder/pairs_all $params.outdir/$outputFolder/all_pairs
        cp -rL ../../../$outputFolder/s4t_merged_pairsam $params.outdir/$outputFolder/cis_trans_pairs
        python ${baseDir}/bin/check_duplications_trans_cis_dist.py --inputdir $params.outdir/$outputFolder/qc_and_stats/stats --resultsdir $params.outdir/$outputFolder/qc_and_stats/qc
        python ${baseDir}/bin/count_double_labeled.py --inputdir ../../../$outputFolder --resultsdir $params.outdir/$outputFolder/qc_and_stats/qc
        """
}

// Copy outputs when stopping before coolers (no cooler artifacts available)
process copy_to_output_iseq_precooler{
    input:
        val (out) from CH_dedup_for_copy_iseq.collect()
    when:
        params.machinetype=='Iseq' && params.stopBeforeCoolers
    script:
        """
        mkdir -p $params.outdir/$outputFolder/qc_and_stats/qc
        cp -rL ../../../$outputFolder/merged_fastq $params.outdir/$outputFolder/merged_fastq
        cp -rL ../../../$outputFolder/dedup_pairsam/stats $params.outdir/$outputFolder/qc_and_stats/stats
        cp -rL ../../../$outputFolder/fastqc $params.outdir/$outputFolder/qc_and_stats/fastqc
        cp -rL ../../../$outputFolder/pairs_all $params.outdir/$outputFolder/all_pairs
        cp -rL ../../../$outputFolder/s4t_merged_pairsam $params.outdir/$outputFolder/cis_trans_pairs
        python ${baseDir}/bin/check_duplications_trans_cis_dist.py --inputdir $params.outdir/$outputFolder/qc_and_stats/stats --resultsdir $params.outdir/$outputFolder/qc_and_stats/qc
        python ${baseDir}/bin/count_double_labeled.py --inputdir ../../../$outputFolder --resultsdir $params.outdir/$outputFolder/qc_and_stats/qc
        """
}

process copy_to_output_novaseq_precooler{
    input:
        val (out) from CH_dedup_for_copy_novaseq.collect()
    when:
        params.machinetype=='Novaseq' && params.stopBeforeCoolers
    script:
        """
        mkdir -p $params.outdir/$outputFolder/qc_and_stats/qc
        cp -rL ../../../$outputFolder/merged_fastq $params.outdir/$outputFolder/merged_fastq
        cp -rL ../../../$outputFolder/dedup_pairsam/stats $params.outdir/$outputFolder/qc_and_stats/stats
        cp -rL ../../../$outputFolder/fastqc $params.outdir/$outputFolder/qc_and_stats/fastqc
        cp -rL ../../../$outputFolder/pairs_all $params.outdir/$outputFolder/all_pairs
        cp -rL ../../../$outputFolder/s4t_merged_pairsam $params.outdir/$outputFolder/cis_trans_pairs
        python ${baseDir}/bin/check_duplications_trans_cis_dist.py --inputdir $params.outdir/$outputFolder/qc_and_stats/stats --resultsdir $params.outdir/$outputFolder/qc_and_stats/qc
        python ${baseDir}/bin/count_double_labeled.py --inputdir ../../../$outputFolder --resultsdir $params.outdir/$outputFolder/qc_and_stats/qc
        """
}
