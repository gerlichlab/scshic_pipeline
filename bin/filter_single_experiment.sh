SAMPLE=$1
OUTPUT_FOLDER=$2
INPUT_S4T_PAIRS=$3

MIN_MAPQ=$4
MIN_RIGHT_MUTS=$5
MAX_WRONG_MUTS=$6

FILTER_SUMMARY="mapq_${MIN_MAPQ}.phred_30.right_muts_${MIN_RIGHT_MUTS}.wrong_muts_${MAX_WRONG_MUTS}"

MAPQ_FILTER="((mapq1 >= ${MIN_MAPQ}) and (mapq2 >= ${MIN_MAPQ}))"

SIDE1_HAS_AG="(int(n_AG_muts_phred30_1) >= ${MIN_RIGHT_MUTS})" 
SIDE1_HAS_TC="(int(n_TC_muts_phred30_1) >= ${MIN_RIGHT_MUTS})" 
SIDE2_HAS_AG="(int(n_AG_muts_phred30_2) >= ${MIN_RIGHT_MUTS})" 
SIDE2_HAS_TC="(int(n_TC_muts_phred30_2) >= ${MIN_RIGHT_MUTS})" 
SIDE1_NO_AG="(int(n_AG_muts_phred30_1) <= ${MAX_WRONG_MUTS})" 
SIDE1_NO_TC="(int(n_TC_muts_phred30_1) <= ${MAX_WRONG_MUTS})" 
SIDE2_NO_AG="(int(n_AG_muts_phred30_2) <= ${MAX_WRONG_MUTS})" 
SIDE2_NO_TC="(int(n_TC_muts_phred30_2) <= ${MAX_WRONG_MUTS})" 

SIDE1_IS_REF="(${SIDE1_HAS_TC} and ${SIDE1_NO_AG})"
SIDE2_IS_REF="(${SIDE2_HAS_TC} and ${SIDE2_NO_AG})"
SIDE1_IS_COMP="(${SIDE1_HAS_AG} and ${SIDE1_NO_TC})"
SIDE2_IS_COMP="(${SIDE2_HAS_AG} and ${SIDE2_NO_TC})"

IS_CIS_REF="(${SIDE1_IS_REF} and ${SIDE2_IS_REF})"
IS_CIS_COMP="(${SIDE1_IS_COMP} and ${SIDE2_IS_COMP})"
IS_TRANS_REF_COMP="(${SIDE1_IS_REF} and ${SIDE2_IS_COMP})"
IS_TRANS_COMP_REF="(${SIDE1_IS_COMP} and ${SIDE2_IS_REF})"
        
IS_CIS="(${IS_CIS_REF} or ${IS_CIS_COMP})"
IS_TRANS="(${IS_TRANS_REF_COMP} or ${IS_TRANS_COMP_REF})"

echo "filtering..."

pairtools select "(${MAPQ_FILTER} and ${IS_CIS_REF})" ${INPUT_S4T_PAIRS} \
    -o ${OUTPUT_FOLDER}/${SAMPLE}.${FILTER_SUMMARY}.cis_ref.pairs.gz

pairtools select "(${MAPQ_FILTER} and ${IS_CIS_COMP})" ${INPUT_S4T_PAIRS} \
    -o ${OUTPUT_FOLDER}/${SAMPLE}.${FILTER_SUMMARY}.cis_comp.pairs.gz
    
pairtools select "(${MAPQ_FILTER} and ${IS_TRANS_REF_COMP})" ${INPUT_S4T_PAIRS} \
    -o ${OUTPUT_FOLDER}/${SAMPLE}.${FILTER_SUMMARY}.trans_ref_comp.pairs.gz
    
pairtools select "(${MAPQ_FILTER} and ${IS_TRANS_COMP_REF})" ${INPUT_S4T_PAIRS} \
    -o ${OUTPUT_FOLDER}/${SAMPLE}.${FILTER_SUMMARY}.trans_comp_ref.pairs.gz
