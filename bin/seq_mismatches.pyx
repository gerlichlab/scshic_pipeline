from cpython cimport array
import cython
cimport cython

cpdef list get_mismatches_c(str seq, array.array quals, list aligned_pairs):
    '''
    This function takes a SAM alignment and, for every mismatch between the read and reference sequences,
    returns a tuple (the reference bp, the read bp, PHRED quality of the bp).
    '''

    cdef cython.int read_pos, ref_pos
    cdef str orig_bp, orig_bp_upper
    cdef list mismatches = []
    
    for read_pos, ref_pos, orig_bp in aligned_pairs:
        orig_bp_upper = orig_bp.upper()
        if (seq[read_pos]!=orig_bp_upper):

            mismatches.append(
                (orig_bp_upper,
                 seq[read_pos],
                 quals[read_pos])
            )
        
    return mismatches 
