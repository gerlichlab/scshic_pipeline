#!/usr/bin/env python

from seq_mismatches import get_mismatches_c
import sys
import pandas as pd
import pysam
import click
import pairtools
import pyximport
pyximport.install()


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('PAIRSAM_PATH', nargs=1)
@click.option(
    '-o',
    '--output',
    type=str,
    default=None,
    show_default=True,
    help='path to the output pairs file (by default, printed into stdout)')
@click.option(
    '--chunksize',
    type=int,
    default=int(1e6),
    show_default=True,
    help='the number of pairs per processed chunk')
def detect_s4t_mutations(
    pairsam_path,
    output,
    chunksize
):
    """Add columns with the number of s4T-induced mutations.
    """
    fill_n_s4t_pairs(pairsam_path, output, chunksize=chunksize)


def fill_n_s4t_pairs(in_path, out_path, chunksize=int(1e6)):
    header, pairs_body = pairtools._headerops.get_header(
        pairtools._fileio.auto_open(in_path, 'r'))

    cols = pairtools._headerops.extract_column_names(header)

    new_cols = ([col for col in cols if not col.startswith('sam')]
                + ['n_AG_muts_phred0_1', 'n_AG_muts_phred0_2', 'n_AG_muts_phred30_1', 'n_AG_muts_phred30_2',
                   'n_TC_muts_phred0_1', 'n_TC_muts_phred0_2', 'n_TC_muts_phred30_1', 'n_TC_muts_phred30_2',
                   'n_matches_1', 'n_matches_2'
                   ])
    new_header = pairtools._headerops.extract_fields(
        header, 'columns', save_rest=True)[1]
    new_header.append('#columns: '+' '.join(new_cols))

    samheader_lines = pairtools._headerops.extract_fields(header, 'samheader')
    samheader = pysam.AlignmentHeader.from_text('\n'.join(samheader_lines))

    out_stream = (
        pairtools._fileio.auto_open(out_path, nproc=3, mode='w')
        if out_path
        else sys.stdout)

    out_stream.writelines([line+'\n' for line in new_header])

    chunk_iterator = pd.read_csv(
        pairs_body,
        header=None,
        names=cols,
        chunksize=chunksize,
        sep='\t'
    )

    for chunk_df in chunk_iterator:
        fill_n_s4t_df(chunk_df, samheader)

        chunk_df.to_csv(
            out_stream,
            sep='\t',
            columns=new_cols,
            header=False,
            index=False,
        )

    out_stream.close()


def fill_n_s4t_df(pairs_df, samheader):
    pairs_df['n_AG_muts_phred0_1'] = 0
    pairs_df['n_AG_muts_phred0_2'] = 0
    pairs_df['n_AG_muts_phred30_1'] = 0
    pairs_df['n_AG_muts_phred30_2'] = 0

    pairs_df['n_TC_muts_phred0_1'] = 0
    pairs_df['n_TC_muts_phred0_2'] = 0
    pairs_df['n_TC_muts_phred30_1'] = 0
    pairs_df['n_TC_muts_phred30_2'] = 0

    pairs_df['n_matches_1'] = 0
    pairs_df['n_matches_2'] = 0

    for i in pairs_df.index:
        for side in ['1', '2']:
            mms, n_matches = get_mismatches_samtxt(pairs_df.at[i, 'sam'+side],
                                                   samheader=samheader)
            n_AG_muts_phred0 = 0
            n_AG_muts_phred30 = 0
            n_TC_muts_phred0 = 0
            n_TC_muts_phred30 = 0
            for mm in mms:
                if (mm[0] == 'A' and mm[1] == 'G'):
                    n_AG_muts_phred0 += 1
                    if mm[2] >= 30:
                        n_AG_muts_phred30 += 1
                elif (mm[0] == 'T' and mm[1] == 'C'):
                    n_TC_muts_phred0 += 1
                    if mm[2] >= 30:
                        n_TC_muts_phred30 += 1

            pairs_df.at[i, 'n_AG_muts_phred0_'+side] = n_AG_muts_phred0
            pairs_df.at[i, 'n_AG_muts_phred30_'+side] = n_AG_muts_phred30

            pairs_df.at[i, 'n_TC_muts_phred0_'+side] = n_TC_muts_phred0
            pairs_df.at[i, 'n_TC_muts_phred30_'+side] = n_TC_muts_phred30
            pairs_df.at[i, 'n_matches_'+side] = n_matches


def get_mismatches_samtxt(sam_col, samheader):
    sam_strs = [sam_str.replace(pairtools._pairsam_format.SAM_SEP, '\t')
                for sam_str in sam_col.split(pairtools._pairsam_format.INTER_SAM_SEP)
                ]
    sam_algns = [pysam.AlignedSegment.fromstring(sam_str, samheader)
                 for sam_str in sam_strs]
    algn = sorted(sam_algns, key=lambda _: (_.query_length - _.query_alignment_end
                                            if _.is_reverse else _.query_alignment_start))[0]
    if not algn.has_tag('MD'):
        return [], algn

    seq = algn.query_sequence.upper()
    quals = algn.query_qualities
    aligned_pairs = algn.get_aligned_pairs(with_seq=True, matches_only=True)
    # aligned_pairs = algn.get_aligned_pairs(with_seq=True, matches_only=False)
    # aligned_pairs = [i for i in aligned_pairs if i[0] is not None and i[2] is not None]
    mismatches = get_mismatches_c(seq, quals, aligned_pairs)
    n_matches = len(aligned_pairs)

    return mismatches, n_matches


if __name__ == '__main__':
    detect_s4t_mutations()
