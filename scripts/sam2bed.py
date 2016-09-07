from __future__ import print_function
from pysam import Samfile
from argparse import ArgumentParser
from sys import stdout, stdin

def parse_args():
    false_arg = dict(action='store_true', default=False)
    #parser.usage()
    parser = ArgumentParser(description="Convert SAM format to BED format (for"
                            " both paired-end and single-end data)")
    parser.add_argument('--uniq', '-u', **false_arg)
    parser.add_argument('--use-RNA-strand', '-r', **false_arg)
    parser.add_argument('--separate-bed', '-s', **false_arg)
    parser.add_argument('--verbose', '-v', **false_arg)
    parser.add_argument('in_sam', help="SAM or BAM formatted file")
    parser.add_argument('out_bed1')
    parser.add_argument('out_bed2', nargs='?', default=None)
    args = parser.parse_args()

    args.in_sam = Samfile('-', 'r') if args.in_sam == '-' else Samfile(args.in_sam)
    assert args.out_bed1 != args.out_bed2
    args.out_bed1 = stdout if args.out_bed1 == '-' else open(args.out_bed1, 'w')
    if args.out_bed2:
        args.out_bed2 = stdout if args.out_bed2 == '-' else open(args.out_bed2, 'w')
        args.separate_bed = True

    return args

def sam_to_bed(read, use_rna_strand=False):
    if read.is_unmapped:
        return 0

    strand = '-' if read.is_reverse else '+'
    if use_rna_strand:
        try:
            strand = read.get_tag('XS')
        except KeyError:
            pass
    try:
        score = read.get_tag('NM')
    except KeyError:
        score = 0

    chrom_start = read.pos
    curr_len = 0
    extend_block = False
    block_sizes = []
    block_starts = []
    for cig_type, cig_len in read.cigartuples:
        if cig_type > 3:
            continue
        if cig_type == 1 or cig_type == 2: #I or D
            extend_block = True
            if cig_type == 2: #D
                if block_sizes:
                    block_sizes[-1] += cig_len
                    curr_len += cig_len
                else:
                    chrom_start += cig_len
        elif cig_type == 0: #M
            if extend_block and block_sizes:
                block_sizes[-1] += cig_len
            else:
                block_sizes.append(cig_len)
                block_starts.append(curr_len)
            extend_block = False
        curr_len += cig_len

    chrom_end = chrom_start + block_starts[-1] + block_sizes[-1] -1
    return dict(
        chrom=read.reference_name,
        chrom_start=chrom_start,
        chrom_end=chrom_end,
        name=read.qname,
        score=score,
        strand=strand,
        thick_start=chrom_start,
        thick_end=chrom_end,
        block_count=len(block_sizes),
        itemRgb=0,
        block_sizes=','.join(str(b) for b in block_sizes),
        block_starts=','.join(str(b) for b in block_starts),
    )

def bed_to_line(bed):
    colnames = [
        'chrom', 'chrom_start', 'chrom_end', 'name', 'score', 'strand',
        'thick_start', 'thick_end', 'itemRgb', 'block_count', 'block_sizes',
        'block_starts'
    ]
    outstr = '\t'.join(str(bed[key]) for key in colnames)
    return outstr

if __name__ == "__main__":
    args = parse_args()

    for read in args.in_sam:
        if read.is_unmapped: continue
        bed = sam_to_bed(read, args.use_RNA_strand)

        try:
            uniq = (read.get_tag('XT') == 'U')
        except KeyError:
            uniq = False

        if not args.uniq or uniq:
            if args.separate_bed and read.is_read2:
                print(bed_to_line(bed), file=args.out_bed2)
            else:
                print(bed_to_line(bed), file=args.out_bed1)
    args.out_bed1.close()
    if args.out_bed2:
        args.out_bed2.close()





