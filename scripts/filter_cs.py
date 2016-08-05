import sys
from pysam import AlignmentFile
from argparse import ArgumentParser

valid_spliced_reads=0
problem_reads=0

parser = ArgumentParser()
parser.add_argument('infile', nargs='?', default='-')
parser.add_argument('outfile', nargs='?', default='-')
args = parser.parse_args()
infile = AlignmentFile(args.infile, 'r')
outfile = AlignmentFile(args.outfile, 'wh', template=infile)

for read in infile:
    splice_len = 0
    min_edge = 1e6
    if read.mapping_quality < 10: continue
    for cig_op, cig_len in read.cigartuples:
        if cig_op == 3: # N
            splice_len += cig_len
        elif cig_op == 0:
            min_edge = min(min_edge, cig_len)
    if splice_len > 50 and min_edge >= 6:
        outfile.write(read)
        valid_spliced_reads += 1
        if valid_spliced_reads % 100000 == 0:
            sys.stderr.write("%d valid, %d problematic spliced reads\n" % (valid_spliced_reads, problem_reads) )
sys.stderr.write("%d valid, %d problematic spliced reads\n" % (valid_spliced_reads, problem_reads) )

