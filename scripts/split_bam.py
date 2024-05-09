#!/usr/bin/env python3

import sys, os
import argparse
import pysam

from TRHelper.Aux import read_genome_sequence_lengths, split_bam

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Split BAM by read type: discordant pairs, proper pairs, '
        'clipped, supplementary')
parser.add_argument('input', help='input file')
parser.add_argument('-c', '--min-clip-length', type=int, default=1,
    help='minimal length of the clipped end of the read')
parser.add_argument('-d', '--dist-disc', type=int, default=1e6,
    help='minimal distance for discordant reads on the same chr')
parser.add_argument('-D', '--dist-conc', type=int, default=1e4,
    help='maximal distance for concordant reads on the same chr')
parser.add_argument('-o', '--out-basename', help='basename for the output '
    'files, default: basename of the input file')
parser.add_argument('-q', '--min-end-qual', type=int, default=3,
    help='minimal quality of the end base of the clipped region')
parser.add_argument('-r', '--reference-fa-gz',
    help='fa.gz reference file')
parser.add_argument('-t', '--threads', type=int, default=1,
    help='number of threads to use')
parser.add_argument('-w', '--write-drp', action='store_true',
    help='identify and write DRPs to .disc.bam')


args = parser.parse_args()
out_bn = (args.out_basename if args.out_basename is not None 
          else os.path.splitext(args.input)[0])


# open the input/output sam files
mode = 'r'
input_fn_extension = os.path.splitext(args.input)[1].lower()
if input_fn_extension == "bam":
    mode = 'rb'
elif input_fn_extension == "cram":
    mode = 'rc'

infile = pysam.AlignmentFile(args.input, mode, reference_filename =
    args.reference_fa_gz if args.reference_fa_gz else None,
    threads=args.threads)
outfile_disc = pysam.AlignmentFile(out_bn + ".disc.bam" if args.write_drp
    else os.devnull, "wb", template=infile, threads=args.threads)
outfile_clip = pysam.AlignmentFile(out_bn + ".clip.bam", 
    "wb", template=infile, threads=args.threads)

# check that the reference file matches input
if args.reference_fa_gz:
    glh = read_genome_sequence_lengths(args.reference_fa_gz)
    for rname in infile.references:
        try:
            rlen = glh[rname]
        except KeyError:
            sys.stdout.write(f"WARNING: contig {rname} is missing from "
                f"the reference file {args.reference_fa_gz}\n")
            continue
        if rlen != infile.get_reference_length(rname):
            raise ValueError(f"Non-matching contig {rname} length in "
                f"{args.reference_fa_gz}: {rlen} vs "
                f"{infile.get_reference_length(rname)} in {args.input}")

# iterate through alignments in the input file
nBases, tlen_dist, rlen_dist = split_bam(
    infile, outfile_clip, outfile_disc, args)

outfile_disc.close()
outfile_clip.close()

#write stats
with open(out_bn + ".nbases", "w") as fid:
    fid.write("{}\n".format(nBases))

with open(out_bn + ".tlen", "w") as fid:
    # TODO: explain what we use this info for
    if tlen_dist:
        for k in sorted(tlen_dist.keys()):
            fid.write("{}\t{}\n".format(k, tlen_dist[k]))
    else:
        raise ValueError("No properly paired reads; SE reads perhaps?")

with open(out_bn + ".rlen", "w") as fid:
    for k in sorted(rlen_dist.keys()):
        fid.write("{}\t{}\n".format(k, rlen_dist[k]))
