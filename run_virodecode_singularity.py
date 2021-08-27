# =================================================
# virodecode
# Version: 1.0 using Singularity
#
# Author: Sisi Huang
# Made at: Jiang Lab
# Github: https://github.com/Hc1023/virodecode
#
# =================================================

__author__ = 'sisihuang'
import sys, os
import argparse

## Set up arguments

parser = argparse.ArgumentParser()

# python run_virodecode_singularity.py -i virodecode_1.sif -1 fq1 -2 fq2
# Required Arguments
parser.add_argument("-i", "--image", type=str, required=True,
                    help="iamge path")
parser.add_argument("-1", "--read1", type=str, required=True, 
                    help="Provide read1 file path")
parser.add_argument("-2", "--read2", type=str, required=True,
                    help="Provide read2 file path")


# NOT Required Arguments
parser.add_argument("-o", "--outDir", type=str, default=".", required=False,
                    help="Provide path to a directory to place the output")
parser.add_argument("-r", "--ref", type=str, default="", required=False,
                    help="Reference genome")
parser.add_argument("-g", "--gff", type=str, default="", required=False,
                    help="gff file")
parser.add_argument("-p", "--para", type=str, default="sc", required=False,
                    help="spades.py parameter")
parser.add_argument("-c", "--coverage", type=str, default="5000", required=False,
                    help="subsample reads if the mean coverage is higher than xxx")
parser.add_argument("-m", "--mode", type=str, default="3", required=False,
                    help="1:variants calling 2:assembly 3:variants calling and assembly")
parser.add_argument("-t", "--thread", type=str, default="8", required=False,
                    help="thread")
parser.add_argument("-G", "--genome_version", type=str, default="", required=False,
                    help="thread")

args = parser.parse_args()



# requirement: read1 and read2 in the same directory
image = os.path.abspath(args.image)
fq1 = os.path.basename(args.read1)
fq2 = os.path.basename(args.read2)
OutDir = os.path.abspath(args.outDir)
InDir = os.path.dirname(os.path.abspath(args.read1))
thread=args.thread
para=args.para
coverage=args.coverage
mode=args.mode

if mode == "1":
    run = "Script1.sh"
if mode == "2":
    run = "Script2.sh"
if mode == "3":
    run = "Script3.sh"

if args.ref == "":
    ref = "/ref/NC_045512.2.fna"
    gff = "/ref/sequence.gff3"
    genome_version = "NC_045512.2"
    os.system("singularity exec \
-B " + args.read1 + ":" + args.read1 + " \
-B " + args.read2 + ":" + args.read2 + " \
-B " + OutDir + ":" + OutDir + " \
" + image + " bash /workfolder/" + run + " \
" + InDir + " " + OutDir + " " + fq1 + " " + fq2 + " " 
+ ref + " " + thread + " " + para + " " + coverage + " " + gff + " " + genome_version)
else:
    os.system("singularity exec \
-B " + args.read1 + ":" + args.read1 + " \
-B " + args.read2 + ":" + args.read2 + " \
-B " + OutDir + ":" + OutDir + " \
-B " + args.ref + ":" + args.ref + " \
-B " + args.gff + ":" + args.gff + " \
" + image + " bash /workfolder/" + run + " \
" + InDir + " " + OutDir + " " + fq1 + " " + fq2 + " " 
+ args.ref + " " + thread + " " + para + " " + coverage + " " + args.gff + " " + args.genome_version)

# singularity installed