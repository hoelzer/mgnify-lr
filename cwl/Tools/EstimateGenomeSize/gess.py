#!/usr/bin/env python3

#docstring
"""
gess: a tool to estimate genome size from Illumina reads
Uses k-mer counting: requires kmc
input: illumina fastq reads - R1 and R2
output: estimated genome size
"""


#global variables
MIN_KMER = 16
MAX_KMER = 127
MIN_BINS = 10
MAX_BINS = 300
MIN_CUTOFF = 3
MAX_CUTOFF = 20

#import modules
import subprocess
import sys
import os
import logging
import re #regular expressions
import argparse
from argparse import RawTextHelpFormatter


def main():



    # make parser, add arguments, read in args from user
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description='estimates genome size', usage='\n %(prog)s args.fastq_file', prog='gess')
    parser.add_argument('--version', '-v', action='version', version='%(prog)s 0.1')
    parser.add_argument("fastq_file")
    parser.add_argument('--kmer', '-k', type=int, default=25, help='kmer size for KMC') 
    #kmer size can affect the estimate; depends on noise and depth 
    parser.add_argument('--nbins', '-n', type=int, default=200, help='number of bins')  
    parser.add_argument('--cutoff', '-ci', type=int, default=3, help='discard kmers if frequency < cutoff')
    parser.add_argument('--tmpdir', '-t', default='/tmp')
    parser.add_argument('--quiet', '-q', action="store_true") #if user specifies this arg, it has value True
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--citation', action="store_true")
    parser.add_argument('--check') ##TODO - check kmc runs
    args = parser.parse_args()

    #set logging level
    if args.quiet:
        logging.basicConfig(level=logging.WARNING)
    else:
        logging.basicConfig(level=logging.DEBUG)

    #citation
    if args.citation:
        print("Citation: Gess, a tool for genome size estimation. GitHub")


    #check input is correct 
    if args.kmer<MIN_KMER or args.kmer>MAX_KMER:
        logging.error("kmer value should be between {} and {}".format(MIN_KMER, MAX_KMER)) 
        #TODO find out kmc's allowed kmers
        sys.exit(1)
    if args.nbins<MIN_BINS or args.nbins>MAX_BINS: #TODO check suitable bin sizes 
        logging.error("nbins should be between {} and {}".format(MIN_BINS, MAX_BINS))
        sys.exit(1)
    if args.cutoff<MIN_CUTOFF or args.cutoff>MAX_CUTOFF: #TODO check suitable cutoffs
        logging.error("cutoff should be between {} and {}".format(MIN_CUTOFF, MAXCUTOFF))
        sys.exit(1)
    if args.threads < 1:
        logging.error("threads should be greater than 0")
        sys.exit(1)

    # use process ID to name temp files
    pid = str(os.getpid())
    logging.info("the pid is: " + pid)


    #TODO: check that kmc works
    #get the version of kmc
    #kmc_cmd = ["kmc"]
    #kmc_check = subprocess.check_output(kmc_cmd)
    #kmc_pattern = re.search(r'ver.\s(d+)', kmc_check)
    # check kmc runs at all
    # try / except


    # talk to the user
    logging.info("You provided this fastq file: {}". format(args.fastq_file))
    logging.info("Now running KMC with these commands:")
    cmd = ["kmc", "-k"+str(args.kmer), "-n"+str(args.nbins), "-ci"+str(args.cutoff), "-t"+str(args.threads), args.fastq_file, pid, args.tmpdir]
    logging.debug(cmd)
    # -ci3 option discards kmers with freq <3
    # -n200 sets number of bins to 200 because you can only have 256 files open at once
    # could also change this with ulimit
    #/tmp is the working directory

    # check that kmc ran properly
    try:
        result = subprocess.check_output(cmd) #this makes a bytes thing
        #subprocess.check_output is equivalent to subprocess.run with the default arg check=True
    except:
        logging.exception("error running kmc")
        sys.exit(1)

    #delete unnecessary files
    logging.info("removing unnecessary files")
    os.remove(pid + ".kmc_pre")
    os.remove(pid + ".kmc_suf")

    #Report the output
    #print(type(result))
    #print(result)
    converted = str(result, 'utf-8') #converts bytes to string
    #print(type(converted))
    #print(converted)
    #No. of unique counted k-mers       :      2784627
    pattern = re.search(r'unique counted k-mers\s*:\s*(\d+)', converted)
    #use re.search not re.match
    genome_size = (pattern.groups())
    #print(genome_size[0])  #extract first item of tuple
    #can exclude this step now
    #split = converted.splitlines() # a list of strings
    #print(split)
    #for line in split:
    #    if "unique counted" in line:
    #        #print(line)
    #        numbers = line
    #        numbers = ''.join(n for n in numbers if n in '0123456789')
    genome_size_inMbp = int(genome_size[0])/1_000_000
    logging.info("Estimated genome size is {0:.2f} Mbp".format(genome_size_inMbp))
    print("Estimated genome size is {0:.2f} Mbp".format(genome_size_inMbp))
    #print(genome_size[0]) #TODO make stdout?



if __name__=="__main__":
    main()
