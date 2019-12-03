#!/usr/bin/env python

import re
import argparse
import random
def get_args():
    """ get file names of sam file, name of the known umis, and the character to indicate what to do with pcr duplicates """
    parser = argparse.ArgumentParser(description=' File name to be used',add_help=False)
    parser.add_argument('-f', '--file', type=str, help='File Name of sorted Sam File', required=True)
    parser.add_argument('-u', '--umi', type=str, help='File Name of known umis if available', required=False, default=False)
    parser.add_argument('-p', '--paired', type=str, help='True if paired end reads are used, else False ', required=False, default=False)
    parser.add_argument('-k', '--keep', type=str, help='Defines how to choose which duplicate is kept. Options are highest quality(q), first in file(1), or random(r)',required=False, default=1)
    parser.add_argument('-h', '--help', action='help', help='a USEFUL help message(see argparse docs)', default=argparse.SUPPRESS)
    return parser.parse_args()
args = get_args()
f = args.file
k = args.keep
u = args.umi

def umi(f):
    """Takes a file of known UMI's and adds them to a list"""
    with open(f, 'r') as f1:
        umis = []
        for line in f1:
            line = line.strip()
            umis.append(line)
        umis = tuple(umis)
    return umis
# Test case: A file containing AAAAAAAA, output: ('AAAAAAAA')


umi = umi(u)


def dflag(f):
    """Decodes bitwise flag for strand, first in pair, and returns true or false for each needed bit """
    if ((f & 16) == 16):
        return True
# Test case: 16  output: True

def dcigar(c,p):
    """Takes a cigar string and leftmost position and returns rightmost position using regex"""
    x = re.sub("\dI","",c)
    # removes insertions
    y = re.split("[A-Z]", x)
    #removes letters and puts each number into a list
    y.pop(-1)
    # removes the empty item
    if re.findall("^\d*[a-z]",c) != None:
        y.pop(0)
    #removes softclipping if the cigar string starts with it
    for item in y:
        p += int(item)
        # Iterates through the list and adds each number to p returning the 5' start
    return p
#Test case: '3S66M2D3I24M2S', 3  output 97


def ec(umi,umis):
    "Checks to see an UMI that doesn't directly match can be error corrected, and corrects it"
    for u in umis:
        if 2 >= sum(a != b for a, b in zip(umi, u)):
            return u
    return None


def dehead(h,umis):
    """Uses sam headers, finds the 5' start sites, compares umi's and returns True if they match and false if not"""
    umi, flag, chrom, pos, cig = h.split("\t")[0], int(h.split("\t")[1]), h.split("\t")[2], int(h.split("\t")[3]), h.split("\t")[5]
    umi = umi.split(":")[7]
    if umis is not False:
        if dflag(flag) is True:
            pos = dcigar(cig,pos)
            flag = "t"
        elif re.search("^\d*S",cig):
            pos += int(re.split("[A-Z]",cig)[0])
            flag = "f"
        else:
            pos = pos
            flag = "f"
        if umi in umis:
            return (str(flag) + str(pos) + umi), str(chrom)
        else:
            ec1 = ec(umi, umis)
            if ec1 != None:
                return (str(flag) + str(pos) + ec1), str(chrom)
            else:
                return None, chrom
    else:
        if dflag(flag) is True:
            pos = dcigar(cig,pos)
            flag = "t"
        elif re.search("^\d*S",cig):
            pos += int(re.split("[A-Z]",cig)[0])
            flag = "f"
        else:
            pos = pos
            flag = "f"
    return (str(flag) + str(pos) + str(umi)), str(chrom)
#Test Case:
# 'NS500451:154:HWKTMBGXX:1:11101:8846:5235-CAACTGGT^CAGTACTG;0^0	83	2	108949756	36	66M	=	108949700	-122'
# 'NS500451:154:HWKTMBGXX:1:11101:8846:5235-CAACTGGT^CAGTACTG;0^0	83	2	108949756	36	66M	=	108949700	-122'
# output: True


def qual(seq1,seq2):
    """takes a quality score and checks if the mean phred score of the score line is greater

     q -- quality score cutoff
     seq -- encoded quality score line
     returns value -- true if mean score is greater and false if not"""
    sum1=0
    sum2=0
    for score in seq1:
        # sums all quality scores
        n = (((ord(score) - 33)))
        sum1 += n

    for score in seq2:
        # sums all quality scores
        n = (((ord(score) - 33)))
        sum2 += n
    if (int(sum1)/len(seq1))>=(int(sum2)/len(seq2)):
        return False
    else:
        return True


def dedupe(f1,umis,k):
    """Takes the initial Sam file, swaps leftmost position for 5' start and places leftmost at the end of the header
    f1 -- name of the original sam file
    returns new file """
    with open(f1,'r') as f, open(f1+"_deduped",'a') as o:
        chrom = None
        chrom1 = "n"
        reads = {}
        for line in f:
            if re.search("^@", line) is not None:
                o.write(line)
            else:
                head, chrom1 = dehead(line, umis)
                if chrom1 == chrom:
                    if head not in reads and head is not None:
                        reads[head] = line
                    elif head is None:
                        None
                    else:
                        if k == "q":
                            if qual(reads[head].split("\t")[10],line.split("\t")[10]) is True:
                                reads[head] = line
                        elif k == "r":
                            if random.random() > .5:
                                reads[head] = line
                        else:
                            reads[heads] = reads[heads]
                else:
                    for read in reads:
                        o.write(reads[read])
                    chrom = chrom1
                    reads = {head: line}
        for read in reads:
            o.write(reads[read])
# Testcase Sam File:
# Header 'NS500451:154:HWKTMBGXX:1:11101:8846:5235-CAACTGGT^CAGTACTG;0^0	83	2	108949756	36	66M	=	108949700	-122'
# Output 'NS500452:154:HWKTMBGXX:1:11101:8846:5235-CAACTGGT^CAGTACTG;0^0	83	2	108949822	36	66M	=	108949700	-122'


dedupe(f, umi, k)
