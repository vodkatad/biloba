#!/usr/bin/env python
from __future__ import with_statement

from sys import stdin, stderr
from optparse import OptionParser

class CNV:
    def __init__(self, chr, begin, end, cn):
        self.chr = chr
        self.begin = int(begin)
        self.end = int(end)
        self.cn = int(cn)

    # return a tuple (bool, integer)
    # if bool is true then the intereg represents of bases in common between self and other
    # if bool is false it's the integer represents their relative position: self.begin - other.begin (XXX reason corner cases)
    # float('nan') for different chr
    def relative_pos(self, other):
        #python check class TODO
        if (self.chr == other.chr):
            if (self.begin >= other.begin and self.begin < other.end) or (self.end > other.begin and self.end <= other.end) or (self.begin < other.begin and self.end >= other.end):
                start_ov = self.begin
                end_ov = self.end
                if self.begin < other.begin:
                    start_ov = other.begin
                if self.end > other.end:
                    end_ov = other.end
                return(True, end_ov - start_ov)
            else:
                return(False, self.begin - other.begin)
        else:
            return (False, float('nan'))

    def length(self):
        return self.end-self.begin
    
    # TODO rename with ide
    def get_cn(self):
        #print("fuck")
        return self.cn


UPPER = 10
MIN_OV = 0.5
class CNV_equivalence:
    def __init__(self, upper_cn, min_ov):
        self.UPPER = upper_cn
        self.MIN_OV = min_ov

    def compare_cn(self, cn_a, cn_b):
        #print("comparing {:s} {:s}".format(cn_a, cn_b));
        #print("comparing {:d} {:d}".format(cn_a.get_cn(), cn_b.get_cn()));
        if cn_a.get_cn() == cn_b.get_cn():
            return True
        elif cn_a.get_cn() > self.UPPER and cn_b.get_cn() > self.UPPER:
            return True
        else:
            return False
 
    def overlap(self, len_a, len_b, ov):
        return (ov / (len_a+len_b-ov)) >= self.MIN_OV


def main():
    usage = '''
            jaccard_cnv -o min_ov -m max_cnv cnv1 cnv2
            '''
    parser = OptionParser(usage=usage)
    parser.add_option('-o', '--minov', type=float, dest='min_ov', default=MIN_OV, help='minimum requested overlap to consider CNVs == [default:'+ str(MIN_OV) + ']')
    parser.add_option('-m', '--max_cnv', dest='upper', default=UPPER, help='upper level to consider cn == [default:' + str(UPPER) + ']')

    options, args = parser.parse_args()

    if len(args) != 2:
        exit('Unexpected argument number.\n' + usage)

    tot_cnv1 = 0
    tot_cnv2 = 0
    common_cnvs = 0
    config = CNV_equivalence(options.upper, options.min_ov)
    with file(args[0], 'r') as fcnv1:
        with file(args[1], 'r') as fcnv2:
            done = False
            line1 = fcnv1.readline()
            line2 = fcnv2.readline()
            while (not done):
                #print("infinite loop {:s} {:s}".format(line1, line2));
                skip = False
                if line1:   
                    tot_cnv1 += 1  
                else:
                    line2 = fcnv2.readline()
                    skip = True
                    # we cycle to count all cnvs but do not look for overlaps between smt and nothing
                if line2:
                    tot_cnv2 += 1
                else:
                    line1 = fcnv1.readline()
                    skip = True
                if not line1 and not line2:
                    done = True
                    break
                if skip:
                    continue
                cnvs1 = line1.rstrip().split('\t')
                cnvs2 = line2.rstrip().split('\t')
                cnv1 = CNV(cnvs1[0], cnvs1[1], cnvs1[2], cnvs1[4])
                cnv2 = CNV(cnvs2[0], cnvs2[1], cnvs2[2], cnvs2[4])
                ov = cnv1.relative_pos(cnv2)
                if ov[0] and config.overlap(cnv1.length(), cnv2.length(), ov[1]):
                    #print("almost found common! {:s} {:s}".format(cnvs1[1], cnvs2[1]));
                    if (config.compare_cn(cnv1, cnv2)):
                        #print("found common! {:s} {:s}".format(cnvs1[1], cnvs2[1]));
                        common_cnvs += 1
                    line1 = fcnv1.readline()
                    line2 = fcnv2.readline() # we consider the first overlap that we find...lots of corner cases
                else:
                    #print("not overlapping {:s} {:s}".format(cnvs1[1], cnvs2[1]));
                    # refactor - give ordering of chr also in Cnv class after calling relative_pos
                    if ov[1] > 0: # cnv1.begin - cnv2.begin
                        line2 = fcnv2.readline() # cnv1 is after cnv2, we go on reading from file2
                    elif ov[1] < 0: # check no == 0 corner case TODO
                        line1 = fcnv1.readline() # cnv1 is before cnv2, we go on reading from file1
                    else:
                        # different chrs, we go on the smaller one
                        chr1 = line1[0];
                        chr2 = line2[0];
                        if (chr1 > chr2): # XXX ensure same ordering with sort!
                            line2 = fcnv2.readline()
                        else:
                            line1 = fcnv1.readline()
            print('{:d}\t{:d}\t{:d}\t{:f}'.format(common_cnvs, tot_cnv1, tot_cnv2, common_cnvs/(tot_cnv1+tot_cnv2)))


if __name__ == '__main__':
    main()
