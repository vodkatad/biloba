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
        return self.cn


UPPER = 10
MIN_OV = 0.5
class CNV_equivalence:
    def __init__(self, upper_cn, min_ov):
        self.UPPER = upper_cn
        self.MIN_OV = min_ov

    def compare_cn(self, cn_a, cn_b):
        #print("comparing {:s} {:s}".format(cn_a, cn_b));
        #print("comparing {:d} {:d} {:d}".format(cn_a.get_cn(), cn_b.get_cn(), self.UPPER));
        if cn_a.get_cn() == cn_b.get_cn():
            return True
        elif cn_a.get_cn() > self.UPPER and cn_b.get_cn() > self.UPPER:
            return True
        else:
            return False
 
    def overlap(self, len_a, len_b, ov):
        #print("{:d} {:d} {:d} {:.15f}".format(len_a, len_b, ov, float(ov) / (len_a+len_b-ov)))
        return (float(ov) / (len_a+len_b-ov)) >= self.MIN_OV # ouch float!

def main():
    usage = '''
            jaccard_cnv -o min_ov -m max_cnv cnv1 cnv2
            '''
    parser = OptionParser(usage=usage)
    parser.add_option('-o', '--minov', type=float, dest='min_ov', default=MIN_OV, help='minimum requested overlap to consider CNVs == [default:'+ str(MIN_OV) + ']')
    parser.add_option('-m', '--max_cnv', type=int, dest='upper', default=UPPER, help='upper level to consider cn == [default:' + str(UPPER) + ']')

    options, args = parser.parse_args()

    if len(args) != 2:
        exit('Unexpected argument number.\n' + usage)

    config = CNV_equivalence(options.upper, options.min_ov)

    cells1 = {}
    cells2 = {}

    def load_cells(dictio, toload, filename):
        with file(filename, 'r') as fh:
            for line in fh:
                l = line.rstrip().split('\t')
                if len(toload.keys()) == 0 or toload.get(l[3]):
                    if dictio.get(l[3]):
                        dictio[l[3]].append(l)
                    else:
                        dictio[l[3]] = [l]

    load_cells(cells1, {}, args[0])
    load_cells(cells2, cells1, args[1])
    #print(cells1)
    #print(cells2)

    def next_line_count(cns, n):
        if n+1 <= len(cns):
            return((cns[n], n+1))
        else:
            return((None, n))

    def process_cell(fcnv1, fcnv2, config):
        tot_cnv1 = 0
        tot_cnv2 = 0
        common_cnvs = 0
        done = False
        (line1, tot_cnv1) = next_line_count(fcnv1, tot_cnv1)
        (line2, tot_cnv2) = next_line_count(fcnv2, tot_cnv2)
        while (not done):
            #print("infinite loop {:s} {:s}".format(line1, line2))
            skip = False
            if not line1 and line2: 
                (line2, tot_cnv2) = next_line_count(fcnv2, tot_cnv2)
                skip = True
                # we cycle to count all cnvs but do not look for overlaps between smt and nothing
            if not line2 and line1:
                (line1, tot_cnv1) = next_line_count(fcnv1, tot_cnv1)
                skip = True
            if not line1 and not line2:
                done = True
                break
            if skip:
                continue
            cnvs1 = line1#.rstrip().split('\t')
            cnvs2 = line2#.rstrip().split('\t')
            cnv1 = CNV(cnvs1[0], cnvs1[1], cnvs1[2], cnvs1[4])
            cnv2 = CNV(cnvs2[0], cnvs2[1], cnvs2[2], cnvs2[4])
            ov = cnv1.relative_pos(cnv2)
            if ov[0] and config.overlap(cnv1.length(), cnv2.length(), ov[1]):
                #print("almost found common! {:s} {:s}".format(cnvs1[1], cnvs2[1]));
                if (config.compare_cn(cnv1, cnv2)):
                    #print("found common! {:s} {:s}".format(cnvs1[1], cnvs2[1]));
                    common_cnvs += 1
                    (line1, tot_cnv1) = next_line_count(fcnv1, tot_cnv1)
                    (line2, tot_cnv2) = next_line_count(fcnv2, tot_cnv2)
                else:
                    #print("not overlapping {:s} {:s}".format(cnvs1[1], cnvs2[1]));
                    # refactor - give ordering of chr also in Cnv class after calling relative_pos
                    if ov[1] > 0: # cnv1.begin - cnv2.begin
                        (line2, tot_cnv2) = next_line_count(fcnv2, tot_cnv2)
                    elif ov[1] < 0: # check no == 0 corner case TODO
                        (line1, tot_cnv1) = next_line_count(fcnv1, tot_cnv1)
            else:
                # different chrs, we go on the smaller one
                chr1 = line1[0];
                chr2 = line2[0];
                if (chr1 > chr2): # XXX ensure same ordering with sort!
                    (line2, tot_cnv2) = next_line_count(fcnv2, tot_cnv2)
                else:
                    (line1, tot_cnv1) = next_line_count(fcnv1, tot_cnv1)
        print('{:d}\t{:d}\t{:d}\t{:f}'.format(common_cnvs, tot_cnv1, tot_cnv2, float(common_cnvs)/(tot_cnv1+tot_cnv2-common_cnvs)))

    for k in cells1.keys():
        if cells2.get(k):
            process_cell(cells1[k], cells2[k], config)

if __name__ == '__main__':
    main()
