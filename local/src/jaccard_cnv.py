#!/usr/bin/env python
from __future__ import with_statement

from sys import stdin, stderr
from optparse import OptionParser

class CNV:
    def __init__(self, chr, begin, end, cn):
        self.chr = chr
        self.begin = begin
        self.end = end
        self.cn = cn

    # return a tuple (bool, integer)
    # if bool is true then the intereg represents of bases in common between self and other
    # if bool is false it's the integer represents their relative position: self.begin - other.begin (XXX reason corner cases)
    # float('nan') for different chr
    def relative_pos(self, other):
        #python check class TODO
        if (self.chr == other.chr):
            if (self.begin >= other.begin and self.begin < other.end) or
                (self.end > other.begin and self.end <= other.end) or
                (self.begin < other.begin and self.end >= other.end):
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

    def same_cn(self, other, config):
        if config.compare_cn(self.cn, other.cn): # reason about efficiency, which comparison should be done before? probably coords FIXME
            ov = self.relative_pos(other)
            if ov[0]:
                return config.overlap(self.length(), other.length(), ov[1])
            else:
                return False
        else:
            return False
        
    def length(self):
        return self.end-self.begin


UPPER = 10
MIN_OV = 0.5
class CNV_equivalence:
    def __init__(self, upper_cn, min_ov):
        self.UPPER = upper_cn
        self.MIN_OV = min_ov

    def compare_cn(cn_a, cn_b):
        if cn_a == cn_b:
            return True
        else if cn_a > self.UPPER and cn_b > self.UPPER:
            return True
        else:
            return False
 
    def overlap(len_a, len_b, ov):
        return (ov / (len_a+len_b-ov)) > self.MIN_OV


def main():
    usage = format_usage('''
        %prog PARAM1 PARAM2 < STDIN
    ''')
    parser = OptionParser(usage=usage)
    parser.add_option('-c', '--cutoff', type=float, dest='cutoff', default=0.05, help='some help CUTOFF [default: %default]', metavar='CUTOFF')
    parser.add_option('-p', '--dump-params', dest='params_file', help='some help FILE [default: %default]', metavar='FILE')
    parser.add_option('-v', '--verbose', dest='verbose', action='store_true', default=False, help='some help [default: %default]')

    options, args = parser.parse_args()

    if len(args) != 2:
        exit('Unexpected argument number.')

    #for id, sample, raw, norm in Reader(stdin, '0u,1s,2i,3f', False):
    with file(FILENAME, 'r') as fd:
        for line in fd:
            tokens = safe_rstrip(line).split('\t')
            #tokens = line.rstrip().split('\t')
            assert len(tokens) == NUM

if __name__ == '__main__':
    main()
