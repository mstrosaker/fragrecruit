#!/usr/bin/env python
"""
Parse the fr-hit output file format.
"""

from __future__ import division
import sys, getopt

class FRHITFragment:
    """
    Represents a fragment match to a reference genome.
    """

    def __init__(self, frhit_string):
        fields = frhit_string.rstrip().split()

        self.name = fields[0]
        self.refseq = fields[8].split('|')[3]
        self.location = int(fields[9])
        self.identity = float(fields[7][:-1])
        self.length = int(fields[3])

    def __repr__(self):
        out = []
        out.append("Name:      %s\n" % self.name)
        out.append("Location:  %d\n" % self.location)
        out.append("Identity:  %4.1f\n" % self.identity)
        return ''.join(out)


def FRHITFile(filename):
    """
    Parses the fr-hit format file, returning FRHITFragment objects.

    This is a generator, so that there is no need to store a full list of
    fragments in memory if it's not necessary to do so.
    """

    f = open(filename, 'rb')

    for line in f:
        yield FRHITFragment(line)

    f.close()


if __name__ == '__main__':
    opts, args = getopt.getopt(sys.argv[1:], 'h')

    for o, a in opts:
        if o == '-h':
            print 'Usage:'
            print '  %s <fr-hit_output_file>' % sys.argv[0]
            sys.exit(0)
        else:
            print 'unhandled option'
            sys.exit(1)

    if len(args) == 0:
        print 'Specify an FR-HIT output file as an argument'
        sys.exit(1)

    for frag in FRHITFile(args[0]):
        print frag

