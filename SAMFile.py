#!/usr/bin/env python
"""
Parse something resembling the SAM file format.

This currently has a few caveats:
 - the fragplot preprocessing steps eliminate the header lines from SAM
   files; therefore, this module assumes those lines are already eliminated
 - the refseq field is assumed to look like:
     gi|228472346|ref|NZ_ACLQ01000011.1| 
   which is common in HMP accessions to SRA.  In this case, the string
   starting with "NZ_" is treated as the refseq
 - due to a preprocessing error, some of these SAM-like files include a
   first line starting with the string "Your job" that must be skipped 
"""

from __future__ import division
import sys, getopt

class SAMFragment:
    """
    Represents a fragment match to a reference genome.
    """

    def __init__(self, sam_string):
        fields = sam_string.rstrip().split()

        self.name = fields[0]
        self.refseq = fields[2].split('|')[3]
        self.location = int(fields[3])
        self.cigar = fields[5]
        self.sequence = fields[9]
        self.length = len(self.sequence)

        for field in fields[11:]:
            if field[0:5] == "MD:Z:":
                md = field[5:]
                # warning: finite state machine ahead!
                # state = 0: nothing to remember
                # state = 1: read a number, looking to see if next is a number
                # state = 2: read a "^", looking for A, T, C, and G for deletion
                state = 0
                tempstr = ''
                matches = 0
                mismatches = 0
                deletions = 0
                for char in md:
                    if char in ('0', '1', '2', '3', '4', '5', '6', '7',
                                '8', '9'):
                        if state == 2:
                            deletions += len(tempstr)
                            tempstr = ''
                        state = 1
                        tempstr += char
                    elif char == '^':
                        if state == 1:
                            matches += int(tempstr)
                            tempstr = ''
                        state = 2
                    elif char in ('A', 'T', 'C', 'G'):
                        if state == 1:
                            matches += int(tempstr)
                            tempstr = ''
                            state = 0
                            mismatches += 1
                        elif state == 2:
                            tempstr += char
                        else:
                            mismatches += 1
                if state == 1:
                    matches += int(tempstr)
                elif state == 2:
                    deletions += len(tempstr)

        self.matches = matches
        self.mismatches = mismatches
        self.deletions = deletions
        self.identity = 100 * (matches / (matches + mismatches))

    def __repr__(self):
        out = []
        out.append("Name:      %s\n" % self.name)
        out.append("Location:  %d\n" % self.location)
        out.append("Identity:  %4.1f\n" % self.identity)
        return ''.join(out)


def SAMFile(filename):
    """
    Parses a SAM format file, returning SAMFragment objects.

    This is a generator, so that there is no need to store a full list of
    fragments in memory if it's not necessary to do so.
    """

    f = open(filename, 'rb')

    for line in f:
        if not line.startswith("Your job"):  # preprocessing error on my part
            yield SAMFragment(line)

    f.close()


if __name__ == '__main__':
    opts, args = getopt.getopt(sys.argv[1:], 'h')

    for o, a in opts:
        if o == '-h':
            print 'Usage:'
            print '  %s <sam_file>' % sys.argv[0]
            sys.exit(0)
        else:
            print 'unhandled option'
            sys.exit(1)

    if len(args) == 0:
        print 'Specify a SAM file as an argument'
        sys.exit(1)

    for frag in SAMFile(args[0]):
        print frag

