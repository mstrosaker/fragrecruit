#!/usr/bin/python

# Use cases:
#  - invoke as:
#       $ decode_sequences <filename>
#  - invoke with no arguments and pass input via stdin
#  - import as a module and call SequenceFile(filename)
#       Any file-like object can be passed, including a StringIO object
#
#  The input can be gzipped or bzip2'd, and must be in either FASTA or
#  FASTQ format.  Bzip2'd data can only be handled by the first use case
#  (passing a filename as command-line argument) because of limitations of
#  handling streaming bzip2 data in Python 2.x (fixed in 3.3; see
#  http://bugs.python.org/issue5863).
#
# Features and assumptions:
#  - filenames do not matter; gzip and bzip2 streams are autodetected, and
#    the input format (FASTA or FASTQ) is autodetected
#  - both DOS and UNIX newline sequences are handled
#  - the first sequence is assumed to start on the first line of the input
#  - no blank lines or leading whitespace may exist, and no whitespace may
#    exist in sequence or quality score strings
#  - the sequence and quality score strings may be spread across multiple
#    lines (in either format)

import sys, getopt, gzip, bz2, StringIO

class SequenceInputError(Exception):
    """Used to raise error conditions."""


class Sequence:
    """
    Represents a protein, DNA, or RNA sequence.

    Contains 3 members:
        name: string
        seq: string containing the full, unbroken sequence
        qual: array of ints, one per character in seq (FASTQ only)
    """

    def __init__(self, name, seq, qual=None):
        self.name = name
        self.seq = seq
        self.qual = qual

    def __repr__(self):
        out = []
        out.append("Name: %s\n" % self.name)
        out.append("Seq:  %s" % self.seq)
        if self.qual:
            out.append("\nQual:\n")
            for q in self.qual:
                out.append("%d " % q)
        return ''.join(out)


def supportCompression(stream):
    """
    Detect if a stream is compressed and modify to handle, if possible.
    """

    # gzip magic number: 0x1f8b
    # bzip magic number: 0x425a ('BZ')

    no_filename = False

    try:
        filename = stream.name
    except AttributeError:
        # probably a StringIO object
        no_filename = True

    if filename == '<stdin>':
        stream = StringIO.StringIO(sys.stdin.read())

    pos = stream.tell()
    magic = stream.read(2)
    stream.seek(pos)
    if magic == '\037\213':
        stream = gzip.GzipFile(None, 'rb', None, stream)
    elif magic == 'BZ':
        if no_filename:
            raise SequenceInputError('Cannot handle bzip2 compressed ' +
                                     'StringIO data.')
        if filename == '<stdin>':
            raise SequenceInputError('Cannot handle bzip2 compressed data ' +
                                     'over stdin.')
        stream = bz2.BZ2File(filename, 'rb')

    return stream


def SequenceFile(stream):
    """
    Parses FASTA or FASTQ data, returning Sequence objects.

    This is a generator, so that there is no need to store a full list of
    sequences in memory if the user does not want to.
    """

    f = supportCompression(stream)

    l = f.readline().rstrip()

    # autodetect which format; the first character knows all
    if l[0] == '>':
        format = 'fasta'
    elif l[0] == '@':
        format = 'fastq'
    else:
        raise SequenceInputError('Could not determine data type ' +
                                 '(fasta or fastq)')

    # l == '' indicates end-of-file
    while l != '':
        name = l[1:]
        seq_list = []
        qual_list = None

       # aannnndd... here comes the ugly parsing
        l = f.readline().rstrip()

        if format == 'fasta':
            while l != '' and l[0] != '>':
                seq_list.append(l)
                l = f.readline().rstrip()
            s = ''.join(seq_list)

        elif format == 'fastq':
            # read in the sequence
            while l != '' and l[0] != '+':
                seq_list.append(l)
                l = f.readline().rstrip()
            s = ''.join(seq_list)

            # check if the second sequence name (if existing) matches the first
            if len(l) > 1:
                if name != l[1:]:
                    raise SequenceInputError('Sequence and quality score ' +
                                             'identifiers do not match ' +
                                             '(\"%s\" != \"%s\")' % (name, l))
            l = f.readline().rstrip()

            # read in the quality scores
            # quality score lines can begin with an '@' character (ugh...), so
            # we can't just read until the next line that begins with an '@'
            q = ''
            while l != '' and len(q) < len(s):
                q += l
                l = f.readline().rstrip()

            # ensure that the quality and sequence strings are the same length
            if len(q) != len(s):
                raise SequenceInputError('Sequence and quality string ' +
                                         'lengths do not match (\"%s\")' % name)
            qual_list = []
            for c in q:
                qual_list.append(ord(c) - 33)

        yield Sequence(name, s, qual_list)

    f.close()


class Progress:
    """
    Simple progress indicator; indicates that long tasks are still running.
    """

    def __init__(self):
        self.num = 0
        self.display()

    def display(self):
        string = '%d records processed' % self.num
        sys.stdout.write(string)
        sys.stdout.flush()
        sys.stdout.write('\b' * len(string))

    def progress(self, num):
        self.num = num
        self.display()

    def erase(self):
        string = '%d records processed' % self.num
        sys.stdout.write(' ' * len(string))
        sys.stdout.write('\b' * len(string))

 
if __name__ == '__main__':
    seqs = 0
    seqlen = 0
    show_progress = False

    opts, args = getopt.getopt(sys.argv[1:], 'hp')

    for o, a in opts:
        if o == '-h':
            print 'Usage:'
            print '  %s [-p] <filename>' % sys.argv[0]
            print '  %s [-p] < <filename>' % sys.argv[0]
            print '  cat <filename> | %s [-p]\n' % sys.argv[0]
            print '  <filename> must be in fasta or fastq format, and can be'
            print '  gzipped, or (with the first usage style) bzip2\'d.'
            sys.exit(0)
        elif o == '-p':
            show_progress = True
        else:
            print 'unhandled option'
            sys.exit(1)

    if len(args) > 0:
        f = open(args[0], 'rb')
    else:
        f = sys.stdin

    if show_progress:
        p = Progress()

    try:
        for s in SequenceFile(f):
            #print s  # DEBUG
            seqs += 1
            seqlen += len(s.seq)
            if show_progress:
                p.progress(seqs)

        if show_progress:
            p.erase()
        print 'Total sequences found: %d' % seqs
        print 'Total residues found:  %d' % seqlen
    except SequenceInputError, e:
        print >> sys.stderr, '\nFaulty input (record %d): %s' % \
                             (seqs+1, e.args[0])

    f.close()

