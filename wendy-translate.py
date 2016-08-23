#!/usr/bin/env python

import sys
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

alignment = AlignIO.read(sys.stdin, 'fasta')

record = alignment[0]
seq = str(record.seq)

print()

substring = str(seq[20:41])
print('length of substring is', len(substring))
print('substring is %r' % substring)

section = Seq(seq[20:41], IUPAC.unambiguous_dna)
print('section', section)
print('translated', section.translate())
