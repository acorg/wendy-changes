#!/usr/bin/env python

import sys
import argparse
from Bio import AlignIO


def main():
    parser = argparse.ArgumentParser(
        description=('Look for significant changes in an alignment. '
                     'Specify your reference sequence ids with '
                     '--reference and give an alignment on standard '
                     'input'))

    parser.add_argument(
        '--reference', action='append', required=True,
        help='Give the name(s) of the reference sequence(s).')

    parser.add_argument(
        '--region', action='append', required=True,
        help='Give the range and name of the wanted translated region(s).')

    args = parser.parse_args()

    alignment = AlignIO.read(sys.stdin, 'fasta')

    if len(args.reference) > len(alignment):
        print('You specified %d reference sequence ids, but your '
              'alignment only has %d sequences.' % (
                  len(args.reference), len(alignment)),
              file=sys.stderr)
        sys.exit(1)
    elif len(args.reference) == len(alignment):
        print('You specified %d reference sequence ids, but you '
              'only have that many sequences.' % (
                  len(args.reference)),
              file=sys.stderr)
        sys.exit(2)

    referencesWanted = set(args.reference)
    referenceSequences = []
    otherSequences = []

    for record in alignment:
        if record.id in args.reference:
            referenceSequences.append(str(record.seq))
            referencesWanted.remove(record.id)
        else:
            otherSequences.append(str(record.seq))

    if referencesWanted:
        print('The following reference sequences were not found:',
              ', '.join(sorted(referencesWanted)), file=sys.stderr)
        sys.exit(3)

    seqLength = len(referenceSequences[0])

    # Make a set of sequence indices for which all the other sequences have
    # an identical nucleotide.
    identicalIndices = set()

    for i in range(seqLength):
        nucleotides = set()
        for sequence in otherSequences:
            nucleotides.add(sequence[i])

        # If all the nucleotides at index i are the same, remember this
        # index.
        if len(nucleotides) == 1:
            identicalIndices.add(i)

    # For all the indices where the other sequences agree on the
    # nucleotide, print out the ones where none of the reference sequences
    # have that nucleotide.
    for i in identicalIndices:
        referenceNucleotides = set()
        for sequence in referenceSequences:
            referenceNucleotides.add(sequence[i])

        # We can safely use otherSequences[0][i] because 1) The
        # otherSequences list cannot be empty and 2) because all the other
        # sequences agree at position i.
        if otherSequences[0][i] not in referenceNucleotides:
            print('Nucleotide %s at position %d occurs in all other '
                  'sequences but not in any reference sequence. The '
                  'reference sequences have nucleotides %s.' %
                  (otherSequences[0][i], i + 1,
                   ', '.join(referenceNucleotides)))

if __name__ == '__main__':
    main()
