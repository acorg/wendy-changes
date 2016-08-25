#!/usr/bin/env python

import sys
import argparse
from Bio import AlignIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError


def oldTranslate(seq):
    newseq = Seq(seq, IUPAC.unambiguous_dna)
    try:
        return newseq.translate()
    except KeyError:
        return None


def translate(seq):
    """
    Translate a nucleotide sequence to amino acids, converting
    any unknown codon (due to 'N' characters) into an 'X'.
    """
    aa = ''
    for start in range(0, len(seq), 3):
        newseq = Seq(seq[start:start + 3], IUPAC.unambiguous_dna)
        try:
            nextAA = newseq.translate()
        except TranslationError:
            nextAA = 'X'
        aa += nextAA
    return aa


parser = argparse.ArgumentParser(
    description=('Look for significant changes in an alignment. '
                 'Specify your reference sequence ids with '
                 '--reference and give an alignment on standard '
                 'input'))

parser.add_argument(
    '--reference', action='append', required=True,
    help='Give the name(s) of the reference sequence(s).')

args = parser.parse_args()

#print("The references given are:", args.reference)
#print("Total number of references given:", len(args.reference))

alignment = AlignIO.read(sys.stdin, 'fasta')
# print(alignment)
# print("Total number of sequences in alignment:", len(alignment))

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
# sets are different than list:
# no matter how many times we put 'a' in a set, it is just 'a': [a],
# whereas in a list, there are a number of a's: [a, a, a, ...n]
# print("References wanted are:", referencesWanted)

referenceSequences = []  # this makes it empty
otherSequences = []

print('\nsequences evaluated are:')

for record in alignment:
    #    print("record id:", record.id)
    #    print("record seq:", record.seq)
    if record.id in args.reference:
        print('reference', record.id)
        referenceSequences.append(str(record.seq))
        referencesWanted.remove(record.id)
        # I need to remove the id of my reference from this list, so
        # that I know all the references are in the list (like checklist)
#        print("referenceSequences are:", referenceSequences)
    else:
        otherSequences.append(str(record.seq))
        print('others', record.id)
#        print("other sequences:", otherSequences)

#in case the reference specified is not in the alignment
if referencesWanted:
    # if len(referencesWanted) != 0:
    print('The following reference sequences were not found:',
          ', '.join(sorted(referencesWanted)), file=sys.stderr)
    sys.exit(3)
"""
Nucleocapsid = alignment[:, 107:1679]
Phosphoprotein = alignment[:, 1800:3324]
Matrix = alignment[:, 3431:4439]
Fusion = alignment[:, 4934:6923]
Hemmaglutinin = alignment[:, 7078:8902]
LargeProt = alignment[:, 9029:15584]
"""
refseqNucleocapsid = []
refseqPhosphoprotein = []
refseqMatrix = []
refseqFusion = []
refseqHemmaglutinin = []
refseqLargeProtein = []

i = 0
for i in range(len(referenceSequences)):
    # record = alignment[i]
    record = referenceSequences[i]
    seq = str(record)

    # print(record.id, 'Nucleocapsid')

    # Nucleocapsid = str(seq[107:1679])
    # coding_Nucleocapsid = Seq(seq[107:1679], IUPAC.unambiguous_dna)
    # p_Nucleocapsid = coding_Nucleocapsid.translate()
    # print(p_Nucleocapsid)
    # sequencesNucleocapsid.append(p_Nucleocapsid)

    Nucleocapsid = str(seq[107:1679])
    # print('length of Nucleocapsid is', len(Nucleocapsid))
    p_Nucleocapsid = translate(Nucleocapsid)
    if p_Nucleocapsid is None:
        print('Ignoring sequence', record.id, file=sys.stderr)
        continue
    # print(p_Nucleocapsid)
    refseqNucleocapsid.append(p_Nucleocapsid)

    Phosphoprotein = str(seq[1800:3324])
    p_Phosphoprotein = translate(Phosphoprotein)
    if p_Phosphoprotein is None:
        print('Ignoring sequence', record.id, file=sys.stderr)
        continue
    refseqPhosphoprotein.append(p_Phosphoprotein)

    Matrix = str(seq[3431:4439])
    p_Matrix = translate(Matrix)
    if p_Matrix is None:
        print('Ignoring sequence', record.id, file=sys.stderr)
        continue
    refseqMatrix.append(p_Matrix)

    Fusion = str(seq[4934:6923])
    p_Fusion = translate(Fusion)
    if p_Fusion is None:
        print('Ignoring sequence', record.id, file=sys.stderr)
        continue
    refseqFusion.append(p_Fusion)

    Hemmaglutinin = str(seq[7078:8902])
    p_Hemmaglutinin = translate(Hemmaglutinin)
    if p_Hemmaglutinin is None:
        print('Ignoring sequence', record.id, file=sys.stderr)
        continue
    refseqHemmaglutinin.append(p_Hemmaglutinin)

    # LargeProtein = str(seq[9029:15584])
    # p_LargeProtein = translate(LargeProtein)
    # if p_LargeProtein is None:
    #     print('Ignoring sequence', record.id, file=sys.stderr)
    #     continue
    # refseqLargeProtein.append(p_LargeProtein)

    i = i + 1
"""
print ('\nreference aminoacid sequences:')
print ('Nucleocapsid\n', refseqNucleocapsid)
print ('Phosphoprotein\n', refseqPhosphoprotein)
print ('Matrix\n', refseqMatrix)
print ('Fusion\n', refseqFusion)
print ('Hemmaglutinin\n', refseqHemmaglutinin)
#print ('Large Protein\n', reseqLargeProtein)
"""
otseqNucleocapsid = []
otseqPhosphoprotein = []
otseqMatrix = []
otseqFusion = []
otseqHemmaglutinin = []
otseqLargeProtein = []

for i in range(len(otherSequences)):
#    record = alignment[i]
    record = otherSequences[i]
    seq = str(record)

#    print(record.id, 'Nucleocapsid')

#    Nucleocapsid = str(seq[107:1679])
#    coding_Nucleocapsid = Seq(seq[107:1679], IUPAC.unambiguous_dna)
#    p_Nucleocapsid = coding_Nucleocapsid.translate()
#    print(p_Nucleocapsid)
#    sequencesNucleocapsid.append(p_Nucleocapsid)

    Nucleocapsid = str(seq[107:1679])
#    print('length of Nucleocapsid is', len(Nucleocapsid))
    p_Nucleocapsid = translate(Nucleocapsid)
    if p_Nucleocapsid is None:
        print('Ignoring sequence', record.id, file=sys.stderr)
        continue
#    print(p_Nucleocapsid)
    otseqNucleocapsid.append(p_Nucleocapsid)

    Phosphoprotein = str(seq[1800:3324])
    p_Phosphoprotein = translate(Phosphoprotein)
    if p_Phosphoprotein is None:
        print('Ignoring sequence', record.id, file=sys.stderr)
        continue
    otseqPhosphoprotein.append(p_Phosphoprotein)

    Matrix = str(seq[3431:4439])
    p_Matrix = translate(Matrix)
    if p_Matrix is None:
        print('Ignoring sequence', record.id, file=sys.stderr)
        continue
    otseqMatrix.append(p_Matrix)

    Fusion = str(seq[4934:6923])
    p_Fusion = translate(Fusion)
    if p_Fusion is None:
        print('Ignoring sequence', record.id, file=sys.stderr)
        continue
    otseqFusion.append(p_Fusion)

    Hemmaglutinin = str(seq[7078:8902])
    p_Hemmaglutinin = translate(Hemmaglutinin)
    if p_Hemmaglutinin is None:
        print('Ignoring sequence', record.id, file=sys.stderr)
        continue
    otseqHemmaglutinin.append(p_Hemmaglutinin)

    LargeProtein = str(seq[9029:15584])
    p_LargeProtein = translate(LargeProtein)
    if p_LargeProtein is None:
        print('Ignoring sequence', record.id, file=sys.stderr)
        continue
    otseqLargeProtein.append(p_LargeProtein)

    i = i + 1
"""
print ('\nother aminoacid sequences:')
print ('Nucleocapsid\n', otseqNucleocapsid)
print ('Phosphoprotein\n', otseqPhosphoprotein)
print ('Matrix\n', otseqMatrix)
print ('Fusion\n', otseqFusion)
print ('Hemmaglutinin\n', otseqHemmaglutinin)
#print ('Large Protein\n', sequencesLargeProtein)
"""
#otseqProtein = (otseqNucleocapsid, otseqPhosphoprotein, otseqMatrix, otseqFusion, otseqHemmaglutinin)
#print ('\nlist of aa sequences in others:\n', otseqProtein)

N_aaLength = len(refseqNucleocapsid[0])
P_aaLength = len(refseqPhosphoprotein[0])
M_aaLength = len(refseqMatrix[0])
F_aaLength = len(refseqFusion[0])
H_aaLength = len(refseqHemmaglutinin[0])
#L_aaLength = len(refseqLargeProtein[0])
# [0] is the first reference sequence specified in the reference list

print('\nAA length of each protein:')
print("Nucleocapsid", N_aaLength)
print("Phosphoprotein:", P_aaLength)
print("Matrix:", M_aaLength)
print("Fusion:", F_aaLength)
print("Hemmaglutinin:", H_aaLength)

# Make a set of sequence positions for which all the other sequences have
# an identical aa.

N_identicalIndices = set()
for p in range(N_aaLength): #for each position in the aa length from reference
    aa = set()
    for sequence in otseqNucleocapsid:
        aa.add(sequence[p])
    #        print(aa)
    #        print(len(aa))
    # If all the aa at position p are the same,
    # remember this position.
    if len(aa) == 1:

# if it is not equal to 1, then there are more than 1 aa in each position
        N_identicalIndices.add(p)
#print (N_identicalIndices)

P_identicalIndices = set()
for p in range(P_aaLength):
    aa = set()
    for sequence in otseqPhosphoprotein:
        aa.add(sequence[p])
    if len(aa) == 1:
        P_identicalIndices.add(p)

M_identicalIndices = set()
for p in range(M_aaLength):
    aa = set()
    for sequence in otseqMatrix:
        aa.add(sequence[p])
    if len(aa) == 1:
        M_identicalIndices.add(p)

F_identicalIndices = set()
for p in range(F_aaLength):
    aa = set()
    for sequence in otseqFusion:
        aa.add(sequence[p])
    if len(aa) == 1:
        F_identicalIndices.add(p)

H_identicalIndices = set()
for p in range(H_aaLength):
    aa = set()
    for sequence in otseqHemmaglutinin:
        aa.add(sequence[p])
    if len(aa) == 1:
        H_identicalIndices.add(p)

    # For all the positions where the other sequences agree on the
    # aminoacid, print out the ones where none of the reference sequences
    # have that aminoacid.
print('\nAA substitutions: (America-1|POS|caspica.2000)')

for p in N_identicalIndices:
    N_refaminoacids = set()
    for sequence in refseqNucleocapsid:
        N_refaminoacids.add(sequence[p])
#        print(N_refaminoacids)
#    print(otseqNucleocapsid[0])
#    print(N_refaminoacids)
    if otseqNucleocapsid[0][p] not in N_refaminoacids:
        print('%s%d%s Nucleocapsid'
            % (otseqNucleocapsid[0][p], p + 1, ', '.join(N_refaminoacids)))

for p in P_identicalIndices:
    P_refaminoacids = set()
    for sequence in refseqPhosphoprotein:
        P_refaminoacids.add(sequence[p])
    if otseqPhosphoprotein[0][p] not in P_refaminoacids:
        print('%s%d%s Phosphoprotein'
            % (otseqPhosphoprotein[0][p], p + 1, ', '.join(P_refaminoacids)))

for p in M_identicalIndices:
    M_refaminoacids = set()
    for sequence in refseqMatrix:
        M_refaminoacids.add(sequence[p])
    if otseqMatrix[0][p] not in M_refaminoacids:
        print('%s%d%s Matrix'
            % (otseqMatrix[0][p], p + 1, ', '.join(M_refaminoacids)))

for p in F_identicalIndices:
    F_refaminoacids = set()
    for sequence in refseqFusion:
        F_refaminoacids.add(sequence[p])
    if otseqFusion[0][p] not in F_refaminoacids:
        print('%s%d%s Fusion'
            % (otseqFusion[0][p], p + 1, ', '.join(F_refaminoacids)))

for p in H_identicalIndices:
    H_refaminoacids = set()
    for sequence in refseqHemmaglutinin:
        H_refaminoacids.add(sequence[p])
    if otseqHemmaglutinin[0][p] not in H_refaminoacids:
        print('%s%d%s Hemmaglutinin'
            % (otseqHemmaglutinin[0][p], p + 1, ', '.join(H_refaminoacids)))
