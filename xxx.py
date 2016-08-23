#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(
    description=('Look for significant changes in an alignment. '
                 'Specify your reference sequence ids with '
                 '--reference and give an alignment on standard '
                 'input'))

parser.add_argument(
    '--reference', action='append', required=True,
    help='Give the name(s) of the reference sequence(s).')

args = parser.parse_args()

print('references are', args.reference)
