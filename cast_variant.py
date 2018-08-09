#!/usr/bin/env python2.7
import argparse
from os import path
import sys
from helpers.exceptions import *
from helpers.IO import *

# Define the sequence files as those previously pickled
# If you changed the filenames output from 
parser = argparse.ArgumentParser(description = 
    'Cast a list of variants onto protein structures. '\
    'Variants must have been run through an effect predictor first. '\
    'For best results, use VEP38. Only missense variants are auto-run. '\
    'However, you can \"expand\" frameshift, early stop, and deletions.')

parser.add_argument('variants', type=str,
                    help='List of variants. Each line is tab separated. '\
                    'Required columns:\n1: consequence\n'\
                    '2: Transcript ID\n3: Transcript position\n'\
                    '4: Mutation as refaa/altaa\n5: Protein ID\n'\
                    '6: Uniprot ID\n7: Uniprot-isoform\n'\
                    '[8]: varcode')
parser.add_argument('--expand', '-e', action='store_true',
                    help='expand frameshift (X follows position), '\
                    'early stops (- follows position) and deletions (-)')
parser.add_argument('--nopdb','-p', action='store_true',
                    help='do not cast to PDB structures')
parser.add_argument('--nomodel','-m', action='store_true',
                    help='do not cast to SWISSMODEL structures')

args = parser.parse_args()

assert not (args.nomodel and args.nopdb), \
    "You've selected not to cast to PDB nor models, my job is done."

# Check if necessary sequence files are there
check_seqs(args.nomodel,args.nopdb)

define_output(args.variants)
variants = process_variants(args.variants,args.expand)
print variants
# Cast to PDBs
#if not args.nopdb:
#    cast_to_
