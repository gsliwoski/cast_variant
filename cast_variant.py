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
parser.add_argument('--expand', '-e', action='store_true', default=True,
                    help='expand frameshift (X follows position), '\
                    'early stops (- follows position) and deletions (-)')
parser.add_argument('--nopdb','-p', action='store_true',
                    help='do not cast to PDB structures')
parser.add_argument('--nomodel','-m', action='store_true',
                    help='do not cast to SWISSMODEL structures')
parser.add_argument('--completed','-c',action='store_true',default=True,
                    help='allow for continuing later with completed file')                                        

args = parser.parse_args()

assert not (args.nomodel and args.nopdb), \
    "You've selected not to cast to PDB nor models, my job is done."

# Check if necessary sequence files are there
check_seqs(args.nomodel,args.nopdb)
# Define the names of the output files based on input filename
define_output(args.variants)
# Read in the variant dict
variants = process_variants(args.variants,args.expand,args.completed)
try:
    if len(variants.keys())==0:
        raise ParseException("variants_file","No usable variants, check for a skipped file.")
except ParseException as e:
    sys.exit(e.fullmsg)        
#print variants

# Load the necessary sequence sets
print "loading necessary datasets"
transcripts = load_transcripts()
if not args.nopdb:
    uniprot = load_uniprot()
    sifts = load_sifts()

# Load the desired model tables if necessary
if not args.nomodel:
    models = load_models("swissmodel")
#x = 0
#for a,b in zip(transcripts,uniprot):
#    if x==10: break
#    print a,transcripts[a]
#    print b,uniprot[b]
#    x += 1
#x = 0
#for a,b in zip(sifts,models):
#    if x==10: break
#    print a,sifts[a]
#    print b,models[b]
#    x += 1

    
# Cast to PDBs
#if not args.nopdb:
#    cast_to_
