#!/usr/bin/env python2.7
import argparse
from os import path
import sys
from helpers.exceptions import *
from helpers.IO import *
from helpers.alignments import filter_sequences
from helpers.workers import cast_variants

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
                    '[8]: varcode\nCan also be raw VEP output; in'\
                    'this case, must set -v to load it')
parser.add_argument('--expand', '-e', action='store_true', default=False,
                    help='expand frameshift (X follows position), '\
                    'early stops (- follows position) and deletions (-)')
parser.add_argument('--nopdb','-p', action='store_true',
                    help='do not cast to PDB structures')
parser.add_argument('--nomodel','-m', action='store_true',
                    help='do not cast to SWISSMODEL structures')
parser.add_argument('--nouniprot','-u',action='store_true',
                    help = 'do not cast to uniprot sequence. Only matters with models.')                    
parser.add_argument('--completed','-c',action='store_false',default=True,
                    help='suppress continuing later through completed file')
parser.add_argument('--num_procs','-n',type=int,default=1,
                    help='number of processes to run at once (default 1)')
parser.add_argument('--noalign','-a',action='store_true',
                    help='suppress creation of full alignment files')
parser.add_argument('--vep','-v',action='store_true',
                    help='input varfile is raw VEP output, select one '\
                         'consequence per coding variant and then cast')
parser.add_argument('--descriptors','-d',type=str,default="",
                    help='Positional descriptors, comma separated:\n'\
                         'all = all possible descriptors without artifacts;\n'\
                         'dssp = descriptors from DSSP (SS, SASA, isolated SASA);\n'\
                         'ds = distance to surface (0 if SASA>0.1);\n'\
                         'unp = uniprot annotations, see README for list;\n'\
                         'ligand = distance from ligand (-1 if no ligand);\n'\
                         'nucleotide = distance from dna/rna (-1 if neither);\n'\
                         'peptide = distance from binding peptide/protein;\n'\
                         "artifacts = don't filter ligand/complex artifacts;\n"\
                         " By default, ligands or biounits previously identified as"\
                         " artifacts in other databases are filtered.")                                           
parser.add_argument('--debug','-x',action='store_true',
                    help='Debug mode, prints all steps')

args = parser.parse_args()

assert not (args.nomodel and args.nopdb), \
    "You've selected neither PDB nor models, my job is done."

# Process the descriptors
if len(args.descriptors)>0:
    if args.debug:
        print "DEBUG: cast_variant: Processing descriptors"
    args.descriptors = args.descriptors.lower().split(",")
    if 'all' in args.descriptors:
        args.descriptors = ['dssp','ligand','nucleotide',
                            'peptide','unp'] #TODO: update as more descriptors implemented

# Load the required config settings
load_config(args)

# Check if necessary sequence files are there
check_seqs(args)

# Check for necessary descriptor applications/files
if len(args.descriptors)>0:
    check_applications(args)
    # If filter artifacts is set, load the artifacts file
    if "artifacts" not in args.descriptors:
        if args.debug:
            print "DEBUG: cast_variant: loading artifacts"
        try:
            load_artifacts()
        except ParseWarning as e:
            print e.fullmsg                        
else:
    if args.debug:
        print "DEBUG: cast_variant: descriptor list defined as empty"
    args.descriptors = list()    
# Define the names of the output files based on input filename
define_output(args)
# Read in the variant dict
variants = load_variants(args)
        
# Filter out completed variants if set
if args.completed:
    variants = filter_complete(variants,args.debug)

try:
    if len(variants.keys())==0:
        raise ParseException("variants_file","No usable variants, check for a skipped/completed file.")
except ParseException as e:
    sys.exit(e.fullmsg)        

if args.debug and False: #Probably unnecessary
    print "DEBUG: cast_variant: writing final variant list to file debugtmp.final_variants"
    with open("debugtmp.final_variants","w") as outfile:
        for t in variants:
            outfile.write(t+"\n")
            outfile.write("\n".join("{}".format(variants[t]) for t in variants))
            outfile.write("\n")
            
# Load the necessary sequence sets
# After loading each sequence set filter variants that
# don't have the necessary sequence
print "loading necessary datasets"
if args.debug:
    print "DEBUG: cast_variant: loading transcripts"
transcripts = load_transcripts()

datasets = {'transcripts':transcripts,
            'uniprots':None,
            'models':None,
            'sifts':None}

if args.debug:
    print "DEBUG: cast_variant: filtering variants for those without transcripts"
variants,earlycomp = filter_sequences(variants,'transcripts',datasets,args.debug)

try:
    if len(variants.keys())==0:
        raise ParseException("variants_file","No usable variants, check for a skipped file.")
except ParseException as e:
    sys.exit(e.fullmsg)  

if not (args.nouniprot and args.nopdb):
    datasets['uniprots'] = load_uniprot()
    variants,ec = filter_sequences(variants,'uniprots',datasets,args.debug)
    earlycomp += ec
try:
    if len(variants.keys())==0:
        raise ParseException("variants_file","No usable variants, check for a skipped file.")
except ParseException as e:
    sys.exit(e.fullmsg)  
    
if args.completed and len(earlycomp)>0:
    write_complete(earlycomp)

if not (args.nopdb):
    datasets['sifts'] = load_sifts()
if not (args.nomodel):
    datasets['models'] = load_models("swissmodel")
# Clear the previous files if continue from completed turned off
if not args.completed:
    remove_previous_files()
    
succ = cast_variants(variants,datasets,args)
print "Complete: {} transcripts successfully aligned".format(succ)


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
