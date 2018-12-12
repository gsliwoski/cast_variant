#!/usr/bin/env python2.7
import argparse
from os import path
import sys
import pickle

SEQ_PATH = "./sequences/"
SIFTS = "sifts.pickle"
SWISSMODEL = "swissmodel.pickle"
TRANS = "trans.pickle"
UNP = "unp.pickle"

parser = argparse.ArgumentParser(description = 
    'Quick check for presence of ENST or NP/XP in datasets and return fasta. '\
    'Also check for uniprot sequence given uniprot ID and return fasta. '\
    'Also do quick check for PDBs and/or models in dataset given a uniprot ID. '\
    'Note: dont provide isoform (-1, -2, etc) or version (.1, .2) to identifiers')

parser.add_argument('identifier', type=str,
                    help='Identifier to check. '\
                    'If starts with ENST: looks for transcript sequence. '\
                    'If starts with NP/XP: looks for protein sequence. '\
                    'If 6 characters: looks for uniprot sequence. '\
                    'Anything else do nothing.')
parser.add_argument('--pdb', '-p', action='store_true', default=False,
                    help='Check for PDB structures (only if uniprot ID provided)')
parser.add_argument('--model','-m', action='store_true',
                    help='Check for model structures (only if uniprot ID provided)')
parser.add_argument('--structure','-s', action='store_true',
                    help='Check for both PDB and model structures if uniprot ID provided')

args = parser.parse_args()

pdbs = dict()
models = dict()
transcripts = dict()
uniprot = dict()

def load_dataset(filename):
    assert path.isfile(filename), "Missing required file {}".format(filename)
    try:
        dataset = pickle.load(open(filename))
    except Exception as e:
        sys.exit("Failed to open required dataset with error: {}".format(e))
    return dataset

def check_structures(sources,identifier):
    if 'pdb' in sources:
        struct_hits = pdbs.get(identifier,set())
    if 'model' in sources:
        struct_hits = models.get(identifier,list())
    return struct_hits
    
def check_sequences(source,identifier):
    if source=='Uniprot':
        sequence = uniprot.get(identifier,"")
    else:
        sequence = transcripts.get(identifier,"")    
    return sequence
            
# Check the identifier type
if len(args.identifier)==6:
    print "loading uniprot dataset"
    uniprot = load_dataset(SEQ_PATH+UNP)
    idtype = "Uniprot"
elif args.identifier.startswith("ENST") or \
     args.identifier.startswith("XP") or \
     args.identifier.startswith("NP"):
        print "loading transcript dataset"
        transcripts = load_dataset(SEQ_PATH+TRANS)
        idtype = "Transcript"
else:
    sys.exit("Unrecognized identifier type")             

# Load only necessary datasets
if args.structure or args.pdb and idtype=="unp":
    print "loading PDB dataset"
    pdbs = load_dataset(SEQ_PATH+SIFTS)
if args.structure or args.model and idtype=="unp":
    print "loading model dataset"
    models = load_dataset(SEQ_PATH+SWISSMODEL)

# Check for sequence
fasta = check_sequences(idtype,args.identifier)
if len(fasta)==0:
    print "{} is not found in {} dataset".format(args.identifier,idtype)
else:
    print "{} found:".format(idtype)
    print "> {}\n{}".format(args.identifier,fasta)

if not args.pdb and not args.model and not args.structure:
    sys.exit()

# Check for structure
if idtype=="Uniprot" and len(fasta)>0:
    if args.structure or args.pdb:
        pdb_hits = check_structures("pdb",args.identifier)
        if len(pdb_hits)==0:
            print "No PDBs found for {}".format(args.identifier)
        else:
            print "{} pdbs found: {}".format(len(pdb_hits),", ".join(["_".join(x) for x in pdb_hits]))
    if args.structure or args.model:
        model_hits = check_structures("model",args.identifier)
        if len(model_hits)==0:
            print "No models found for {}".format(args.identifier)
        else:
            print "{} models found".format(len(model_hits))
            for m in sorted(model_hits,key = lambda(x): (x[0],-float(x[1]))):
                print "{} (identity = {}): {}".format(m[0],m[1],m[2])
