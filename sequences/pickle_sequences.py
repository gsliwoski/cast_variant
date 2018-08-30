#!/usr/bin/env python2.7
from os.path import isfile
import gzip
import sys
import re
import pickle
import glob

# Gathers uniprot, ENST, and refseq fasta's into a dict and pickles it
# for faster loading.
# Required files in this directory (can be gz or not):
#   uniprot_sprot.fasta = all uniprot fasta's
#       ID format:
#       >sp|UNIPROT|Description (only keeps those containing OX=9606)
#
#   refseq_mammalian.fasta = all refseq fasta's
#       ID format:
#       >NP_#######.# description
#
#   ensembl_human.fasta = all ensembl fasta's
#       ID format:
#       >ENSP####.# pep loci gene:ENSG###.# transcript:ENST####.# desc
#
#   sifts file pdb_chain_uniprot.tsv
#   also need a path to the xml folder of sifts data
#   see /capra_lab/data_clean/sifts/ for the latest files
#
#   swissmodel_index.tab = INDEX from swissmodel database

# If a pickle is already present, it will skip
# Unless -f is passed, in which case it will create it again

# Define the required files/paths                
datapath = "/dors/capra_lab/data_clean/"
unpfile = "uniprot_sprot.fasta.gz"
enstfile = "Homo_sapiens.GRCh38.pep.all.fa.gz"
refseqfile = datapath + "refseq/2018-07-24/refseq_mammalian.fasta.gz"
siftsfile = "pdb_chain_uniprot.tsv.gz"
swissmodelpath = datapath + "swissmodel/2018-07-23/SWISS-MODEL_Repository/"
swissmodelindex = swissmodelpath + "INDEX"
# Define output files
unpout = "unp.pickle"
transout = "trans.pickle"
swissout = "swissmodel.pickle"
siftsout = "sifts.pickle"

try:
    force = sys.argv[1]
except:
    force = ""
overwrite = force=="-f"
        
def get_id(line,source):
    '''Gets the id of the sequence depending on source.'''
    if source=="ensembl":
        line = line.strip().split()
        seqid = [x for x in line if x.startswith("transcript:")][0]
        seqid = re.split("[:.]",seqid)[1]
        assert seqid.startswith("ENST"),\
            "unrecognized enst id {}".format(seqid)
    elif source=="refseq":
        seqid = line.strip().split()[0].lstrip(">")
    elif source=="uniprot":
        seqid = line.strip().split("|")[1]
    else:
        sys.exit("unrecognized source of fasta: {}".format(source))
    return seqid

def check_species(line,source):
    '''checks if species is human. Only relevant to uniprot
       or if refseq, checks if it's a NP or XP sequence'''
    if source=="uniprot":
        return "OX=9606" in line
    elif source=="refseq":
        return line.startswith("XP") or line.startswith("NP")
    
def parse_fasta(fastafile,target,source):
    '''Parses a standard fasta file.
       Into dict of key = identifier, entry = fasta sequence
       Optionally takes an id_filter that stores only containing ids
       optional argument checkid is bool to check id (uniprot)'''
    current_id = None
    if source=="uniprot":
        passed = False
    else:
        passed = True        
    current_fasta = ""
    with gzip.open(fastafile) as infile:
        for line in infile:
            if line.strip()=="" or line.strip()=="*": continue
            if line.strip().startswith(">"):
                if current_id in target:
                    print "warning: {} found multiple times".format(current_id)
                if current_fasta != "" and passed:
                    if source=="refseq":
                        passed = check_species(current_id,source)
                    if passed:               
                        target[current_id] = current_fasta.rstrip("*")
                current_id = get_id(line,source)
                if source=="uniprot":
                    passed = check_species(line,source)
                current_fasta = ""
            else:
                current_fasta += line.strip()
        if source=="refseq":
            passed = check_species(current_id,source)
        if current_fasta!="" and passed:
            target[current_id] = current_fasta.rstrip("*")
        
def parse_special(fastafile,target):
    '''Parses the tab delim fasta format into
       dict with key = identifier (col1) and
       entry = fasta seq (col2)'''
    with gzip.open(fastafile) as infile:
        for line in infile:
            line = line.strip().split()
            if len(line)<2: continue
            if line[0] in transcripts:
                print "Warning {} seen multiple times.".format(line[0])
            else:
                target[line[0]] = line[1]                

def parse_swissmodel(target,index,path):
    '''Make a dictionary of all available swissmodels
       and their basic info
       dict has entries:
       key = uniprot
       entry = list of models with format:
        [isoform,identity,filename.pdb]
       '''
    models = dict()
    with open(index) as infile:
            for line in infile:
                if line.startswith("#") or line.startswith("UniProtKB_ac"):
                    continue
                line = line.strip().split("\t")
                if len(line)<8: continue
                unp = line[0].split("-")[0]
                start,end = line[5:7]
                key = line[3]
                acid,isoid = line[:2]
                if isoid.strip()=="":
                    isoid = "UNK"
                if unp in models:
                    models[unp].append([acid,isoid,key,start,end])
                else:
                    models[unp] = [[acid,isoid,key,start,end]]
    # Get the necessary information for all the models
    i = 0
    t = len(models.keys())
    for unp in models.keys():
        if i>0 and i % 500 == 0:
            sys.stdout.write("\r{} of {} uniprots parsed".format(i,t))
            sys.stdout.flush()
        i += 1                    
        current = models[unp]
        srcpath = (unp[:2],unp[2:4],unp[4:])
        for x in current:
            acid,isoid,key,start,end = x
            filename = glob.glob("{}{}/{}/{}/swissmodel/*{}*".format(
                                            swissmodelpath,
                                            srcpath[0],
                                            srcpath[1],
                                            srcpath[2],
                                            key))
            if len(filename)==0: continue
            elif len(filename)>1:
                print "Error key {} in unp {} matched multiple models".format(key,unp)
                continue
            seqid = None
            with open(filename[0]) as infile:
                for line in infile:
                    if line.strip().startswith("REMARK"):
                        line = line.strip().split()
                        if len(line)<4: continue
                        if line[2]=="SID":
                            seqid = line[3]
                            break
            if seqid is None:
                print "Warning: no identity for unp {} file {}".format(
                                                            unp,filename[0])
                seqid="0.00"
            if unp in target:
                target[unp].append([isoid,seqid,filename[0]])
            else:
                target[unp] = [[isoid,seqid,filename[0]]]
    print ""    
def parse_sifts(siftsfile,target):
    '''returns dict with entries
       key = uniprot ID
       entry = [pdbid,chain]'''
    with gzip.open(siftsfile) as infile:
        for line in infile:
            if line.startswith("#") or line.startswith("PDB"):
                continue
            line = line.strip().split()
            if len(line)<3: continue
            pdbid,chain,unp = line[:3]
            if unp in target:
                target[unp].add((pdbid,chain))
            else:
                target[unp] = set([(pdbid,chain)])

def check_skip(destination):
    '''Checks if pickle is already there.
       If it is, skips unless -f flag used'''
    if isfile(destination) and not overwrite:
        print "{} pickle already found, skipping".format(destination)
        print "to prevent skipping, run with -f"
        return True
    else:
        return False

# Make sure files are there
for x in [unpfile,enstfile,refseqfile,siftsfile,swissmodelindex]:
    assert isfile(x), "{} is missing.".format(x)

if not check_skip(unpout):
    print "pickling uniprot seqs"
    uniprot = dict()
    parse_fasta(unpfile,uniprot,"uniprot")
    if len(uniprot) > 0:
        pickle.dump(uniprot,open("unp.pickle",'w'))
    print "{} uniprots pickled.".format(len(uniprot.keys()))

if not check_skip(transout):
    print "pickling transcript seqs"
    transcripts = dict()
    print "RefSeq"
    parse_fasta(refseqfile,transcripts,"refseq")
    print "Ensembl"
    parse_fasta(enstfile,transcripts,"ensembl")
    if len(transcripts)>0:
        pickle.dump(transcripts,open("trans.pickle",'w'))
    print "{} transcripts pickled.".format(len(transcripts.keys()))

if not check_skip(swissout):
    print "pickling swissmodel table"
    modelmap = dict()
    parse_swissmodel(modelmap,swissmodelindex,swissmodelpath)
    if len(modelmap)>0:
        pickle.dump(modelmap,open("swissmodel.pickle",'w'))
    print "{} unp modelsets pickled".format(len(modelmap.keys()))

if not check_skip(siftsout):
    sift = dict()
    print "pickling sifts table"
    parse_sifts(siftsfile,sift)
    if len(sift)>0:
        pickle.dump(sift,open("sifts.pickle",'w'))
    print "{} unp sifts sets pickled".format(len(sift.keys()))

print "Done pickling"
    
