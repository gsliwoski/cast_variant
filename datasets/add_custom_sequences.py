#!/usr/bin/env python2.7
# This is meant to be run after setting up uniprot and transcript fasta's to allow users' to provide
# Custom transcripts

import sys
from os.path import isfile
import pickle

try:
    fastafile = sys.argv[1]
    destfile = sys.argv[2]
except:
    sys.exit("add_custom_sequences.py name_of_file_containing_custom_fastas.fasta destination_pickle(unp.pickle or trans.pickle)")

try:
    skipper = sys.argv[3]
except:
    skipper = 0
        
assert isfile(fastafile), "{} not found".format(fastafile)
assert isfile(destfile), "{} not found, make sure to run pickle_sequences.py first".format(destfile)

new = dict()
current_id = None
current_fasta = ""
with open(fastafile) as infile:
    for line in infile:
        if line.strip()=="" or line.strip()=="*": continue
        if line.strip().startswith(">"):
            if current_id in new:
                print "warning: {} found multiple times".format(current_id)
            if current_fasta != "":
                new[current_id] = current_fasta.rstrip("*")
            current_id = line.lstrip(">").strip()
            current_fasta = ""
        else:
            current_fasta += line.strip()
    if current_fasta!="":
        new[current_id] = current_fasta.rstrip("*")            

assert len(new.keys()) > 0, "No new fasta files loaded!"
print "loading"     
destdict = pickle.load(open(destfile))
print "adding"
n = 0
for newfasta in new:
    if newfasta in destdict:
        if skipper == 0:
            print "Warning, you overwrote a fasta already in the dict ({})".format(newfasta)
            print "To skip anything already in destination, add a 1 to the end of the command"
            destdict[newfasta] = new[newfasta]
            n += 1
        else:
            print "Skipping {} already in destination".format(newfasta)
    else:
        destdict[newfasta] = new[newfasta]
        n += 1
print "saving"
pickle.dump(destdict,open(destfile,'w'))

print "Added {} new sequence id's to {}".format(n,destfile)
