#!/usr/bin/env python2.7
import gzip,re,sys

unp_file = "uniprot_sprot.dat.gz"
unp_outfile = "uniprot_sprot_human.dat"
canonical_outfile = "uniprot_canonical_isoforms.tab"

# Known errors
# P61812 - doesn't conform to standard displayed= section
human = list()
# Step 1: filter for human proteins only
print "Isolating human"
with gzip.open(unp_file,'r') as fin:
    data_buffer = []
    confirmed_human=False
    for line in fin:
        row = line.strip().split('   ')
        if row[0] == '//':
            # Found new ID, if the previous entry was human,
            # flush the buffer
            if confirmed_human:
                human.extend(data_buffer)
            # Wait for confirmation that next entry is human
            confirmed_human = False
            # Clear the data buffer for the next entry
            data_buffer = []
        elif row[0] == 'OS' and row[1] == 'Homo sapiens (Human).':
            # The current entry is human, flush when finished
            confirmed_human = True
        # Store the row in the data buffer in case it is
        # human and needs to be printed
        data_buffer.append(line)

with open(unp_outfile,"w") as outfile:
    outfile.write("".join(human))

# Step 2: create canonical isoform table
print "Creating canonical isoform table"		
output = open(canonical_outfile,"w")
unps = list()
iso = dict()
trans = dict()
p = re.compile("\[[A-Z0-9]*-[A-Z0-9]*\]")

output.write("\t".join(["Uniprot","Transcript","Protein","Isoform","Canonical","nIsoforms"]) + "\n")

def addunps(u,i,t):
    global output
    i["Unassigned"]="Unassigned"
    for x in t:
        if t[x] not in i:
            i[t[x]] = "NO"
    numisos = len(i.keys())-1
    if len(t)==0 and len(u)>0:
        if len(i)>1:
            pass
            #print "Warning, {} had no transcripts but had isoforms".format(u)
        else:
            pass
            #print "Warning, {} had no transcripts".format(u)
        for x in u:
            output.write("{}\tNA\tNA\tNA\tNA\tnumisos\n".format(x))                          
    else:
        for x in u:
            for y in t:
                output.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(x,y[0],y[1],t[y],i[t[y]],numisos))           
    return dict(),dict(),list()

ln,pln = 0,0
for line in human:
    ln += 1
    line = line.strip().split()
    if len(line)==1: continue
    if line[0] == "AC":
        if ln>pln+1:
            iso,trans,unps = addunps(unps,iso,trans)
        unps+=[x.rstrip(";") for x in line[1:]]
        pln = ln
    elif line[0] == "CC" and line[1].startswith("IsoId"):
        i = line[1].split("=")[1].rstrip(";")
        if "=" in line[2]:
            s = "YES" if line[2].split("=")[1].rstrip(";")=="Displayed" else "NO"
            if i in iso:
                print "Warning {} has a repeated isoform".format(unps)
            iso[i] = s     
    elif line[0]=="DR":
        if line[1]=="RefSeq;" or line[1]=="Ensembl;":
            if p.match(line[-1])>-1:
                assignment = line[-1].strip("[]")
            else:
                assignment = "Unassigned"                               
            if line[1]=="RefSeq;":               
                ct = (line[3][:-1],line[2][:-1])
            else:
                ct = (line[2][:-1],line[3][:-1])                
            if ct in trans:
                print "Warning {} from {} appears multiple times".format(ct,unps)
            trans[ct]=assignment
output.close()
print "Finished processing, the following files can be removed:\n"
print "uniprot_sprot.dat.gz"
print "uniprot_sprot_human.dat"
