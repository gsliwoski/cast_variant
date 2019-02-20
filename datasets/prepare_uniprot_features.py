#!/usr/bin/env python2.7
#import xml.etree.ElementTree as et
import copy
from lxml import etree
import psutil
import os,re,sys
#from pandas import DataFrame
import pandas as pd
import pickle
pid = os.getpid()

# Annotations are grouped and features are True/False depending on
# whether one or more features in the group is found at that residue
# Groups are primarily based on uniprot groupings:
# https://www.uniprot.org/help/sequence_annotation
# But some features are currently ignored and appear in the 'unused' set
# 
# Any feature with unknown status at start or end is ignored

# Defintions:
#  top: Topological domain (M = Transmembrane, E = extracellular,
#                           C = Cytoplasmic, L=lumenal, LM = "Lumenal, melanosome", I = intramembranem,
#                           MI = Mitochondrial Intermembrane, MM = mitochondrial matrix,
#                           N = Nuclear, NP = Perinuclear space, V = Vesicular, VL = "lumenal, vesicle", 
#                           PM = Peroxisomal matrix, P = Peroxisomal, X = Exoplasmic loop,
#                           O = other) [only stores first occurrence] (region)

# Change location of config file if using unexpected hierarchy
#CONFIG_FILE = "../config.sys"

feature_groups = {'cleaved': ['initiator methionine',
                             'propeptide'],
                  'sorting': ['signal',
                              'transit peptide',
                              'signal peptide'],
                  'membrane': ['transmembrane',
                               'intramembrane',
                               'transmembrane region',
                               'intramembrane region'],
                  'binding region': ['calcium binding',
                                     'DNA binding',
                                     'DNA-binding region',
                                     'nucleotide binding',
                                     'calcium-binding region',
                                     'phosphate-binding region',
                                     'nucleotide phosphate-binding region',
                                     'lipid moiety-binding region'],
                  'motif': ['motif',
                            'short sequence motif'],
                  'binding site': ['active site',
                                   'metal binding',
                                   'binding site',
                                   'site',
                                   'metal ion-binding site'],
                  'ptm': ['modified residue',
                          'lipidation',
                          'glycosylation',
                          'glycosylation site'],
                  'covalent': ['disulfide bond',
                               'cross-link'],
                  'unused': ['chain',
                             'peptide',
                             'domain',
                             'repeat',
                             'coiled coil',
                             'coiled-coil region',
                             'compositional bias',
                             'non-standard residue',
                             'alternative sequence',
                             'natural variant',
                             'mutagenesis',
                             'mutagenesis site',
                             'sequence uncertainty',
                             'sequence conflict',
                             'non-adjacent residues',
                             'non-terminal residues',
                             'helix',
                             'turn',
                             'beta strand',
                             'strand',
                             'sequence variant',
                             'splice variant',
                             'compositionally biased region',
                             'zinc finger region',
                             'non-standard amino acid',
                             'non-terminal residue',
                             'topological domain',
                             'region',
                             'region of interest',
                             'non-consecutive residues']}

#assert os.path.isfile(CONFIG_FILE), "{} not found".format(CONFIG_FILE)

DATA_FILE = "./uniprot_sprot.xml"
OUT_FILE = "unp_features.pickle"

#with open(CONFIG_FILE) as infile:
#    for rawline in infile:
#        line = rawline.strip().split("#")[0].split()
#        if len(line)<2:continue
#        if line[0]=="UNP_RAW":
#            DATA_FILE = line[1]
#        elif line[0]=="UNP_FEATURES":
#            OUT_FILE = line[1]                        
#        if DATA_FILE is not None and OUT_FILE is not None:
#            break

#assert DATA_FILE is not None, "Missing UNP_RAW line in {}".format(CONFIG_FILE)                             
#assert OUT_FILE is not None, "Missing UNP_FEATURES line in {}".format(CONFIG_FILE)

WARNINGS_FILE = OUT_FILE+".warnings"

empty_res = {'cleaved': False,
             'sorting': False,
             'membrane': False,
             'binding region': False,
             'motif': False,
             'binding site': False,
             'ptm': False,
             'covalent': False}

all_residues = {'uniprot':list(),
                'uniprot_position':list(),
                'cleaved':list(),
                'sorting':list(),
                'membrane':list(),
                'binding region':list(),
                'motif':list(),
                'binding site':list(),
                'ptm':list(),
                'covalent':list()}
  
def append_entries(df,d):
    ddf = pd.DataFrame.from_dict(d)
    return pd.concat([df,ddf])
                          
def translate_feature(fnode,annotations):
    global fucounts,fcounts
    warnings = ""
    ftype = fnode.get('type')
    fdesc = fnode.get('description')
    if ftype is None:
#        print "Warning %s had a none-type feature" % annotations['name']
        warnings+="Warning {} had a none-type feature\n".format(annotations['uniprot'])
        return warnings
    if fdesc is None:
        fdesc = "NoDesc"
    begin = xbegin(fnode)
    end = xend(fnode)
    singlepos = xpos(fnode)
    positions = set()
    #Break out all the sequence positions
    if len(begin)>0 and len(end)>0:
        #Disulfide bonds aren't ranges but partners
        if fnode.get('type')=='disulfide bond':
            positions.add(int(begin[0]))
            positions.add(int(end[0]))
        else:
            positions.update(range(int(begin[0]),int(end[0])+1))
    if len(singlepos)>0:
        positions.add(int(singlepos[0]))
    if len(positions)==0:
        #Don't bother warning if begin/end annotated as 'unknown'
        us = xus(fnode)
        if 'unknown' in us:
            return warnings
#        print "Warning %s type %s has no positions" % (annotations['name'],ftype)
        warnings+="Warning {} type {} has no positions\n".format(annotations['uniprot'],ftype)
        return warnings
    #initialize necessary positions
    for x in positions:
        if x not in annotations:
            annotations[x] = copy.deepcopy(empty_res)
    Known = False
    #Process topology domain   
    if ftype=='topological domain':
        Known = True
        if fdesc == 'Extracellular':
            curtop = 'E'
        elif fdesc == 'Cytoplasmic':
            curtop = 'C'
        elif fdesc=='Lumenal':
            curtop = 'L'
        elif fdesc=='Lumenal, melanosome':
            curtop = 'LM'
        elif fdesc == 'Mitochondrial intermembrane':
            curtop = 'MI'
        elif fdesc == 'Mitochondrial matrix':
            curtop = 'MM'
        elif fdesc == 'Perinuclear space':
            curtop = 'NP'
        elif fdesc == 'Nuclear':
            curtop = 'N'
        elif fdesc == 'Vesicular':
            curtop = 'V'
        elif fdesc == 'Lumenal, vesicle':
            curtop = 'VL'
        elif fdesc == 'Peroxisomal matrix':
            curtop = 'PM'
        elif fdesc == 'Peroxisomal':
            curtop = 'P'
        elif fdesc == 'Exoplasmic loop':
            curtop = 'X'
        else:
            #print "Warning: %s has unrecognized top description" % annotations['name']
            warnings += "Warning: {} has unrecognized topology description\n".format(annotations['uniprot'])
            curtop='?'
        # Not going to use this and restrict descriptors to True/False but could add in future
#        for x in positions:
#            if annotations[x]['top']==empty_res['top']:
#                annotations[x]['top']=curtop
    
    # Process feature groups
    for f in feature_groups:
        if ftype in feature_groups[f]:
            Known = True
            fucounts[f].add(annotations['uniprot'])
            for x in positions:
                fcounts[f] += 1
                annotations[x][f] = True
            break
    #Unrecognized feature                            
    if not Known:
        print "Warning: {} has unrecognized entry type {}\n".format(annotations['uniprot'],ftype)            
        warnings+="Warning: {} has unrecognized entry type {}\n".format(annotations['uniprot'],ftype)        
        return warnings

    return warnings
    
def assign_aa(snode,annotations):
    if len(snode)==0:
        return False
    seq = re.sub(r'\s+','',snode[0].text.strip())
    if len(seq)==0:
        return False
    for i,c in enumerate(seq):
        pos = i+1
        if pos in annotations:
            annotations[pos]['aa'] = c
    return True

# Iterative tree to avoid loading large XML into memory
context = etree.iterparse(DATA_FILE,events=('start','end'))
context = iter(context)
event,root = context.next()

#Define xpath queries
xname = etree.XPath("./*[local-name()='name']/text()")
xacc = etree.XPath("./*[local-name()='accession']/text()")
h0 = etree.XPath("./*[local-name()='dbReference']/@id")
h1 = etree.XPath("./*[local-name()='organism']/*[local-name()='dbReference']/@id")
h2 = etree.XPath("./*[local-name()='organism']/*[local-name()='name']/text()")
xfeatstring = "./*[local-name()='feature'"
for f in feature_groups['unused']:
    xfeatstring += " and @type!='{}'".format(f)
xfeatstring += "]"
xfeat = etree.XPath(xfeatstring)
xbegin = etree.XPath("./*[local-name()='location']/*[local-name()='begin']/@position")
xend = etree.XPath("./*[local-name()='location']/*[local-name()='end']/@position")
xpos = etree.XPath("./*[local-name()='location']/*[local-name()='position']/@position")
xus =  etree.XPath("./*[local-name()='location']/*[local-name()='begin' or local-name()='end']/@status")
                             
# Initialize empty dataframes and dicts
final_df = None
n = 0
annotations_set = dict()
for x in all_residues:
    annotations_set[x] = list()
warnings_set = ""
nn = 1
print "Extracting feature groups from all human residues in {}".format(DATA_FILE)

# To see what features may be too over/under-represented, keep track of counts
fcounts = dict()
fucounts = dict()
for k in empty_res.keys():
    fcounts[k] = 0
    fucounts[k] = set()
    
for event,elem in context:
    # Append every 500 uniprots
    if n>=(500*nn):
        print "{} uniprots complete".format(n)
        if warnings_set!="":
            with open(WARNINGS_FILE,'a') as outfile:
                outfile.write(warnings_set)
            warnings_set = ""           
        if final_df is None:
            final_df = pd.DataFrame.from_dict(annotations_set)
        else:
            final_df = append_entries(final_df,annotations_set)
#        final_df.to_csv("testunp.tab",columns=['uniprot','uniprot_position','cleaved','sorting','membrane','motif','binding region','binding site','ptm','covalent'],sep="\t",index=False)
        annotations_set = dict()
        for x in all_residues:
            annotations_set[x] = list()
        nn += 1
    if event=='end' and elem.xpath("local-name()")=="entry":
        accession = ",".join(xacc(elem))
        unp = xacc(elem)[0]
        name = ",".join(xname(elem))
        # Only use human uniprots
        human = name.find("HUMAN")>0\
            or "9606" in h0(elem)\
            or len(set(h1(elem) + h2(elem)).intersection(set(["Human","Homo sapiens"]))) > 0
        if human:
            annotations = {'uniprot':unp}
            feature_nodes = xfeat(elem)
            current_warnings = ""
            for fnode in feature_nodes:
                current_warnings += translate_feature(fnode,annotations)
#                assign_aa(elem.xpath("./*[local-name()='sequence']"),annotations)
#                warnings_set += current_warnings
            for r in annotations:
                if r == 'uniprot': continue
                for f in annotations[r]:
                    annotations_set[f].append(annotations[r][f])
                annotations_set['uniprot'].append(annotations['uniprot'])
                annotations_set['uniprot_position'].append(r)                        
            n += 1
        root.clear()
        elem.clear()
        while elem.getprevious() is not None:
            del elem.getparent()[0]

if final_df is None:
    final_df = pd.DataFrame.from_dict(annotations_set)
else:
    final_df = append_entries(final_df,annotations_set)

print "Finished compiling uniprot annotations across {} uniprots and {} residues".format(n,len(final_df.index))
print "Total residues flagged for each feature:"
for k in fcounts:
    print "{}: {} uniprots, {} residues".format(k,len(fucounts[k]),fcounts[k])
print "Pickled to {}".format(OUT_FILE)
final_df.to_pickle(OUT_FILE)                 
