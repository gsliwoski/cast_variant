import multiprocessing as mp
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from exceptions import *
import IO
import traceback
import pandas as pd
import copy

# Alignment parameters
GAP_OPEN = -10
GAP_EXTEND = -0.5
MATRIX = matlist.blosum62

# Main alignment class
class Alignment(object):
    '''
    Main alignment object starts with 2 sequences and aligns them
    Then can populate either side of alignment with different sequences
    such as numbering and icodes that are adjusted to account for gaps
    Can also get identity and final dataframe
    '''
    def __init__(self, seq1, seq2, debug, name1='raw1', name2='raw2'):
        self.identical = seq1==seq2
        self.left_name = name1
        self.right_name = name2
        self.left_sequence = dict()
        self.right_sequence = dict()
        self.matches = 0
        self._df = None
        self._dict = None
        self.gap_open = GAP_OPEN
        self.gap_extend = GAP_EXTEND
        self.matrix = MATRIX
        self.current = False # df up-to-date?
        self.debug = debug
        self.debug_head = "DEBUG: alignments: Alignment: "
        if self.left_name == self.right_name:
            raise AlignException("pairwise aln failed",
                                 "provided identifiers can't be identical")
       
        self.left_sequence['nAA'] = len([x for x in seq1 if x!="-"])
        self.right_sequence['nAA'] = len([x for x in seq2 if x!="-"])
        seq1 = seq1.replace("-","X")
        seq2 = seq2.replace("-","X")
        if self.debug:
            print self.debug_head+"generated align class with {} ({} AA) and {} ({} AA)".format(
                        self.left_name, 
                        self.left_sequence['nAA'], 
                        self.right_name, 
                        self.right_sequence['nAA'])
        try:
            if self.identical:
                self.left_sequence['aligned'] = [x for x in seq1]
                self.right_sequence['aligned'] = [x for x in seq2]
                self.matches = len(seq1)
                if self.debug:
                    print self.debug_head+ "Sequences are identical"
            else:
                if self.debug:
                    print self.debug_head+"doing pairwise alignment"
                l,r = pairwise2.align.globalds(
                                             seq1,
                                             seq2,
                                             MATRIX,
                                             GAP_OPEN,
                                             GAP_EXTEND,
                                             one_alignment_only = True)[0][:2]
                self.left_sequence['aligned'] = [x for x in l]
                self.right_sequence['aligned'] = [x for x in r]
                for x,y in zip(l,r):
                    if x==y and x!="-":
                        self.matches += 1
                if self.debug:
                    print self.debug_head+"Total matches: {}".format(self.matches)                        
        except Exception as e:
            if self.debug:
                print self.debug_head+"Alignment failed for {} and {}".format(self.left_name,self.right_name)
            raise AlignException("pairwise aln failed",
                                 ":".join("{}".format(x) for x in [e, e.args]))
        
        self.x_to_dash()
                                    
    def add(self, sequence, side, name=None):
        '''
        Adds a sequence to a given side of the alignment
        accounting for gaps as "-" in the alignment sequence
        Takes sequence that is either string or list
        and a side ("left" or "right")
        and a name for the column header
        and adds to the dict of that column header
        Returns nothing
        '''
        currentseq = self.left_sequence if side=="left" else self.right_sequence
        if self.debug:
            print self.debug_head+"Adding seq {} to {} accounting for gaps".format(name,side)
        naa = currentseq['nAA']

        if naa!=len(sequence):
            if self.debug:
                print self.debug_head+"Failed to add sequence, naa!=len(sequence)"
            raise AlignException("adding to aln failed",
                  "adding alignment with {} AA's a seq of length {}".format(
                                                                      naa,
                                                                      len(sequence)))
        if name is None:
                name = len(currentseq) + 1
        if name in currentseq:
            print "Warning, attempting to add sequence with name already\
                  nothing added"
            return None
        da = 0
        if self.identical:
            if self.debug:
                print self.debug_head+"sequences are identical"
            newseq = [x for x in sequence]
        else:
            newseq = list()
            aln = currentseq['aligned']
            i = -1
            for x in aln:
                if x!="-":
                    i += 1
                    newseq.append(sequence[i])
                else:
                    da += 1
                    newseq.append(" ")                    
        if self.debug:
            print self.debug_head+"Added {} gap positions".format(da)
        currentseq[name] = newseq
        self.current = False
        return None

    def refresh_df(self):
        '''
        Refreshes the dataframe since I believe it's faster
        to populate dict and remake df than to add columns to premade df
        Before making the df, it looks at keys to see if any are identical
        between the two sides. If so, adds a left/right to them
        '''
        if self.debug:
            print self.debug_head+"Refreshing dataframe"
        newdict = dict()
        lk = set(self.left_sequence.keys())
        rk = set(self.right_sequence.keys())
        sharedkeys = lk.intersection(rk)
        if self.debug:
            print self.debug_head+"refreshing left sequence"
        for k in self.left_sequence.keys():                
            if k == 'aligned':
                newkey = "{}_aa".format(self.left_name)
            elif k in sharedkeys:
                newkey = "{}_left".format(k)
            else:
                newkey = k
            newdict[newkey] = self.left_sequence[k]
        if self.debug:
            print self.debug_head+"refreshing right sequence"
        for k in self.right_sequence.keys():
            if k == 'aligned':
                newkey = "{}_aa".format(self.right_name)              
            elif k in sharedkeys:
                newkey = "{}_right".format(k)
            else:
                newkey = k                            
            newdict[newkey] = self.right_sequence[k]
#        for x in newdict:
#            if isinstance(newdict[x],int): continue
#            print x,newdict[x],len(newdict[x])
        self._dict = newdict            
        self._df = pd.DataFrame.from_dict(newdict)
        if self.debug:
            print self.debug_head+"refreshed dataframe to {} rows".format(len(self._df.index))
        self.current = True

    def remove(name,side="both"):
        '''
        removes a given name from the selected side
        takes name of key and side (left,right,both)
        returns nothing
        '''
        if self.debug:
            print self.debug_head+"removing {} from {}".format(name,side)
        if side=="left" or side=="both":
            if name in self.left_sequence.keys():
                self.left_sequence.pop(name)
                self.current = False
        if side=="right" or side=="both":
            if name in self.right_sequence.keys():
                self.right_sequence.pop(name)                                     
                self.current = False
                                            
    def identity(self, side="left"):
        '''
        returns sequence identity depending on how it's defined
        takes side (left[total is len left],right, outer[total is max len],
        inner [total is min len]
        return float
        '''
        if self.debug:
            print self.debug_head+"getting identity with respect to {}".format(side)
        if side=="left":
            m = self.left_sequence['nAA']
        elif side=="right":
            m = self.right_sequence['nAA']
        elif side=="outer":
            m = max(self.left_sequence['nAA'],self.right_sequence['nAA'])
        else:
            m = min(self.left_sequence['nAA'],self.right_sequence['nAA'])
        if self.debug:
            print self.debug_head+"{} matches over {} residues".format(self.matches,m)
        return float(self.matches)/m
        
    @property
    def df(self):
        if self.current:
            return self._df
        else:
            self.refresh_df()
            return self._df
    
    @property
    def dict(self):
        if not self.current:
            self.refresh_df()
        return copy.deepcopy(self._dict)

    def x_to_dash(self):
        '''
        Replace any X in sequence with -
        '''
        if self.debug:
            print self.debug_head+"replacing X with -"
        for i in range(len(self.left_sequence['aligned'])):
            if self.left_sequence['aligned'][i]=="X":
                self.left_sequence['aligned'][i]="-"
            if self.right_sequence['aligned'][i]=="X":
                self.right_sequence['aligned'][i]="-"
                
                       
def filter_sequences(variants,destination,datasets,debug):
    debug_head = "DEBUG: alignments: filter_sequences: "
    filtered = dict()
    completed = list()
    debug_set = set()
    lost = 0
    for var in variants:
        varcodes = [x[-1] for x in variants[var][-1]]
        if destination=='transcripts' and var not in datasets['transcripts']:
            IO.write_failures(
                [var,variants[var][0][0]],"no transcript seq found for {}".format(var))
            completed += varcodes
            if debug:
                debug_set.add(var)
                lost += 1        
        elif destination=='uniprots' and variants[var][0][0] not in datasets['uniprots']:
            IO.write_failures(
                [var,variants[var][0][0]],"no unp seq found for {}".format(variants[var][0][0]))
            # Only transcript sequence is required for custom models, if unp missing, give it transcript
            print variants[var]
            if variants[var] in datasets['custom_models']:
                datasets['uniprot'][variants[var][0][0]] = datasets['custom_models'][var]
                continue
            completed += varcodes
            if debug:
                debug_set.add(var+"_"+variants[var][0][0])
                lost += 1
        else:
            filtered[var] = variants[var]
    if debug:
        print debug_head+"{} variants from missing {} were lost:".format(lost, destination)
        print debug_head+", ".join(debug_set)
        print debug_head+"{} transcripts remain".format(len(filtered))
    return (filtered,completed)
