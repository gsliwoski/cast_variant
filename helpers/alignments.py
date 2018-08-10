import multiprocessing as mp
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from exceptions import *
import IO
import traceback

# Alignment parameters
GAP_OPEN = -10
GEP_EXTEND = -0.5
MATRIX = matlist.blosum62

def align(seq1, seq2):
    '''
    Aligns two sequences and returns one alignment
    '''
    try:
        if seq1==seq2:
            aligns = [seq1,seq2]
        else:
            aligns = pairwise2.align.globalds(
                                             seq1,
                                             seq2,
                                             MATRIX,
                                             GAP_OPEN,
                                             GAP_EXTEND,
                                             one_alignment_only = True)[0][:2]
    except Exception as e:
        raise AlignException("pairwise aln failed",
                             ":".join("{}".format(x) for x in [e,e.args]))
    return aligns

def number_alignment(seq1, seq2, num1=None, num2=None):
    '''
    Aligns sequence numbering with alignment to account for gaps
    Takes two alignments and two residue numbering lists and returns
    two new lists with " " for each gap
    '''
    # Sequences are alignments so must be same length
    if len(seq1)!=len(seq2):
        raise AlignException("number alignment","seq's not same length")
    # If no resnums provided, assume they are 1 - length
    if num1 is None:
        num1 = range(1,len(seq1)+1)
    if num2 is None:
        num2 = range(1,len(seq2)+1)
    if seq1==seq2 and len(seq1)==len(num1) and len(num1)==len(num2):
        return (num1,num2,1.0)
    aindex,bindex = -1,-1
    anumbers = list()
    bnumbers = list()
    matches = 0
    total = 0
    try:
        for x,y in zip(aseq,bseq):
            if x!="-":
                if x==y: matches += 1
                if y!="-": total += 1
                aindex += 1
                current_a = anum[aindex]
            else:
                current_a = " "
            if y!="-":
                bindex += 1
                current_b = bnum[bindex]
            else:
                current_b = " "
            anumbers.append(current_a)
            bnumbers.append(current_b)
    except IndexError as e:
        raise AlignError("number_alignment index error",
                         ":".join("{}".format(x) for x in [e,e.args]))
    return (anumbers,bnumbers,float(matches)/total*100)
    
def aligner(seq1,seq2, num1=None, num2=None):
    '''
    Takes two sequences and residue number lists
    and calls align and number align and returns
    aligned sequences and residue numbers for each and identity
    '''
    aln1, aln2 = align(seq1, seq2)                                                              
    renum1,renum2,identity = number_alignment(aln1,aln2,num1,num2)
    return (aln1,aln2,renum1,renum2,identity)
                                               
def pdb_worker(variant,transcripts,uniprots,sifts,multiproc):
    '''
    Worker process, handles calls to alignments and calls to final output
    Takes a current variant and True/False multiprocessing and returns nothing
    '''
    if multiproc:
        lock = mp.Lock()
    else:
        lock = None
    current_transcript, current_protein, variant_list = variant
    failmsg = list()
    try:
        uniprot,isoform = current_protein
        astring = "{}_{}".format(current_transcript,uniprot)
        trans_seq = transcripts.get(current_transcript,"")
        if trans_seq == "":
            raise AlignException("{}".format(current_transcript),"no seq found")
        unp_seq = uniprots.get(uniprot,"")
        if unp_seq == "":
            raise AlignException("{}".format(uniprot),"no seq found")
        current_unp = {'fasta': unp_seq,
                       'resnums': range(1,len(unp_seq)+1),
                       'uniprot': uniprot,
                       'isoform': isoform,
                       'transcript': current_transcript}
        trans_aln, unp_aln, trans_num, unp_num, identity = aligner(
                                                           trans_seq,
                                                           unp_seq)  
        current_unp['transcript_identity'] = identity
        current_unp['transcript_sequence'] = trans_aln
        current_unp['transcript_position'] = trans_num
        current_unp['uniprot_sequence'] = unp_aln
        current_unp['uniprot_position'] = unp_num
    except AlignException as e:
        IO.write_failures([current_transcript,uniprot],e.fullmsg,lock)
        return False
    except Exception:
        IO.write_failures([current_transcript,uniprot],traceback.format_exc(),lock)
        return False
    return True
def init(l):
    '''Shared lock initializer, used during pool init'''
    global lock
    lock = l

def cast_to_pdbs(variants, transcripts, uniprot, sifts, nproc):
    if nproc>1:
        lock = mp.Lock()
        pool = mp.Pool(processes=nproc,initializer=init,initargs=(l,))
    results = list()
    for x in variants.keys():
        current_var = [x]+variants[x]
        if nproc>1:
            results = pool.apply_async(pdb_worker,
                                       args = (current_var,
                                               transcripts,
                                               uniprot,
                                               sifts,
                                               True))
        else:
            results.append(pdb_worker(current_var,
                                      transcripts,
                                      uniprot,
                                      sifts,
                                      False))
    return (1,1,sum(results))
def cast_to_models(variants, transcripts, models):
    pass
