import multiprocessing as mp
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from exceptions import *
import IO
import traceback
import pandas as pd

##### MISC #####

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

# TODO: Make this a df intersection instead
# TODO: allowing you to get rid of t quickfix
def number_alignment(seq1, seq2, num1=None, num2=None, trailer=None):
    '''
    Aligns sequence numbering with alignment to account for gaps
    Takes two alignments and two residue numbering lists and returns
    two new lists with " " for each gap
    trailer is optional sequence to be adjusted with seq2 (icodes)
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
        return (num1,num2,1.0,trailer)
    aindex,bindex = -1,-1
    anumbers = list()
    bnumbers = list()
    matches = 0
    total = 0
    t = list
    trailer = [x for x in trailer] if trailer else None
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
                t = trailer.pop(0) if trailer else '-'
            else:
                current_b = " "
                current_t = "-"
            anumbers.append(current_a)
            bnumbers.append(current_b)
            t.append(current_t)
    except IndexError as e:
        raise AlignError("number_alignment index error",
                         ":".join("{}".format(x) for x in [e,e.args]))
    t = t if trailer else None
    return (anumbers,bnumbers,float(matches)/total*100,t)
   
def aligner(seq1,seq2, num1=None, num2=None):
    '''
    Takes two sequences and residue number lists
    and calls align and number align and returns
    aligned sequences and residue numbers for each and identity
    '''
    aln1, aln2 = align(seq1, seq2)                                                              
    renum1,renum2,identity,trailer = number_alignment(aln1,aln2,num1,num2)
    return (aln1,aln2,renum1,renum2,identity,trailer)

def expand_variants(variants):
    '''
    Expands non-missense variants into 
    all affected residues
    '''
    pass

def filter_sequences(variants,transcripts,uniprot):
    filtered = dict()
    for var in variants:
        if var not in transcripts:
            IO.write_failures(
                [var,variants[var][0][0]],"no seq found for {}".format(var))
        elif variants[var][0][0] not in uniprot:
            IO.write_failures(
                [var,variants[var][0][0]],"no seq found for {}".format(var[0][0]))
        else:
            filtered[var] = variants[var]
    return filtered            

                                                                        
##### PDB #####
   
def pdb_worker(variant,transcripts,uniprots,sifts,multiproc,expand):
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
        varcodes = set([x[3] for x in variant_list])
        #TODO: these checks can probably be removed
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
        trans_aln, unp_aln, trans_num, unp_num, identity,_ = aligner(
                                                           trans_seq,
                                                           unp_seq)  
        current_unp['transcript_identity'] = identity
        current_unp['transcript_aa'] = trans_aln
        current_unp['transcript_position'] = trans_num
        current_unp['uniprot_aa'] = unp_aln
        current_unp['uniprot_position'] = unp_num
        
        # Get the sifts entries for this uniprot
        pdbs = gather_sifts(uniprot, sifts)
        if pdbs.shape[0]==0:    
            raise AlignException("sifts parsing", "no residues returned")

        #Expand variants if set
        if expand:
            expand_variants(variant)        
        # Compile the output dataframes
        aln_table = IO.generate_alignment_table(current_unp, pdbs)
        
        variant_table = IO.generate_variant_table(variant_list)
        conv_table = IO.generate_conversion_table(aln_table, variant_table)
        
    except AlignException as e:
        IO.write_failures([current_transcript,uniprot],e.fullmsg,lock)
        if multiproc:
            lock.acquire()
        IO.write_complete(varcodes)
        if multiproc:
            lock.release()       
        return False
    except Exception:
        IO.write_failures([current_transcript,uniprot],traceback.format_exc(),lock)
        return False
    # Safe-write the output
    if multiproc:
        lock.acquire()
    IO.write_output(aln_table,conv_table)                                                
    IO.write_complete(varcodes)
    if multiproc:
        lock.release()
    print "{} complete.".format(astring)
    return True
    
def init(l):
    '''Shared lock initializer, used during pool init'''
    global lock
    lock = l

def cast_to_pdbs(variants, transcripts, uniprot, sifts, nproc,expand):
    if nproc>1:
        l = mp.Lock()
        pool = mp.Pool(processes=nproc,initializer=init,initargs=(l,))
    results = list()
    for x in variants.keys():
        current_var = [x]+variants[x]
        if nproc>1:
            results.append(pool.apply_async(pdb_worker,
                                       args = (current_var,
                                               transcripts,
                                               uniprot,
                                               sifts,
                                               True,
                                               expand)))
        else:
            results.append(pdb_worker(current_var,
                                      transcripts,
                                      uniprot,
                                      sifts,
                                      False,
                                      expand))
    if nproc>1:
            pool.close()
            pool.join()
            results = [x.get() for x in results]
        
    return sum(results)

##### SIFTS #####

# holds PDB entries from sifts alignments in case they are 
# needed again to save time
sifts_holder = dict()

def gather_sifts(uniprot,sifts):
    '''
    Get all the sifts alignments for the given uniprot
    Takes a string uniprot from the sifts dict
    Returns pandas df of all the PDB residues for that uniprot
    '''
    pdbs = sifts.get(uniprot,None)

    all_residues = None
    if pdbs is None:
        raise AlignException(
              "sifts_lookup","no sifts entry for {}".format(uniprot))
    #print "{} pdbs found for {}".format(len(pdbs),uniprot)
    for pdb in pdbs:
        residues = sifts_holder.get(pdb[0],None)
        if residues is None:
            residues = IO.parse_sifts(pdb)
            if residues is None:
                print "Warning, None returned when parsing sifts for {}".format(pdb)
                continue
            sifts_holder[pdb[0]] = residues
        # TODO: improve?
        # Store residues as dataframe and filter for relevant chain
        resdf = pd.DataFrame.from_dict(residues)
        # Filter for only the relevant chain
        resdf = resdf[(resdf['chain'] == pdb[1]) &
                      ((resdf["uniprot"] == uniprot) |
                       (resdf["uniprot"] == ""))]
        if len(resdf.index)==0:
            continue
        resdf.drop("uniprot",inplace=True,axis=1)
        if all_residues is None:
           all_residues = resdf
        else:
            all_residues = pd.concat([all_residues,resdf],copy=False)            
    if all_residues is None:
        all_residues = pd.DataFrame.from_dict(None)
    all_residues.to_csv(open("tempstr","w"))
    return all_residues
##### MODELS #####
    
def cast_to_models(variants, transcripts, models):
    pass
