import multiprocessing as mp
import pandas as pd
import alignments
import IO
import traceback
from exceptions import *

VAR_STRUCT_HEADER = ["transcript","uniprot","isoform",
               "transcript_identity","transcript_position",
               "uniprot_position", "structure_position", "icode",
               "transcript_aa", "uniprot_aa", "structure_aa",
               "ref_aa", "alt_aa", "structure_identity","template_identity",
               "structure", "chain","structure_isoform",
               "complex_state","varcode"]

ALN_STRUCT_HEADER = ["transcript","uniprot","isoform",
              "transcript_identity","transcript_position",
              "uniprot_position","structure_position",
              "transcript_aa","uniprot_aa","structure_aa",
              "structure_identity", "template_identity",
              "structure","chain","structure_isoform","complex_state"]

VAR_UNP_HEADER = ["transcript","uniprot","isoform",
                  "transcript_identity","transcript_position",
                  "uniprot_position","transcript_aa","uniprot_aa",
                  "ref_aa","alt_aa","varcode"]

ALN_UNP_HEADER = ["transcript","uniprot","isoform",
                  "transcript_identity","transcript_position",
                  "uniprot_position","transcript_aa","uniprot_aa"]

def init(l):
    '''Shared lock initializer, used during pool init'''
    global lock
    lock = l

def cast_variants(variants, datasets,arguments):
    if arguments.num_procs>1:
        l = mp.Lock()
        pool = mp.Pool(processes=arguments.num_procs,
                       initializer=init,
                       initargs=(l,))
    results = list()                       
    for x in variants.keys():
        current_var = [x]+variants[x]
        if arguments.num_procs>1:
            results.append(pool.apply_async(worker,
                                       args = (current_var,
                                               datasets,
                                               arguments)))
        else:
            results.append(worker(current_var,
                                  datasets,
                                  arguments))
    if arguments.num_procs>1:
            pool.close()
            pool.join()
            results = [x.get() for x in results]
        
    return sum(results)
                                   
def worker(variant, datasets, arguments):
    ''''
    Worker process that does all alignment work
    Takes the variant, a dictionary of datasets, and args object
    Returns nothing
    '''
    alntables = list()
    vartables = list()
    completed = list()
    multiproc = arguments.num_procs > 1
    if multiproc:
        lock = mp.Lock()
    else:
        lock = None
    current_transcript, current_protein, variant_list = variant
    try:
        uniprot,isoform = current_protein
        astring = "{}_{}".format(current_transcript,uniprot)
        transcripts = datasets['transcripts']
        trans_seq = datasets['transcripts'][current_transcript]
        unp_seq = datasets['uniprots'][uniprot]
        #varcodes set is used to write completion table
        varcodes = set([x[4] for x in variant_list])
        
        # Expand the nonmissense variants if set
        if arguments.expand:
            variant_list = expand_variants(variant_list,trans_seq)
        
        variant_table = generate_variant_table(variant_list)
        
        if not arguments.nopdb or not arguments.nouniprot:
            # Generate uniprot tables
            current_aln = alignments.Alignment(trans_seq, 
                                               unp_seq, 
                                               'transcript', 
                                               'uniprot')
            current_aln.add(range(1,1+len(trans_seq)),
                            'left',
                            'transcript_position')
            current_aln.add(range(1,1+len(unp_seq)),
                            'right',
                            'uniprot_position')
            current_unp = current_aln.dict
            current_unp['transcript_identity'] = current_aln.identity("left")
            current_unp['uniprot'] = uniprot
            current_unp['isoform'] = isoform
            current_unp['transcript'] = current_transcript
            unp_df = pd.DataFrame.from_dict(current_unp).round({'transcript_identity':1})
#            var_unp_df = unp_df.merge(variant_table,how='inner',on=['transcript_position'])
#            var_unp_df.sort_values(by=["transcript_position","varcode"],inplace=True)                
            # Unless supressing the just uniprot, also include the any vars that hit a unp res
            if not arguments.nouniprot:
                current_unp['structure'] = "Uniprot"
                current_unp['chain'] = "-"
                current_unp['structure_aa'] = "-"
                current_unp['structure_position'] = 0
                current_unp['icode'] = " "
                current_unp['template_identity'] = 1.0
                current_unp['structure_identity'] = 1.0
                current_unp['structure_isoform'] = isoform
                current_unp['complex_state'] = "Uniprot"
                unp_df_out = pd.DataFrame.from_dict(current_unp).round({'transcript_identity':1})
                unp_df_out.sort_values(by=["structure","chain","transcript_position",
                                           "structure_position","icode"],inplace=True)
                var_unp_df = unp_df_out.merge(variant_table,how='inner',on=['transcript_position'])
                var_unp_df.sort_values(by=['structure','chain','transcript_position',
                                           'structure_position','icode'],inplace=True)                                                          
                alntables.append([unp_df_out,ALN_STRUCT_HEADER])
                vartables.append([var_unp_df,VAR_STRUCT_HEADER])           

        if not arguments.nopdb:
            # Get the PDB table from SIFTS
            pdbs = gather_sifts(uniprot, isoform, datasets['sifts'])
            if pdbs.shape[0] == 0:
                raise WorkerException("sifts parsing", "no residues returned")
            # Generate the PDB tables
            pdb_df = unp_df.merge(pdbs,how="inner",on=['uniprot_position']).round({'structure_identity':1})
            pdb_df.sort_values(by = ["structure","chain","transcript_position",
                                     "structure_position","icode"],inplace=True)
            var_pdb_df = pdb_df.merge(variant_table,how='inner',on=['transcript_position'])
            var_pdb_df = var_pdb_df[var_pdb_df['structure_position']!=' ']
            var_pdb_df.sort_values(by = ["structure","chain","transcript_position",
                                     "structure_position","icode","varcode"],inplace=True)
            alntables.append([pdb_df,ALN_STRUCT_HEADER])
            vartables.append([var_pdb_df,VAR_STRUCT_HEADER])
        for x in ALN_STRUCT_HEADER:
            if x not in list(pdb_df):
                print "{} missing from pdbdf".format(x)
        for x in VAR_STRUCT_HEADER:
           if x not in list(var_pdb_df):
              print "{} missing from varpdbdf".format(x)               
        for x in ALN_STRUCT_HEADER:
            if x not in list(unp_df_out):
                print "{} missing from unpdf".format(x)
        for x in VAR_STRUCT_HEADER:
           if x not in list(var_unp_df):
              print "{} missing from varunpdf".format(x)               

    except WorkerException as e:
        if multiproc:
            lock.acquire()
        IO.write_failures([current_transcript, uniprot], e.fullmsg)          
        if multiproc:
            lock.release()        
    except AlignException as e:
        if multiproc:
            lock.acquire()
        IO.write_failures([current_transcript, uniprot], e.fullmsg)
        if multiproc:
            lock.release()
        return False
    except Exception:
        if multiproc:
            lock.acquire()
        IO.write_failures([current_transcript, uniprot], traceback.format_exc())
        if multiproc:
            lock.release()
        return False                                            

    # Move on to swissmodels
    try:
        if not arguments.nomodel:
            # Generate the model table
            modelids = gather_models(uniprot,isoform,datasets['models'])
            if len(modelids)==0:
                raise WorkerException("model lookup",
                                      "no models for {}".format(uniprot))
            all_residues = None
            for model in modelids:
                current_model = IO.load_model(model[2])
                if len(current_model['fasta']) == 0:
                    print "Warning, no residues in model {}".format(
                                            current_model['filename'])
                    continue
                model_seq = current_model['fasta']
                
                current_aln = alignments.Alignment(trans_seq, 
                                                   model_seq, 
                                                   'transcript', 
                                                   'structure')
                current_aln.add(range(1,1+len(trans_seq)),
                                'left',
                                'transcript_position')
                current_aln.add(current_model['resnums'],
                                'right',
                                'structure_position')
                current_aln.add(current_model['icodes'],
                                'right',
                                'icode')                                               
                current_mod = current_aln.dict
                current_mod.pop('transcript_aa')
                current_mod['structure_identity'] = current_aln.identity("left")
                current_mod['structure'] = current_model['filename']
                current_mod['chain'] = current_model['chain']
                current_mod['complex_state'] = current_model['complex_state']
                current_mod['template_identity'] = model[1]
                current_mod['structure_isoform'] = model[0]
                current_df = pd.DataFrame.from_dict(current_mod)
                if all_residues is None:
                   all_residues = current_df
                else:
                    all_residues = pd.concat([all_residues,current_df],copy=False)            
                
            if all_residues is None:
                raise WorkerException("model parsing","no model residues")                           
#            print "modelkeys: {}".format(list(all_residues))
#            print "unpkeys: {}".format(list(unp_df))
#            print "alnstructhead: {}".format(ALN_STRUCT_HEADER)
#            print "varstructhead: {}".format(VAR_STRUCT_HEADER)
            model_df = all_residues.merge(unp_df,how="left",on=['transcript_position']).round({'structure_identity':1})
            model_df.sort_values(by = ["structure","chain","transcript_position",
                                     "structure_position","icode"],inplace=True)
           
            var_model_df = model_df.merge(variant_table,how='inner',on=['transcript_position'])
            var_model_df = var_model_df[var_model_df['structure_position']!=' ']
            var_model_df.sort_values(by = ["structure","chain","transcript_position",
                                         "structure_position","icode","varcode"],inplace=True)     
            alntables.append([model_df,ALN_STRUCT_HEADER])
            vartables.append([var_model_df,VAR_STRUCT_HEADER])
    except WorkerException as e:
        if multiproc:
            lock.acquire()
        IO.write_failures([current_transcript, uniprot], e.fullmsg)          
        if multiproc:
            lock.release()        
    except AlignException as e:
        if multiproc:
            lock.acquire()
        IO.write_failures([current_transcript, uniprot], e.fullmsg)
        if multiproc:
            lock.release()
        return False
    except Exception:
        if multiproc:
            lock.acquire()
        IO.write_failures([current_transcript, uniprot], traceback.format_exc())
        if multiproc:
            lock.release()
        return False   
                
    # Finished alignments, now write tables
    if multiproc:
        lock.acquire()
    if not arguments.noalign:
        for table in alntables:
            IO.write_table(table,"alignments")
    for table in vartables:
        IO.write_table(table,"variants")
    if arguments.completed:
        IO.write_complete(varcodes)
    if multiproc:
        lock.release()
    return True                                                                                       

def generate_variant_table(variant_list):
    '''
    Generates a variant dataframe
    from a variant list
    '''
    variant_dict = {'transcript_position': [x[0] for x in variant_list],
                    'ref_aa': [x[1] for x in variant_list],
                    'alt_aa': [x[2] for x in variant_list],
                    'varcode': [x[4] for x in variant_list]}
    variant_df = pd.DataFrame.from_dict(variant_dict)
    return variant_df

def gather_models(unp,isoform,models):
    '''Gather the appropriate model files
       If there is a isoform that matches, use those models
       If not, use all models
       Returns a list of models in format:
       [isoform,identity,filename]'''
    current_models = list()
    if unp in models:
        for model in models[unp]:
            if model[0]==isoform:
                current_models.append(model)
            #In case the -1 has been inferred as no isoform designation
            #There shouldn't be a -1
            elif (isoform.endswith("-1") and 
                 len(
                 [x[0] for x in models[unp] if x[0].endswith("-1")]
                 )==0 and
                 "-" not in model[0]):
                current_models.append(model)                 
        if len(current_models) == 0: # If no isoform was found, just use them all
            current_models = models[unp]
    if len(current_models)==0:
        raise WorkerException(
              "swissmodel_lookup","no models entry for {}".format(uniprot))
    else:
        return current_models
        
# holds PDB entries from sifts alignments in case they are 
# needed again to save time
sifts_holder = dict()

def gather_sifts(uniprot,isoform,sifts):
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
        residues['structure_isoform'] = 'PDB'
        residues['template_identity'] = 1.0
        residues['complex_state'] = 'PDB'
        resdf = pd.DataFrame.from_dict(residues)
        # Filter for only the relevant chain
        resdf = resdf[(resdf['chain'] == pdb[1]) &
                      ((resdf["uniprot"] == uniprot) |
                       (resdf["uniprot"] == ""))]
        if len(resdf.index)==0:
            continue
        resdf = resdf.drop("uniprot",inplace=False,axis=1)
        if all_residues is None:
           all_residues = resdf
        else:
            all_residues = pd.concat([all_residues,resdf],copy=False)            
    if all_residues is None:
        all_residues = pd.DataFrame.from_dict(None)
    all_residues.to_csv(open("tempstr","w"))
    return all_residues

def expand_variants(variants,trans_seq):
    '''
    Expands non-missense variants into 
    all affected residues
    Takes a list of variants as
    [position, ref, alt,annotation,varcode]
    and a transcript sequence
    Expands inframe deletions with an - for alt at every pos lost
    Expands early stop and frameshifts with X for alt at every pos
    after start of mutation
    Drops any variants that can not be expanded, writing them to failures
    Returns new list of variants including expanded vars when possible
    '''
    failures = list()
    newvariants = list()
    reason = "Failed to expand"
    for variant in variants:
        position,refaa,altaa,annotation,varcode = variant
        # Don't expand missense        
        if "MISSENSE" in annotation:
            newvariants.append(variant)
            continue
        # Expand in-frame deletions
        elif "-" in position and "INFRAME_DELETION" in annotation:
            newref = refaa
            if len(altaa)>1:
                write.failures(variant,reason)
                continue
            try:                
                pos1,pos2 = [int(x) for x in variant[0].split("-")]
            except IndexError,ValueError:
                write.failures(variant,reason+": parse position")
                continue
            # Inframe deletions are sometimes defined with one AA remaining
            # instead of just a - in altaa
            if altaa==refaa[0]:
                pos1 += 1
                newref = refaa[1:]
            elif altaa==refaa[-1]:
                pos2-=1
                newref = refaa[:-1]
            for i,x in enumerate(range(pos1,pos2+1)):   
                newvariants.append([i,newref,"-",annotation,varcode])
        # Single residue deletions
        elif "INFRAME_DELETION" in annotation and "STOP" not in annotation:
            newvariants.append([int(position),refaa,"-",annotation,varcode])
        # Everything else (frameshift, early stop) are from pos to end of transcript
        else:
            startpos = int(position.split("-")[0])
            for pos in range(startpos-1,len(trans_seq)):
                newvariants.append([pos+1,refaa,"X",annotation,varcode])
    return newvariants                            
                                                                                         
