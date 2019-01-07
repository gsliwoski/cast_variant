import multiprocessing as mp
import pandas as pd
import alignments
import IO
import traceback
from exceptions import *
from descriptors import Structure, add_structure_descriptors, add_uniprot_descriptors, HEADERS
import time
import sys

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
    debug = arguments.debug
    debug_head = "DEBUG: workers: cast_variants: "
    if arguments.num_procs>1:
        if debug:
            print debug_head+"generating pool for {} procs".format(arguments.num_procs)
        l = mp.Lock()
        pool = mp.Pool(processes=arguments.num_procs,
                       initializer=init,
                       initargs=(l,))
    results = list()                       
    print "Casting variants across {} transcripts".format(len(variants.keys()))
    for x in variants.keys():
        current_var = [x]+variants[x]
        if arguments.num_procs>1:
            if debug:
                print debug_head+"running {} in pool".format(current_var[0])
            results.append(pool.apply_async(worker,
                                       args = (current_var,
                                               datasets,
                                               arguments)))
        else:
            if debug:
                print debug_head+"running {} in sequence".format(current_var[0])
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
    worker_start_time = time.time()
    debug = arguments.debug
    debug_head = "DEBUG: workers: worker({}): ".format(variant[0])
    alntables = list()
    vartables = list()
    completed = list()
    multiproc = arguments.num_procs > 1
    if debug:
        ct = time.localtime()
        print debug_head+"Starting worker at {}:{}:{}".format(ct.tm_hour,ct.tm_min,ct.tm_sec)
    if multiproc:
        lock = mp.Lock()
    else:
        lock = None
    current_transcript, current_protein, variant_list = variant

    print "Casting {} variants".format(current_transcript)
    var_df_header = VAR_STRUCT_HEADER[:]
    if debug:
        print debug_head+"Updating header with descriptors start: {}".format(var_df_header)
    for d in arguments.descriptors:
        var_df_header += HEADERS.get(d,[])
    if debug:
        print debug_head+"Final header: {}".format(var_df_header)
    try:
        uniprot,isoform = current_protein
        astring = "{}_{}".format(current_transcript,uniprot)
        if current_transcript == "ENST00000359218": #TODO: Deal with extremely long transcripts like Titin?
            print "Skipping {}, Titin is currently not allowed".format(current_transcript)
            raise WorkerException("loading stage","Titin ({}) is prevented".format(current_transcript))
        if debug:
            ce = round(time.time()-worker_start_time,1)
            print debug_head+"Getting sequences for {} [{} secs elapsed]".format(astring,ce)
        transcripts = datasets['transcripts']
        trans_seq = datasets['transcripts'][current_transcript]
        unp_seq = datasets['uniprots'][uniprot]
        #varcodes set is used to write completion table
        varcodes = set([x[4] for x in variant_list])
        if debug:
            ce = round(time.time()-worker_start_time,1)
            print debug_head+"processing {} variants [{} secs elapsed]".format(len(varcodes),ce)
        # Expand the nonmissense variants if set
        if arguments.expand:
            variant_list = expand_variants(variant_list,trans_seq)
        if debug:
            ce = round(time.time()-worker_start_time,1)        
            print debug_head+"generating variant table from variant_list [{} secs elapsed]".format(ce)
        variant_table = generate_variant_table(variant_list)
        
        if not arguments.nopdb or not arguments.nouniprot:
            # Generate uniprot tables
            if debug:
                ce = round(time.time()-worker_start_time,1)                
                print "aligning sequences for {} [{} secs elapsed]".format(astring,ce)
            current_aln = alignments.Alignment(trans_seq, 
                                               unp_seq, 
                                               debug,
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
            if debug:
                ce = round(time.time()-worker_start_time,1)            
                print debug_head+"Generating dictionary from {} alignment [{} secs elapsed]".format(astring,ce)
            unp_df = pd.DataFrame.from_dict(current_unp).round({'transcript_identity':1})
#            var_unp_df = unp_df.merge(variant_table,how='inner',on=['transcript_position'])
#            var_unp_df.sort_values(by=["transcript_position","varcode"],inplace=True)                
            # Unless supressing the just uniprot, also include the any vars that hit a unp res
            if not arguments.nouniprot:
                if debug:
                    ce = round(time.time()-worker_start_time,1)                
                    print debug_head+"Generating uniprot output table [{} secs elapsed]".format(ce)
                current_unp['structure'] = "Uniprot"
                current_unp['chain'] = "-"
                current_unp['structure_aa'] = "-"
                current_unp['structure_position'] = current_unp['uniprot_position']
                current_unp['icode'] = " "
                current_unp['template_identity'] = 100.0
                current_unp['structure_identity'] = 100.0
                current_unp['structure_isoform'] = isoform
                current_unp['complex_state'] = "Uniprot"
                unp_df_out = pd.DataFrame.from_dict(current_unp).round({'transcript_identity':1})
                unp_df_out.sort_values(by=["structure","chain","transcript_position",
                                           "structure_position","icode"],inplace=True)
                if debug:
                    ce = round(time.time()-worker_start_time,1)                
                    print debug_head+"Alignment df {} rows [{} secs elapsed]".format(len(unp_df_out.index),ce)
                var_unp_df = unp_df_out.merge(variant_table,how='inner',on=['transcript_position'])
                if debug:
                    ce = round(time.time()-worker_start_time,1)
                    print debug_head+"Variant df {} rows [{} secs elapsed]".format(len(var_unp_df.index),ce)
                #Add the descriptors (holders for unp since not structure)
                var_unp_df.sort_values(by=['structure','chain','transcript_position',
                                           'structure_position','icode'],inplace=True)                                                          
       
                alntables.append([unp_df_out,ALN_STRUCT_HEADER])
                vartables.append([var_unp_df,var_df_header])           
                
        if not arguments.nopdb:
            # Get the PDB table from SIFTS
            if debug:
                ce = round(time.time()-worker_start_time,1)            
                print debug_head+"Getting PDBs from sifts [{} secs elapsed]".format(ce)
            pdbs = gather_sifts(uniprot, isoform, datasets['sifts'])
            if len(pdbs.index) == 0:
                raise WorkerException("sifts parsing", "no residues returned")
            # Generate the PDB tables
            if debug:
                ce = round(time.time()-worker_start_time,1)
                print debug_head+"Retrieved {} rows from sifts [{} secs elapsed]".format(len(pdbs.index),ce)
            pdb_df = unp_df.merge(pdbs,how="inner",on=['uniprot_position']).round({'structure_identity':1})
            pdb_df.loc[pdb_df.icode == "", 'icode'] = " "
            pdb_df.sort_values(by = ["structure","chain","transcript_position",
                                     "structure_position","icode"],inplace=True)
            if debug:
                ce = round(time.time()-worker_start_time,1)            
                print debug_head+"Merged with uniprot for {} rows [{} secs elapsed]".format(len(pdb_df.index),ce)
            var_pdb_df = pdb_df.merge(variant_table,how='inner',on=['transcript_position'])
            var_pdb_df = var_pdb_df[var_pdb_df['structure_position']!=' ']
            var_pdb_df.structure_position = pd.to_numeric(var_pdb_df.structure_position,errors='coerce')
            var_pdb_df.sort_values(by = ["structure","chain","transcript_position",
                                     "structure_position","icode","varcode"],inplace=True)
            if debug:
                ce = round(time.time()-worker_start_time,1)            
                print debug_head+"Merged with variants for {} rows [{} secs elapsed]".format(len(var_pdb_df.index),ce)
            alntables.append([pdb_df,ALN_STRUCT_HEADER])
            vartables.append([var_pdb_df,var_df_header])

    except WorkerException as e:
        if debug:
            print debug_head+"unp/pdb worker exception {} {}: {}".format(current_transcript,uniprot,e.fullmsg)
        if multiproc:
            lock.acquire()
        IO.write_failures([current_transcript, uniprot], e.fullmsg)          
        if multiproc:
            lock.release()        
    except AlignException as e:
        if debug:
            print debug_head+"unp/pdb align exception {} {}: {}".format(current_transcript,uniprot,e.fullmsg)
        if multiproc:
            lock.acquire()
        IO.write_failures([current_transcript, uniprot], e.fullmsg)
        if multiproc:
            lock.release()
        return False
    except Exception as e:
        if debug:
            print debug_head+"unp/pdb uncaught worker exception {} {}: {}".format(current_transcript,uniprot,e)    
        if multiproc:
            lock.acquire()
        IO.write_failures([current_transcript, uniprot], traceback.format_exc())
        if multiproc:
            lock.release()
        return False                                            

    # Move on to swissmodels
    try:
        if not arguments.nomodel: #TODO: Fix the float structure_residue numbering
            # Generate the model table
            modelids = gather_models(uniprot,isoform,datasets['models'],debug)
            if debug:
                ce = round(time.time()-worker_start_time,1)            
                print debug_head+"{} total models for {} [{} secs elapsed]".format(len(modelids),uniprot,ce)
            if len(modelids)==0:
                raise WorkerException("model lookup",
                                      "no models for {}".format(uniprot))
            all_residues = None
            for model in modelids:
                current_model = IO.load_model(model[2],debug)
                if len(current_model['fasta']) == 0:
                    print "Warning, no residues in model {}".format(
                                            current_model['filename'])
                    continue
                model_seq = current_model['fasta']
                if debug:
                    ce = round(time.time()-worker_start_time,1)    
                    print debug_head+"aligning models to {} [{} secs elapsed]".format(uniprot,ce)
                current_aln = alignments.Alignment(trans_seq, 
                                                   model_seq,
                                                   debug, 
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
                if debug:
                    ce = round(time.time()-worker_start_time,1)
                    print debug_head+"generating model df for {} [{} secs elapsed]".format(uniprot,ce)
                current_mod = current_aln.dict
                current_mod.pop('transcript_aa')
                current_mod['structure_identity'] = current_aln.identity("left")*100
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
                if debug:
                    ce = round(time.time()-worker_start_time,1)                
                    print debug_head+"overall model df for {} has {} rows [{} secs elapsed]".format(uniprot,len(all_residues.index),ce)
            if all_residues is None:
                raise WorkerException("model parsing","no model residues")                           
#            print "modelkeys: {}".format(list(all_residues))
#            print "unpkeys: {}".format(list(unp_df))
#            print "alnstructhead: {}".format(ALN_STRUCT_HEADER)
#            print "varstructhead: {}".format(VAR_STRUCT_HEADER)
            if debug:
                ce = round(time.time()-worker_start_time,1)            
                print debug_head+"merging models with tables for unp {} [{} secs elapsed]".format(uniprot,ce)
            model_df = all_residues.merge(unp_df,how="left",on=['transcript_position']).round({'structure_identity':1})
            model_df.sort_values(by = ["structure","chain","transcript_position",
                                     "structure_position","icode"],inplace=True)
            if debug:
                ce = round(time.time()-worker_start_time,1)            
                print debug_head+"final model with {} has {} rows [{} secs elapsed]".format(uniprot,len(model_df.index),ce)
            var_model_df = model_df.merge(variant_table,how='inner',on=['transcript_position'])
            var_model_df = var_model_df[var_model_df['structure_position']!=' ']
            var_model_df.sort_values(by = ["structure","chain","transcript_position",
                                         "structure_position","icode","varcode"],inplace=True)     
            if debug:
                ce = round(time.time()-worker_start_time,1)            
                print debug_head+"final model var with {} has {} rows [{} secs elapsed]".format(uniprot,len(var_model_df.index),ce)
            alntables.append([model_df,ALN_STRUCT_HEADER])
            vartables.append([var_model_df,var_df_header])
                           
    except WorkerException as e:
        if debug:
            print debug_head+"model worker exception {} {}: {}".format(current_transcript,uniprot,e.fullmsg)
        if multiproc:
            lock.acquire()
        IO.write_failures([current_transcript, uniprot], e.fullmsg)          
        if multiproc:
            lock.release()        
    except AlignException as e:
        if debug:
            print debug_head+"model align exception {} {}: {}".format(current_transcript,uniprot,e.fullmsg)    
        if multiproc:
            lock.acquire()
        IO.write_failures([current_transcript, uniprot], e.fullmsg)
        if multiproc:
            lock.release()
        return False
    except Exception as e:
        if debug:
            print debug_head+"model worker uncaught exception {} {}: {}".format(current_transcript,uniprot,e)
        if multiproc:
            lock.acquire()
        IO.write_failures([current_transcript, uniprot], traceback.format_exc())
        if multiproc:
            lock.release()
        return False   
                
    # Finished alignments, now add descriptors and write tables
    if len(arguments.descriptors)>0:
        if debug:
            ce = round(time.time()-worker_start_time,1)        
            print "Adding descriptors to {} {} [{} secs elapsed]".format(current_transcript,uniprot,ce)
        old_vartables = vartables
        vartables = list()
        dtimes = [time.time()]
        for table in old_vartables:
            if table[0].empty:
                if debug:
                    print debug_head+"Empty table encountered while adding descriptors to {} {}".format(current_transcript,uniprot)
                dtimes.append(time.time()-dtimes[0])                    
                continue
            # Add the structure-based descriptors
            fulltable = add_structure_descriptors(table[0],arguments.descriptors,debug)
            # Add the uniprot feature group descriptors which only rely on uniprot position
            if 'unp' in arguments.descriptors:
                fulltable = add_uniprot_descriptors(fulltable,debug)
            fulltable.sort_values(by = ["structure","chain","transcript_position",
                                        "structure_position","icode","varcode"],inplace=True)
            vartables.append([fulltable,var_df_header])
            dtimes.append(time.time()-dtimes[0])                                                                
        if debug:
            print "{} tables got descriptors with times taking: {}".format(len(old_vartables),dtimes[1:])
    if debug:
        print debug_head+"writing final tables"
    if multiproc:
        lock.acquire()
    if not arguments.noalign:
        for table in alntables:
            IO.write_table(table,"alignments")
    asdf = 0
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

def gather_models(unp,isoform,models,debug):
    '''Gather the appropriate model files
       If there is a isoform that matches, use those models
       If not, use all models
       Returns a list of models in format:
       [isoform,identity,filename]'''
    debug_head = "DEBUG: workers: gather_models: "
    if debug:
        print debug_head+"Gathering models for {} {}".format(unp,isoform)
    current_models = list()
    if unp in models:   
        if debug:
            print debug_head+"{} has models, selecting".format(unp)
        for model in models[unp]:
            if model[0]==isoform:
                if debug:
                    print debug_head+"model matched isoform {}, adding".format(isoform)
                current_models.append(model)
            #In case the -1 has been inferred as no isoform designation
            #There shouldn't be a -1
            elif (isoform.endswith("-1") and 
                 len(
                 [x[0] for x in models[unp] if x[0].endswith("-1")]
                 )==0 and
                 "-" not in model[0]):
                if debug:
                    print debug_head+"isoform {} and model {}, missing -1, adding".format(isoform,model[0])     
                current_models.append(model)                 
        if len(current_models) == 0: # If no isoform was found, just use them all
            if debug:
                print debug_head+"No isoform found for {}, adding all {} models".format(isoform,len(models[unp]))
            current_models = models[unp]
    if len(current_models)==0:
        if debug:
            print debug_head+"No models found for {}".format(unp)
        raise WorkerException(
              "swissmodel_lookup","no models entry for {}".format(unp))
    else:
        if debug:
            print debug_head+"Returning models: {}".format(current_models)
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
        raise WorkerException(
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
        residues['template_identity'] = 100.0
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
#    all_residues.to_csv(open("tempstr","w"))
    return all_residues

def expand_variants(variants,trans_seq,debug):
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
    debug_head="DEBUG: workers: expand_variants"
    if debug:
        print debug_head+"Expanding {} variants".format(len(variants))
    failures = list()
    newvariants = list()
    reason = "Failed to expand"
    for variant in variants:
        position,refaa,altaa,annotation,varcode = variant
        # Don't expand missense        
        if "MISSENSE" in annotation:
            if debug:
                print debug_head+"missense, skipping: {}".format(variant)
            newvariants.append(variant)
            continue
        # Expand in-frame deletions
        elif "-" in position and "INFRAME_DELETION" in annotation:
            if debug:
                print debug_head+"in-frame variant: {}".format(variant)
            newref = refaa
            if len(altaa)>1:
                if debug:
                    print debug_head+"failed expand, altaa>1"
                write.failures(variant,reason)
                continue
            try:                
                pos1,pos2 = [int(x) for x in variant[0].split("-")]
            except IndexError,ValueError:
                if debug:
                    print debug_head+"failed expand, bad position {}".format(variant[0])
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
            if debug:
                print "adding altaa as - for residues {} thru {}".format(pos1,pos2)
            for i,x in enumerate(range(pos1,pos2+1)):   
                newvariants.append([i,newref,"-",annotation,varcode])
        # Single residue deletions
        elif "INFRAME_DELETION" in annotation and "STOP" not in annotation:
            if debug:
                print debug_head+"inframe_del, nonstop: {}".format(variant)
                print debug_head+"adding - altaa to position {}".format(position)                
            newvariants.append([int(position),refaa,"-",annotation,varcode])
        # Everything else (frameshift, early stop) are from pos to end of transcript
        else:
            startpos = int(position.split("-")[0])
            if debug:
                print debug_head+"misc: {}".format(variant)
                print "adding altaa X for residues {} thru {}".format(startpos,len(trans_seq))            
            for pos in range(startpos-1,len(trans_seq)):
                newvariants.append([pos+1,refaa,"X",annotation,varcode])
    if debug:
        print "Final variant set contains {} positions".format(len(newvariants))
    return newvariants                            
                                                                                         
