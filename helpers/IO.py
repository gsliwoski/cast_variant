from os import path
from exceptions import *
import gzip
import sys
import pickle
import xml.etree.ElementTree as ET
import traceback
from AA import *
import pandas as pd

##### INITIALIZATION FUNCTIONS #####

# Define what the sequence files are
SEQ_PATH = "./sequences/" # Where pickles are
SWISS_SEQ = "swissmodel.pickle" # Model mapping
UNP_SEQ = "unp.pickle" # Uniprot sequences
TRANS_SEQ = "trans.pickle" # Transcript sequences
SIFTS_SEQ = "sifts.pickle" # uniprot - pdb mapping
SIFTS_PATH = "/dors/capra_lab/data_clean/sifts/2018-08-08/xml/" # Where sifts alignments are

# Ensure that the required sequence pickles are available
def check_seqs(nomodel,nopdb):
    for x in [SWISS_SEQ, UNP_SEQ, TRANS_SEQ, SIFTS_SEQ]:
        if (x=="SWISS_SEQ" and not nomodel) or\
           (x!="SWISS_SEQ" and not nopdb):\
            assert path.isfile(SEQ_PATH+x),\
            "{}{} not found, run {}pickle_sequences.py".format(
                                                        SEQ_PATH,
                                                        x,
                                                        SEQ_PATH)
    if not nopdb:
        assert path.isdir(SIFTS_PATH), "{} not found. "\
                                   "sifts dataset required for PDB".format(
                                                                    SIFTS_PATH)
                                                                               

# These are the names of the all files that will be written to
outfiles = dict()
def define_output(varfilename):
    global outfiles
    basefile = path.splitext(path.basename(varfilename))[0]
    outfiles = {
        'skipped': basefile+".skipped",
        'alignments': basefile+".alignments",
        'variants': basefile+".variants",
        'failures': basefile+".failures",
        'completed': basefile+".completed"}

# These are the header and column formats for output
# For models, the unp columns will be NA and therefore
# filtering out models will require uniprot_position=="NA" filter
CONV_HEADER = ["transcript","uniprot","isoform",
               "transcript_identity","transcript_position",
               "uniprot_position", "structure_position", "icode",
               "transcript_aa", "uniprot_aa", "structure_aa",
               "ref_aa", "alt_aa", "structure_identity",
               "structure", "chain","varcode"]

ALN_HEADER = ["transcript","uniprot","isoform",
              "transcript_identity","transcript_position",
              "uniprot_position","structure_position",
              "transcript_aa","uniprot_aa","structure_aa",
              "structure_identity","structure","chain"]
                             

##### VARIANT PARSING #####

# Read in the variant file
def process_variants(varfilename,expand,continue_flag):
    # Try to open the file
    try:
        infile = open(varfilename)
    except IOError as e:
        sys.exit("Unable to open {}: {}".format(varfilename,e))
    # Try to parse it
    try:
        variants = parse_varfile(infile,expand,continue_flag)
    except ParseException as e:
        print "Critical failure:"
        sys.exit(e.fullmsg)
    return variants

# Process a single variant line
def process_variant(variant,expand):
    '''
    If expand is set:
    Expands usable non-missense variants to list all residues affected.
    Deletions appear as - for all residues affected
    Frameshift appear as X from start of frameshift to end of transcript
    Early stop appear as X from position of stop to end of transcript
    Takes the list of columns of the current variant line and returns
    [position, ref, alt,annotation(for expanding if necessary)]
    '''
    if variant[-2].upper() in [""," ","-", ".","?","NA","NONE","UNASSIGNED"]:
        return ["no uniprot assigned"]
    annotation = variant[0].upper().replace(" ","_")
    # Anything can pair with NMD so skip those immediately
    if "NMD" in annotation or "NONSENSE_MEDIATED_DECAY" in annotation:
        return ["NMD variant"]
    try:
        position = variant[2]
        refaa,altaa = variant[3].split("/")
        # TODO: Remove this line if one day desire synonymous mapping 
        if refaa==altaa:
            return ["refaa equals altaa"]
    # Check the refaa and altaa as we go        
    # Missense is always kept
        if "missense" in variant[0]:
            # Your standard missense nothing needs to change
            if len(refaa)==1 and len(altaa)==1:
                position=int(position)
            # Sometimes VEP will give AA/AB or AA/BA as missense
            # with position 1-2 instead of 1 or 2
            elif len(refaa)==2 and len(altaa)==2:
                position1,position2 = position.split("-")
                #TODO: Adjust these checks if one day desire synonymous mapping
                if refaa[0]==altaa[0]:
                    refaa,altaa = refaa[1],altaa[1]
                    position = int(position2)
                #TODO: Adjust this part if one day desire DNP or TNP
                else:
                    refaa,altaa = refaa[0],altaa[0]
                    position = int(position1)
            # Anything else is a bad missense and outputs warning and skips
            else:
                raise ValueError
        # Stop losses and changes to start can't be mapped
        elif "STOP_LOST" in annotation or "STOPLOST" in annotation:
            return ["stop loss variant"]
        elif "START" in annotation:
            return ["variant affects start"]
        # Only keep the rest if expand is set
        elif not expand:
            return ["expand not set"]
        # Make sure the inframe deletion can be used
        elif "INFRAME_DELETION" in annotation:
            if len(altaa)>1:
                raise ValueError
    except ValueError:
        raise ParseWarning(
            "variant file","bad variant column: {}".format(variant[3]))               
    except ParseWarning as e:
        print e.fullmsg
        return [e.fullmsg]
    return [position,refaa,altaa]
    
# Iterate through variant lines    
def parse_varfile(varfile,expand,continue_flag):
    '''
    Parses variant file. Checks for format adherence.
    Returns a dictionary of variants with entries:
    Key = Transcript ID (for ENST) or Protein ID (for NP)
    Entry = {uniprot(str), isoform(str), variants(list_of_variants)}
        variants are [protein_position, refaa, altaa, misc]
    Function takes instream of variant file and dict variable to store vars
    ---
    Required format = tab delim (NO HEADER; SKIPS ANYTHING AFTER #)
    Required columns
    col 1 = Consequence (missense, frameshift, deletion, etc)
    col 2 = Transcript ID
    col 3 = Transcript Position
    col 4 = Mutation as refaa/altaa
    col 5 = protein ID
    col 6 = Uniprot ID
    col 7 = Uniprot-Isoform
    col 8 = optional unique variant code
            useful when handling non-missense variants
    Returns a dict with entries:
    KEY = variant identifier (transcript or protein)
    ENTRY = [[uniprot, isoform],[list of variants]]
    '''
    completed = list()
    if continue_flag and path.isfile(outfiles['completed']):
        with open(outfiles['completed']) as infile:
            completed = [x.strip() for x in infile]
    skipped = list()
    variants = dict()
    count = 0
    for line in varfile:
        # Skip blank lines
        # Remove comments
        # Report incomplete lines
        # Skip bad variants
        # The following are "good" variants:
        #   All good variants must have a uniprot assigned
        #   missense in form R/A
        #   missense in the form RR/AA where either the first or second
        #       R,A pair is the same and the other is different.
        #   in-frame deletions if expand is set
        #   frameshift if expand is set
        #   early stop if expand is set
        #
        # All others are written to the skipped file
        line = line.strip().split("#")[0]
        if line=="": continue
        line = line.strip().split("\t")
        if len(line)==0: continue
        annotation = ["unknown reason"]
        try:
            if len(line)<7:
                raise ParseWarning("variant_file","incomplete line: {}".format(
                                               "\t".join(line)))
            consequence,trans,pos,var,protein,unp,iso = line[:7]
            try:
                varcode = line[7]
            except IndexError:
                varcode = "_".join([consequence,
                                   trans,
                                   pos,
                                   var,
                                   protein,
                                   unp])

            
            if continue_flag and varcode in completed: continue                               
            current_var = process_variant(line,expand)           
            #Use ENST identifier if its ensembl otherwise protein identifier
            identifier = line[1] if line[1].startswith("ENST") else line[4]
        except ParseWarning as e:
            print e.fullmsg
            continue      
        if len(current_var)<3:
            skipped.append("\t".join(line+current_var))
        elif identifier in variants:
            count += 1
            variants[identifier][1].append(current_var+[varcode])
        else:
            count += 1
            variants[identifier] = [[unp,iso],[current_var+[varcode]]]
    print variants
    print "{} variants processed".format(count)
    if len(skipped)>0:
        with open(outfiles['skipped'],'w') as outfile:
            outfile.write("\n".join(skipped))
            outfile.write("\n")
        print "{} variants skipped".format(len(skipped))
    return variants

##### PICKLE LOADING #####

def open_seqfile(filename):
    filename = SEQ_PATH+filename
    try:
        infile = open(filename)
    except IOError:
        sys.exit("Failed to open {}".format(filename))
    return infile
            
# Load the transcript sequences
def load_transcripts():
    infile = open_seqfile(TRANS_SEQ)
    trans = pickle.load(infile)
    print "loaded {} transcript sequences".format(len(trans))
    return trans
    
# Load the sequences required for PDB casting
def load_uniprot():
    infile = open_seqfile(UNP_SEQ)
    unp = pickle.load(infile)
    print "loaded {} uniprots".format(len(unp))
    return unp

def load_sifts():
    infile = open_seqfile(SIFTS_SEQ)
    sifts = pickle.load(infile)
    print "loaded {} sifts uniprots".format(len(sifts))
    return sifts

# Load model mappings
def load_models(source):
    if source=="swissmodel":
        infile = open_seqfile(SWISS_SEQ)
        models = pickle.load(infile)
        print "loaded {} swissmodel uniprots".format(len(models))
    else:
        sys.exit("Critical: unrecognized model-type {}".format(source))
    return models

##### SIFTS LOADING #####

def parse_sifts(pdb):
    '''
    Parses a sifts xml for a given pdbid
    Takes a pdb that is a list of [structid,chain]
    Returns a dict of all the residues in the pdb
    '''
#    print pdb
    structure,chain = pdb
    residues = {'structure': structure,
                'chain': list(),
                'structure_position': list(),
                'icode': list(),
                'structure_aa': list(),
                'uniprot_position': list(),
                'secondary_structure': list(),
                'uniprot': list()}
                
    total_res = 0
    id_res = 0
    pdbfile = SIFTS_PATH+structure+".xml.gz"
    try:
        infile = gzip.open(pdbfile)
    except:
        raise AlignException("parse_sifts", "failed to open {}".format(filename))
    try:
        # SIFTS files should be XML so parse into table
        root = ET.parse(infile).getroot()

        # Remove the namespace from tags
        for elem in root.getiterator():
            if not hasattr(elem.tag,'find'): continue
            i = elem.tag.find('}')
            if i >= 0:
                elem.tag = elem.tag[i+1:]
        
        # Iterate over all PDB chains
        for chain in root.findall("entity"):
            if chain.get("type") == "protein":
                # Iterate over annotation segments
                for s in chain.getchildren():
                    # Iterate over segment residues
                    for residue in s.find("listResidue"):
                        # Parse residue annotations
                        res = {'chain': "",
                               'structure_position': None,
                               'icode': "",
                               'structure_aa': "",
                               'uniprot_position':None,
                               'secondary_structure':"",
                               'uniprot':""}
                        nullres = False
                        for db in residue.findall("crossRefDb"):
                            if db.get("dbSource") == "UniProt":
                                res['uniprot'] = db.get("dbAccessionId")                                                                    
                            if db.get("dbSource") == "PDB":
                                res['chain'] = db.get("dbChainId")
                                res['structure_position'] = db.get("dbResNum")
                                if res["structure_position"] == "null":
                                    nullres = True
                                try: # Check for insertion code
                                    int(res["structure_position"][-1])
                                except ValueError:
                                    res["icode"] = res["structure_position"][-1]
                                    res["structure_position"] = res["structure_position"][:-1]
                                try: # Try converting to single letter AA
                                    res["structure_aa"] = AA[db.get("dbResName")]
                                except KeyError:
                                    try: # Try converting noncanonical AA to canonical
                                        res["structure_aa"] = AA[MODRES[db.get("dbResName")]]
                                    except KeyError:
                                        res["structure_aa"] = "X" # Can't resolve AA letter
                            elif db.get("dbSource") == "UniProt":
                                res["uniprot_position"] = db.get("dbResNum")
                                res["uniprot_aa"] = db.get("dbResName")
                        for rd in residue.findall("residueDetail"):
                            if rd.get("property") == "codeSecondaryStructure":
                                res['secondary_structure'] = rd.text
                        if not nullres and \
                               res['structure_position'] is not None and \
                               res['uniprot_position'] is not None:
                            res['uniprot_position'] = int(res['uniprot_position'])
                            total_res += 1
                            if res["uniprot_aa"] == res["structure_aa"]:
                                id_res += 1
                            for x in ['chain',
                                      'structure_position',
                                      'icode',
                                      'structure_aa',
                                      'uniprot_position',
                                      'secondary_structure',
                                      'uniprot']:
                                residues[x].append(res[x])
    except Exception:       
        raise AlignException("sifts parse",traceback.format_exc())
    finally:
        infile.close()
    if len(residues["chain"]) == 0:
        print "Warning, no usable residues in {}".format(pdbfile)
    residues["structure_identity"] = float(id_res)/total_res
    return residues
            
#### OUTPUT WRITING #####

def write_failures(source,msg,lock=None):
    '''
    Writes any failures encountered during variant alignment
    and tries to assign reason with msg
    recieves lock that is not none when multiproc
    '''
    outmsg = "{}\t{}\n".format("\t".join(source),msg)
    if lock:
        lock.acquire()
    with open(outfiles['failures'],'a+') as outfile:
        outfile.write(outmsg)
    if lock:
        lock.release()

def generate_alignment_table(uniprot,structures):
    '''
    Generate the alignment dataframe
    Translates sequences into lists
    Merges the transcript and structure alignments
    '''
    uniprot['uniprot_aa'] = [x for x in uniprot['uniprot_aa']]
    uniprot['transcript_aa'] = [x for x in uniprot['transcript_aa']]
    uniprot.pop('fasta')
    alignment_table = pd.DataFrame.from_dict(uniprot)
#    alignment_table.to_csv(open("tempaln",'a'))
    alignment_table = alignment_table.merge(structures,
                                            how="left",
                                            on=['uniprot_position'])
#    structures.to_csv(open("tempstruct","a"))
#    alignment_table.to_csv(open("tempaln2","a"))
    alignment_table = alignment_table.round({
                                'transcript_identity':1,
                                'structure_identity':1})
    
    alignment_table.sort_values(by=["structure","chain","transcript_position",
                                    "structure_position","icode"],inplace=True)
    return alignment_table

def generate_variant_table(variant_list):
    '''
    Generates a variant dataframe
    from a variant list
    '''
    variant_dict = {'transcript_position': [x[0] for x in variant_list],
                    'ref_aa': [x[1] for x in variant_list],
                    'alt_aa': [x[2] for x in variant_list],
                    'varcode': [x[3] for x in variant_list]}
    variant_df = pd.DataFrame.from_dict(variant_dict)
#    variant_df.to_csv(open("tempvar","a"))
    return variant_df
            
def generate_conversion_table(alignment_table,variant_table):
    '''
    Intersects alignment table with variant table
    '''
    conv_table = alignment_table.merge(variant_table,
                                       how='inner',
                                       on=['transcript_position'])
    # Filter any variants that were not on a structure
    conv_table = conv_table[pd.notnull(conv_table["structure"])]
    # Group structures and sort residues
    conv_table.sort_values(by=["structure","chain","transcript_position",
                               "structure_position","icode","ref_aa"],inplace=True)


    return conv_table

def write_output(align_out,conv_out):
    aout = outfiles['alignments']
    if path.isfile(aout):
        align_out.to_csv(
            open(aout,'a'),
            header=False,
            index=False,
            sep="\t",
            columns=ALN_HEADER)
    else:
        align_out.to_csv(
            open(aout,'w'),
            index=False,
            sep="\t",
            columns=ALN_HEADER)
    cout = outfiles['variants']
    if path.isfile(cout):
        conv_out.to_csv(
            open(cout,'a'),
            header=False,
            index=False,
            sep="\t",
            columns=CONV_HEADER)
    else:
        conv_out.to_csv(
            open(cout,'w'),
            index=False,
            sep="\t",
            columns=CONV_HEADER)

def write_complete(varcodes):
    with open(outfiles["completed"],"a+") as outfile:
        outfile.write("\n".join(varcodes)+"\n")                   
