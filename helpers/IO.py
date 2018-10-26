from os import path
from exceptions import *
import gzip
import sys
import pickle
import xml.etree.ElementTree as ET
import traceback
from AA import *
import pandas as pd
from Bio import PDB
import re
from numpy import where as npwhere

##### INITIALIZATION FUNCTIONS #####

# Define what the dataset files are
SEQ_PATH = "./sequences/" # Where required datasets are
SWISS_SEQ = "swissmodel.pickle" # Model mapping
UNP_SEQ = "unp.pickle" # Uniprot sequences
UNP_CANONICAL = "uniprot_canonical_isoforms.tab" # Canonical uniprot ACs
UNP_MAP = "uniprot_sec2prim_ac.txt" # Uniprot secondary to primary AC map
TRANS_SEQ = "trans.pickle" # Transcript sequences
SIFTS_SEQ = "sifts.pickle" # uniprot - pdb mapping
SIFTS_PATH = "/dors/capra_lab/data_clean/sifts/2018-08-08/xml/" # Where sifts alignments are
DSSP = '/dors/capra_slab/bin/dssp' # DSSP application (required for some descriptors)

# Ensure that the required datasets are available
def check_seqs(nomodel,nopdb):
    for x in [SWISS_SEQ,
              UNP_SEQ,
              UNP_CANONICAL,
              UNP_MAP,
              TRANS_SEQ,
              SIFTS_SEQ]:
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
                                                                               
def check_applications(name):
    name_map = {'DSSP':DSSP}
    assert path.isfile(name_map[name]),\
        "{} not found, remove {} dependent descriptors".format(
            name_map[name],
            name)

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

##### VARIANT PARSING #####

# Read in the variant file
def load_variants(args):
    varfilename = args.variants
    # If a raw VEP was given, select one consequence per coding variant
    # and then parse
    if args.vep:
        vepoutfile = path.splitext(path.basename(varfilename))[0]+".vep_selected"
        if path.isfile(vepoutfile):
            print "Found potential selected file {} already created, overwrite it?".format(vepoutfile)
        try:
            a = raw_input("(y/n/q) ").lower()[0]
        except SyntaxError:
            a = "NA"            
        while a not in ["y","n","q"]:
            try:
                a = raw_input("(y/n/q) ").lower()[0]
            except SyntaxError:
                a = "NA"
        if a=="q":
            sys.exit("Canceled")
        if a=="y":                        
            fullvep = parse_vepfile(varfilename)
            unique_vep = process_vep(fullvep)
            try:
                if len(unique_vep.index) == 0:
                    raise ParseException("VEP file","No variants selected from VEP file")
            except ParseException as e:
                sys.exit(e.fullmsg)  
            # Write selected vars to file and parse them as normal
            unique_vep.to_csv(open(vepoutfile,"w"),
                              sep="\t",
                              columns=["Consequence",
                                       "Transcript",
                                       "Protein_position",
                                       "Amino_acids",
                                       "Protein",
                                       "Uniprot",
                                       "Isoform",
                                       "Varcode"],
                              header=False,
                              index=False)                                   
            print "{} Selected vars written to standard varfile {}".format(
                        len(unique_vep.index),
                        vepoutfile)
        variants = parse_varfile(open(vepoutfile),args)        

    else: # Otherwise, parse it directly
        #Try to open the file
        try:
            infile = open(varfilename)
        except IOError as e:
            sys.exit("Unable to open {}: {}".format(varfilename,e))
        # Try to parse it
        try:
            variants = parse_varfile(infile,args)
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
    [position, ref, alt,annotation(for expanding if necessary),varcode]
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
    return [position,refaa,altaa,annotation]
    
# Iterate through variant lines    
def parse_varfile(varfile,args):
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
    expand = args.expand
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

            
#            if continue_flag and "{}_{}".format(varcode, in completed: continue                               
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
    print "{} variants processed".format(count)
    if len(skipped)>0:
        with open(outfiles['skipped'],'w') as outfile:
            outfile.write("\n".join(skipped))
            outfile.write("\n")
        print "{} variants skipped".format(len(skipped))
    return variants

completed = list()
def filter_complete(variants):
    global completed
    keptvars = dict()
    if path.isfile(outfiles['completed']) and len(completed)==0:
        with open(outfiles['completed']) as infile:
            completed = [x.strip() for x in infile]
    if len(completed)==0:
        return variants
    for trans in variants:
        filtered = list()
        for var in variants[trans][-1]:
            varcode = var[-1]
            if varcode in completed:
                continue
            else:
                filtered.append(var)
        if len(filtered)>0:
            keptvars[trans] = [variants[trans][0],filtered]                
    return keptvars

# Parse raw VEP output file

def parse_vepfile(vepfile):
    '''
    Reads in a raw VEP file
    Takes a filename and first checks if necessary columns are there
    Then, reads it in as a pandas DataFrame, prepares the columns,
    and filters for coding variants
    Returns prepared dataframe
    '''
    # Check to make sure the necessary columns are there
    with open(vepfile) as infile:
        header = infile.readline().strip().split()
    required_cols = ["#Uploaded_variation",
                     "Location",
                     "Allele",
                     "Consequence",
                     "SYMBOL",
                     "Gene",
                     "Feature",
                     "BIOTYPE",
                     "Protein_position",
                     "Amino_acids",
                     "Codons",
                     "ENSP",
                     "SWISSPROT",
                     "TREMBL"]
    for x in required_cols:
        if x not in header:
            raise ParseException(
                "Raw VEP file",
                "required column {} missing from header".format(x))
    df = pd.read_csv(open(vepfile),sep="\t")[required_cols]
    df.rename(columns={"#Uploaded_variation":"Varcode",
                       "Feature":"Transcript",
                       "ENSP":"Protein"},inplace=True)
    df = df[df["Amino_acids"].str.len()>1]
    print "{} unique coding variants loaded".format(df["Varcode"].nunique())
    return df

    
def process_vep(vep):
    '''
    Takes the full VEP dataframe and selects a single
    consequence per variant, order of preference:
    1) Select canonical ENST when possible
    2) Select canonical RefSeq when possible
    3) Select ENST with unassigned canonical
    4) Select Refseq with unassigned canonical
    5) Select ENST specifically not canonical
    6) Select Refseq specifically not canonical
    7) Select XM entries left over
    8) Select entries with uniprot assignments only from VEP
    Anything without a Uniprot assignment is filtered
    In the event of a tie, take highest transcript
    ID if AA and pos are the same, otherwise take highest pos
    returns a list of unique variants in the form of
    the typical variant input file

    Note: In rare cases, the same variant will hit multiple
    primary uniprot entries. In these cases, uniprot has identified
    sufficiently different isoforms to call them different proteins
    Therefore, repeats across different primary uniprot entries may occur.
    
    A minor (major?) drawback to the current approach:
    If there is a variant which affects 2 different uniprot entries and in one
    entry it hits the canonical transcript, while in the other it hits a 
    non-canonical isoform, since the variant is filtered out before getting
    to the non-canonical isoforms, both wont be captured.
    One example is 1:g.156842168C>T which hits INSRR (one transcript) and NTRK1
    isoform 3. Only INSRR is kept since it has 1 transcript and isoform 3 is not
    canonical.

    Another drawback that might need to be addressed:
    As it stands, if it has multiple potentials in the same step,
    it selects one based on max protein position and if still multiple
    max transcript ID. The problem with this is that you can have 2 variants
    assigned different transcripts even though it would probably be easier and better
    if they both had the same one. For example, vars A and B don't hit canonical.
    Var A is position 10, Var B is position 100. There is isoform name ENST1 with 
    100 residues and ENST2 with 98 residues. The missing residues are at position
    11-12
    Var A gets ENST2 since it's 10 in both but ENST2>ENST1
    Var B gets ENST1 since it's at 100 in ENST1 and 98 in ENST2
    It's probably preferable to put them both in ENST1.
    However, the current datasets (and uniprot_sprot_human.tab) don't have transcript
    length. So, the best solution (to assign a transcript to a uniprot at a time before
    selecting with var, and resolving overlaps by using longest transcript can't be done
    without pulling in another dataset. I'm not sure how often this happens and if it
    warrants the extra time required for this better solution.
    One potential solution that would also solve another drawback is to read in the sequence
    sets to get the transcript lengths. This would allow cases (no idea how common or rare)
    where the variant hits a transcript in the fasta set but is assigned a different transcript
    that is not in the sequences set but is, for example, canonical according to uniprot. In this
    case, a transcript present in the sequence set could be preferred.
    '''        
    # canonical identifies the canonical transcript used by uniprot
    canonical = load_canonical()
    # sec2prime allows filtering out any secondary uniprot ACs
    sec2prime = load_sec2prime()
    
    # Attach the uniprot names
    vep = vep.merge(canonical,how="left",on=["Transcript","Protein"]) 

    # Filter out anything that doesn't have a uniprot name
    # Keep those without a uniprot name for use in case they have one
    # assigned by VEP
    vep_nullunp = vep[vep.Uniprot.isnull()]
    vep = vep[vep.Uniprot.notnull()]
    nrows = len(vep.index)
    try:
        if len(vep.index) == 0:
            raise ParseException("VEP file","Failed to extract variants from VEP file")
    except ParseException as e:
        sys.exit(e.fullmsg)  

#    print "{} unique coding variants were mapped to a uniprot".format(
#                                                        vep.Varcode.nunique())

    # Filter any secondary uniprot AC's
    secondary = [x[0] for x in load_sec2prime()]
    vep = vep[~vep.Uniprot.isin(secondary)]
    
    # Now select 1 consequence per variant
    # Note: There may be mulitple consequences if they map to
    # different uniprots (different proteins, same gene)

    def addvars(vep,vep_final,group): #append to final set and filter from initial
        vep_final = pd.concat([vep_final,
                group.apply(lambda x: x.sort_values(["quickpos","Transcript"],
                            ascending=False).head(1))])
        return vep[~vep.Varcode.isin(vep_final.Varcode)],vep_final                                                    

    vep_final = None
    
    #Need to add a quickposition column that is the int of the first
    # position for finding max position instance for repeats
    vep["quickpos"] = vep.Protein_position.str.extract('(\d+)',expand=False).astype(int)
    vep_nullunp["quickpos"] = vep_nullunp.Protein_position.str.extract('(\d+)',expand=False).astype(int)

    print "selecting unique variants"
    # Loops through steps 1-6:
    # Canonical = YES: ENST then NM then XM
    # Canonical = Unassigned: ENST then NM then XM
    # Canonical = NO: ENST then NM then XM
    
    # Typically, when a transcript is unassigned it means it is the only
    # one and therefore no isoform is designated. However, there is an edge
    # case where it is unassigned even though there are assignments
    # So, need to deal with these after and remove them for now
    vep_unassigned = vep[(vep.Canonical=="Unassigned")\
                          & (pd.to_numeric(vep.nIsoforms)>0)]
    vep = vep[~((vep.Canonical=="Unassigned")\
                          & (pd.to_numeric(vep.nIsoforms)>0))]   
    print "vep"
    print vep[vep.Varcode=="1:g.155187269G>A"]
    print "unassigned"
    print vep_unassigned[vep_unassigned.Varcode=="1:g.155187269G>A"]    
    print "nullunp"
    print vep_nullunp[vep_nullunp.Varcode=="1:g.155187269G>A"]    
    
    for outer in ["YES","Unassigned","NO"]:
        for inner in ["ENST","NM","XM"]:
            if len(vep.index)==0: break
            print "{}: {}".format(outer,inner)
            group = vep[(vep.Canonical==outer)\
                    & (vep.Transcript.str.startswith(inner))]\
                    .groupby(["Varcode","Uniprot"])
            if vep_final is None:
                vep_final = group.apply(lambda x: x.sort_values(["quickpos","Transcript"],
                                                        ascending=False).head(1))
                vep = vep[~vep.Varcode.isin(vep_final.Varcode)]
            else:
                vep,vep_final = addvars(vep,vep_final,group)   
    # Deal with any of potential cases where an unassigned transcript is only hit
    # in a uniprot with multiple isoforms
    vep_unassigned = vep_unassigned[~vep_unassigned.Varcode.isin(vep_final.Varcode)]
    for inner in ["ENST","NM","XM"]:
        if len(vep_unassigned.index)==0: break
        group = vep_unassigned[vep_unassigned.Transcript.str.startswith(inner)]\
                                .groupby(["Varcode","Uniprot"])
        current_vars = group.apply(lambda x: x.sort_values(["quickpos","Transcript"],
                                             ascending=False).head(1))
        vep_final = pd.concat([vep_final,current_vars])
        vep_unassigned = vep_unassigned[~vep_unassigned.Varcode.isin(vep_final.Varcode)]        
    
    # Step 8: select any remaining variants that could not
    # be assigned a uniprot entry based on the read in uniprot data.
    # Only consider instances where VEP provided a uniprot entry.
    # In some cases, the VEP uniprot may not match the current uniprot
    # which is why only those in the sequence datasets provided to this
    # program were considered first and anything left over the uniprot is 
    # taken from VEP
    vep_nullunp = vep_nullunp[~vep_nullunp.Varcode.isin(vep_final.Varcode)]
    for outer in ["SWISSPROT","TREMBL"]:
        for inner in ["ENST","NM","XM"]:
            if len(vep_nullunp.index)==0: break
            print "{}: {}".format(outer,inner)
            group = vep_nullunp[(vep_nullunp[outer]!="-")\
                                & (vep_nullunp.Transcript.str.startswith(inner))]\
                                .groupby(["Varcode",outer])
            current_vars = group.apply(lambda x: x.sort_values(["quickpos","Transcript"],
                                                        ascending=False).head(1))
            vep_final = pd.concat([vep_final,current_vars])
            vep_final["Uniprot"] = npwhere(vep_final.Uniprot.isnull(),
                                           vep_final[outer],
                                           vep_final.Uniprot)
            vep_final["Isoform"] = npwhere(vep_final.Isoform.isnull(),
                                           vep_final[outer],
                                           vep_final.Isoform)                                          
            vep_nullunp = vep_nullunp[~vep_nullunp.Varcode.isin(vep_final.Varcode)]

    print "{} unique variants with uniprot assignment selected".format(
                                            vep_final.Varcode.nunique())
    vep = pd.concat([vep,vep_nullunp,vep_unassigned])
    if len(vep.index)>0:
        print "The following {} unique coding variants were unable to be assigned:".format(
            vep.Varcode.nunique())
        print ",".join(set(vep.Varcode))
    return vep_final    
      
##### DATASET LOADING #####

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

# Load sifts uniprot-PDB alignments
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

# Load table of canonical uniprot isoforms
def load_canonical():
    df = pd.read_csv(open(SEQ_PATH+UNP_CANONICAL),sep="\t")
    df.drop_duplicates(inplace=True)
    return df

# Load secondary to primary AC mapping
def load_sec2prime():
    with open(SEQ_PATH+UNP_MAP) as infile:
        sec2prime = [x.strip().split() for x in infile]    
    return sec2prime
    
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

def write_table(df, destination):
    '''
    Writes pandas table to destination
    takes a [[df,headers],destination]
    '''
    if destination not in outfiles:
        return
    outfile = outfiles[destination]
    table,headers = df
    
    if path.isfile(outfile):
        table.to_csv(
            open(outfile,'a'),
            header=False,
            index=False,
            sep="\t",
            columns=headers)
    else:
        table.to_csv(
            open(outfile,"w"),
            index=False,
            sep="\t",
            columns=headers)

def write_complete(varcodes):
    '''
    Record completion, takes list of variant codes
    and a string destination (PDB, SWISSMODEL, etc)
    '''
    with open(outfiles["completed"],"a+") as outfile:
        outfile.write("\n".join(x for x in varcodes)+"\n")

##### MODELS #####

parser = PDB.PDBParser(PERMISSIVE=1)

def load_model(modelfile):
    '''Loads the model file and generates a fasta
       with dashes for any skipped residue #'s
       Takes a filename and returns a dict with
       filename: model filename
       fasta: fastaseq (of chain, with - for skipped res #)
       chain: chainID (only uses first chain if multiple)
       icodes: sequence of icodes (' ' for none)
       resnums: list of residue numbers'''
    structure = parser.get_structure("Model",modelfile)
    chainID,complex_state = model_information(modelfile)
    chain = structure[0][chainID]
    resnum = 0
    fasta = ""
    icodes = ""
    resnums = list()
    for residue in chain:
        #Skip hetatoms
        if not PDB.is_aa(residue): continue
        resid = residue.get_id()
        #Fill in gaps with -
        while resid[1] > resnum+1:
            fasta += '-'
            resnum+=1
        hetflag,resnum,icode = resid
        aa = AA[residue.get_resname()]
        fasta += aa
        icodes += icode
        resnums.append(resnum)
    return {'filename':modelfile,
            'fasta':fasta,
            'chain':chainID,
            'icodes':icodes,
            'resnums':resnums,
            'complex_state':complex_state}

def model_information(modelfile):
    '''Gets specific model information.
       Regarding complex status and current chain.
       Relevant for non-monomers'''
    #Right now, I'm getting the current chain based on the filename
    #However, I think in REMARK 3 the MMCIF also corresponds to chain
    chain = re.search("[.][\w][_][\w]+\.pdb$",modelfile).group(0)[1:2]
    #Get the complex state of the model
    with open(modelfile) as infile:
        complex_state = 'monomer'
        for line in infile:
            if not line.strip().startswith("REMARK   3  OSTAT"): continue
            complex_state = line.strip().split()[-1]
    return chain,complex_state             
