from os import path,remove
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
from numpy import vectorize

# Dictionary of paths and applications
CONFIG = {"SEQ_PATH":"./",
          "SWISS_SEQ":"",
          "UNP_SEQ":"",
          "UNP_CANONICAL":"",
          "UNP_MAP":"",
          "TRANS_SEQ":"",
          "SIFTS_SEQ":"",
          "SIFTS_PATH":"",
          "DSSP":"",
          "PDB_PATH":"",
          "SWISS_PATH":"",
          "ARTIFACTS_FILE":""}
          
# Add each descriptor to applications list
# Which is set in load_config() below
# If descriptor requires no application, then set to None
# Used to check for presence before beginning
DESCRIPTOR_APPLICATIONS = dict()

# Set of artifacts used if artifact filtering is turned on
ARTIFACTS = None

##### INITIALIZATION FUNCTIONS #####

def load_config(arguments):
    debughead = "DEBUG: IO: load_config: "
    global CONFIG
    global DESCRIPTOR_APPLICATIONS
    configfile = "config.sys"
    required = set(["TRANS_SEQ","UNP_SEQ"])
    # Collect all the required settings to ensure they are found
    if not arguments.nopdb:
        required.update(["SIFTS_SEQ","SIFTS_PATH"])
        if 'ligand' in arguments.descriptors or\
           'nucleotide' in arguments.descriptors or\
           'peptide' in arguments.descriptors or\
           'dssp' in arguments:
            required.update(["PDB_PATH"])
    if not arguments.nomodel:
        required.update(["SWISS_SEQ","SWISS_PATH"])
    if arguments.vep:
        required.update(["UNP_CANONICAL","UNP_MAP"])
    if 'dssp' in arguments.descriptors:
        required.update(["DSSP"])
    if 'unp' in arguments.descriptors:
        required.update(["UNP_FEATURES"])
    if "artifacts" not in arguments.descriptors:
        required.update(["ARTIFACTS_FILE"])                                         
    with open(configfile) as infile: #TODO: Make as dict
        for rawline in infile:
            line = rawline.strip().split("#")[0].split()
            if len(line)<2: continue
            curset,curval = line[:2]
            CONFIG[curset] = curval
    required = required - set(CONFIG.keys())
    if len(required)>0:
        sys.exit("Missing required config line(s) for {} in config file").format(",".join(required))
    if arguments.debug:
        for k in CONFIG:
            print debughead+"{} = {}".format(k,CONFIG[k])
    DESCRIPTOR_APPLICATIONS = {'dssp': CONFIG['DSSP'],
                               'pdb': CONFIG['PDB_PATH'],
                               'model': CONFIG['SWISS_PATH'],
                               'ligand':None,
                               'nucleotide':None,
                               'peptide':None}
                

# Ensure that the required datasets are available
def check_seqs(arguments):
    nomodel = arguments.nomodel
    nopdb = arguments.nopdb
    debug = arguments.debug
    if debug:
        print "DEBUG: I/O: checking for sequence pickles"
    for x in ['SWISS_SEQ',
              'UNP_SEQ',
              'UNP_CANONICAL',
              'UNP_MAP',
              'TRANS_SEQ',
              'SIFTS_SEQ']:
        if (x=="SWISS_SEQ" and not nomodel) or\
           (x!="SWISS_SEQ" and not nopdb):\
            assert path.isfile(CONFIG['SEQ_PATH']+CONFIG[x]),\
            "{}{} not found, run {}pickle_sequences.py".format(
                                                        CONFIG['SEQ_PATH'],
                                                        x,
                                                        CONFIG['SEQ_PATH'])
    if not nopdb:
        assert path.isdir(CONFIG['SIFTS_PATH']), "{} not found. "\
                                   "sifts dataset required for PDB".format(
                                                                    CONFIG['SIFTS_PATH'])
                                                                               
def check_applications(arguments):
    '''
    Check paths for required descriptor files/applications
    name = 'DSSP' or 'PDB' or 'SWISS'
    '''
    names = [DESCRIPTOR_APPLICATIONS[x] for x in arguments.descriptors if x in DESCRIPTOR_APPLICATIONS]
    # Models need to be present, TODO: Add ability to retrieve models from url
    if not arguments.nomodel:
        names.append(DESCRIPTOR_APPLICATIONS['model'])
    for name in names:
        if name is None: continue
        if arguments.debug:
            print "DEBUG: IO: check_applications for {}".format(name)
        assert path.isfile(name) or \
               path.isdir(name),\
            "dependent {} not found".format(name)
                

# These are the names of the all files that will be written to
outfiles = dict()
def define_output(arguments):
    global outfiles
    varfilename = arguments.variants
    basefile = path.splitext(path.basename(varfilename))[0]
    if arguments.debug:
        print "DEBUG: IO: define_output for basename {}".format(basefile)
    outfiles = {
        'skipped': basefile+".skipped",
        'alignments': basefile+".alignments",
        'variants': basefile+".variants",
        'failures': basefile+".failures",
        'completed': basefile+".completed"}                           

##### VARIANT PARSING #####

# Read in the variant file
def load_variants(args):
    debug = args.debug
    debug_head = "DEBUG: IO: load_variants: "
    if debug:
        print debug_head+"loading variants from {}".format(args.variants)
    a = "y"
    varfilename = args.variants
    # If a raw VEP was given, select one consequence per coding variant
    # and then parse
    if args.vep or check_vep(args.variants):
        if args.debug:
            print debug_head+"Converting VEP file to variants"
        if not vep_format_check(args.variants):
            sys.exit("VEP output must be in the txt output format selection. The vep format output is not currently supported")
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
            if args.debug:
                print debug_head+"Previous var file not found, proceeding with process VEP"       
            fullvep = parse_vepfile(args)
            unique_vep = process_vep(fullvep,args.debug)
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
        if debug:
            print debug_head+"Calling parse_varfile with created varfile {}".format(vepoutfile)
        variants = parse_varfile(open(vepoutfile),args)        

    else: # Otherwise, parse it directly
        #Try to open the file
        try:
            infile = open(varfilename)
        except IOError as e:
            sys.exit("Unable to open {}: {}".format(varfilename,e))
        # Try to parse it
        try:
            if debug:
                print debug_head+"Calling parase_varfile with previously created {}".format(varfilename)
            variants = parse_varfile(infile,args)
        except ParseException as e:
            print "Critical failure:"
            sys.exit(e.fullmsg)
    return variants

def check_vep(varfilename):
    '''
    Quick checks the variants file to see if it's a vep file
    in case user forgot to pass -v
    '''
    #Try to open the file
    try:
        infile = open(varfilename)
    except IOError as e:
        sys.exit("Unable to open {}: {}".format(varfilename,e))
    v = infile.read(19)=='#Uploaded_variation'
    infile.close()
    return v

def vep_format_check(varfilename):
    '''
    Quick checks the vep file 
    to make sure it's the txt format selection
    '''
    #Try to open the file
    try:
        infile = open(varfilename)
    except IOError as e:
        sys.exit("Unable to open {}: {}".format(varfilename,e))
    headline = [x.strip() for x in infile.readline().split("\t")]
    v = headline[-1]!='Extra'
    infile.close()
    return v
                
# Process a single variant line
def process_variant(variant,expand,debug):
    '''
    If expand is set:
    Expands usable non-missense variants to list all residues affected.
    Deletions appear as - for all residues affected
    Frameshift appear as X from start of frameshift to end of transcript
    Early stop appear as X from position of stop to end of transcript
    Takes the list of columns of the current variant line and returns
    [position, ref, alt,annotation(for expanding if necessary),varcode]
    '''
    debug_head = "DEBUG: IO: process_variant: "
    if debug:
        print debug_head+"processing line {}".format(variant)
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
                if refaa[0]!=altaa[0] and refaa[1]!=altaa[1]:
                    if debug:
                        print debug_head+"skipping double missense"
                    return ["double missense variant"]
                if debug:
                    print debug_head+"fixing bad missense {} {}/{}".format(position,refaa,altaa)
                position1,position2 = position.split("-")
                #TODO: Adjust these checks if one day desire synonymous mapping
                if refaa[0]==altaa[0]:
                    refaa,altaa = refaa[1],altaa[1]
                    position = int(position2)
                #TODO: Adjust this part if one day desire DNP or TNP
                else:
                    refaa,altaa = refaa[0],altaa[0]
                    position = int(position1)
                if debug:
                    print debug_head+"resolved to {} {}/{}".format(position,refaa,altaa)
            # Anything else is a bad missense and outputs warning and skips
            else:
                if debug:
                    print debug_head+"unusual missense failed: {}/{}".format(refaa,altaa)
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
        if debug:
            print debug_head+"parse error for {}".format(variant)
        raise ParseWarning(
            "variant file","bad variant column: {}".format(variant[3]))               
    except ParseWarning as e:
        print e.fullmsg
        return [e.fullmsg]
    if debug:
        print debug_head+"Formatted variant to {} {}/{} {}".format(position,refaa,altaa,annotation)        
    return [position,refaa,altaa,annotation]

def get_uniprot_from_canonical(trans, canonical, debug):
    '''
    Gets all uniprots assigned to given transcript
    Only provides one isoform, based on priority:
    canonical = YES
    canonical = Unassigned
    canonical = NO
    Then lowest in alphabetical order
    '''
    debug_head = "DEBUG: IO: get_uniprot_from_canonical: "    
    sel = canonical[canonical.Transcript==trans]
    unp = list()
    iso = list()
    if debug:
        print "Found {} uniprot entries for {}".format(sel.shape[0],trans)
    if sel.shape[0] > 0:
        # Filter out secondary if a primary remains
        if sel[~sel.Uniprot.isin(secondary)].shape[0] > 0:
            sel = sel[~sel.Uniprot.isin(secondary)]
        for name,group in sel.groupby('Uniprot'):
            unp.append(name)
            iso.append(list(group.sort_values(['score','Isoform']).Isoform)[0])
    if debug:
        print debug_head+"Settled on {}, {} for {}".format(unp[0],iso[0],trans)
    return unp,iso        
        
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
    col 6 = Uniprot ID, assigned if blank
    col 7 = Uniprot-Isoform, assigned if blank
    col 8 = optional unique variant code
            useful when handling non-missense variants
    Returns a dict with entries:
    KEY = variant identifier (transcript or protein)
    ENTRY = [[uniprot, isoform],[list of variants]]
    '''
    debug_head = "DEBUG: IO: parse_varfile: "
    debug = args.debug
    expand = args.expand
    skipped = list()
    variants = dict()
    count = 0
    canonical_uniprot = None
    secondary = None
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
        # Since splitting by tab, in cases where their editor inserts spaces for tabs,
        # replace any whitespace of 3 or more with tab so it can be parsed.
        whitepattern = re.compile("[ ]{3,}")
        line = whitepattern.sub("\t",line.strip()).split("\t")
        if len(line)==0: continue
        try:
            if len(line)<7:
                if debug:
                    print debug_head+"skipping incomplete line with {} columns".format(len(line))
                    print line
                raise ParseWarning("variant_file","incomplete line: {}".format(
                                               "\t".join(line)))
            consequence,trans,pos,var,protein,unp,iso = line[:7]
            if unp.strip() == "":
                if canonical_uniprot is None:
                    canonical_uniprot = load_canonical(debug)
                    canonical_uniprot = canonical_uniprot[~pd.isnull(canonical_uniprot.Canonical)]
                    scoremap = {'YES': 1, 'Unassigned': 2, 'NO': 3}
                    canonical_uniprot['score'] = map(lambda x: scoremap[x], canonical_uniprot.Canonical)
                    secondary = [x[0] for x in load_sec2prime(debug)]                  
                unp,iso = get_uniprot_from_transcript(trans, canonical_uniprot, secondary, debug) 
            else:
                unp = [unp]
                iso = [iso]
            # For now, just use the first unp returned
            unp = unp[0]
            iso = iso[0]                
            try:
                varcode = line[7]
            except IndexError:
                varcode = "_".join([consequence,
                                   trans,
                                   pos,
                                   var,
                                   protein,
                                   unp])
            if len(unp) == 0:
                if debug:
                    print debug_head+"unable to assign unp entry for {}, skipping {}".format(trans,varcode)
                skipped.append("\t".join(line+["unable to assign uniprot"]))
                continue                    

#            if continue_flag and "{}_{}".format(varcode, in completed: continue                               
            current_var = process_variant(line,expand,debug)           
            #Use ENST identifier if its ensembl otherwise protein identifier
            identifier = line[1] if line[1].startswith("ENST") else line[4]
        except ParseWarning as e:
            print e.fullmsg
            continue      
        if len(current_var)<3:
            if debug:
                print debug_head+"skipping bad variant length {}".format(current_var)
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
def filter_complete(variants,debug):
    debug_head = "DEBUG: IO: filter_complete: "
    global completed
    keptvars = dict()
    if path.isfile(outfiles['completed']) and len(completed)==0:
        if debug:
            print debug_head+"Loading completed variants from".format(outfiles['completed'])
        with open(outfiles['completed']) as infile:
            completed = [x.strip() for x in infile]
    if len(completed)==0:
        if debug:
            print debug_head+"completed file was empty"
        return variants
    if debug:
        print debug_head+"Filtering {} completed variants from {} transcripts".format(len(completed),len(variants))
    kv = 0
    for trans in variants:
        filtered = list()
        for var in variants[trans][-1]:
            varcode = var[-1]
            if varcode in completed:
                continue
            else:
                filtered.append(var)
                kv += 1
        if len(filtered)>0:
            keptvars[trans] = [variants[trans][0],filtered]
    print "skipped {} completed variants".format(len(completed))
    if debug:
        print debug_head+"{} variants in {} transcripts remain".format(kv,len(keptvars))
    return keptvars

# Parse raw VEP output file

def parse_vepfile(args):
    '''
    Reads in a raw VEP file
    First checks if necessary columns are there
    Then, reads it in as a pandas DataFrame, prepares the columns,
    and filters for coding variants
    Returns prepared dataframe
    '''
    # Check to make sure the necessary columns are there
    vepfile = args.variants
    if args.debug:
        print "DEBUG: IO: Parsing VEP file {}".format(vepfile)
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
    if args.debug:
        print "DEBUG: IO: Checking for required VEP columns"
    for x in required_cols:
        if x not in header:
            raise ParseException(
                "Raw VEP file",
                "required column {} missing from header".format(x))
    if args.debug:
        print "DEBUG: IO: Reading CSV into dataframe"   
    df = pd.read_csv(vepfile,sep="\t")[required_cols]
    df.rename(columns={"#Uploaded_variation":"Varcode",
                       "Feature":"Transcript",
                       "ENSP":"Protein"},inplace=True)
    df['Varcode'] = df.Location + ":" + df.Allele
    if args.debug:
        print "DEBUG: IO: Filtering for non-synonymous (AA field >1)"
    df = df[df["Amino_acids"].str.len()>1]
    print "{} unique coding variants loaded".format(df["Varcode"].nunique())
    return df

    
def process_vep(vep,debug):
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
    debug_head = "DEBUG: IO: process_vep: "
    if debug:
        print debug_head+"processing raw df with {} rows".format(len(vep.index))
    # canonical identifies the canonical transcript used by uniprot
    canonical = load_canonical(debug)
    # sec2prime allows filtering out any secondary uniprot ACs
    sec2prime = load_sec2prime(debug)
    
    # Attach the uniprot names
    if debug:
        print debug_head+"merging uniprot names to transcripts"
    vep = vep.merge(canonical,how="left",on=["Transcript","Protein"]) 

    try:
        if len(vep.index) == 0:
            raise ParseException("VEP file","Failed to extract variants from VEP file")
    except ParseException as e:
        sys.exit(e.fullmsg)  

    # Filter out anything that doesn't have a uniprot name
    # Keep those without a uniprot name for use in case they have one
    # assigned by VEP
    if debug:
        print debug_head+"extracting null uniprots into separate df"
    vep_nullunp = vep[vep.Uniprot.isnull()]
    vep = vep[vep.Uniprot.notnull()]
    nrows = len(vep.index)
    if debug:
        print debug_head+"result withunp {} rows; nullunp df {} rows".format(len(vep.index),len(vep_nullunp.index))    

    # Filter any secondary uniprot AC's
    if debug:
        print debug_head+"Removing secondary uniprots"
    secondary = [x[0] for x in sec2prime]
    vep = vep[~vep.Uniprot.isin(secondary)]

    if debug:
        print debug_head+"Result df {} rows".format(len(vep.index))

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
    if debug:
        print debug_head+"Adding position column"
    vep["quickpos"] = vep.Protein_position.astype(str).str.extract('(\d+)',expand=False).astype(int)
    vep_nullunp["quickpos"] = vep_nullunp.Protein_position.astype(str).str.extract('(\d+)',expand=False).astype(int)

    print "selecting unique variants"
    # Loops through steps 1-6:
    # Canonical = YES: ENST then NM then XM
    # Canonical = Unassigned: ENST then NM then XM
    # Canonical = NO: ENST then NM then XM
    
    # Typically, when a transcript is unassigned it means it is the only
    # one and therefore no isoform is designated. However, there is an edge
    # case where it is unassigned even though there are assignments
    # So, need to deal with these after and remove them for now
    if debug:
        print debug_head+"Separating unassigned transcripts"
    vep_unassigned = vep[(vep.Canonical=="Unassigned")\
                          & (pd.to_numeric(vep.nIsoforms)>0)]
    vep = vep[~((vep.Canonical=="Unassigned")\
                          & (pd.to_numeric(vep.nIsoforms)>0))]   
    if debug:
        print debug_head+"Result: {} rows assigned; {} rows unassigned".format(len(vep.index),len(vep_unassigned.index))
    for outer in ["YES","Unassigned","NO"]:
        for inner in ["ENST","NM","XM"]:
            if len(vep.index)==0: break
            if len(vep[(vep.Canonical==outer)\
                    & (vep.Transcript.str.startswith(inner))].index)==0: continue
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
            if debug:
                print debug_head+"Selected {} rows; remaining {} rows".format(len(vep_final.index),len(vep.index))                   
    # Deal with any of potential cases where an unassigned transcript is only hit
    # in a uniprot with multiple isoforms
    vep_unassigned = vep_unassigned[~vep_unassigned.Varcode.isin(vep_final.Varcode)]
    if debug:
        print debug_head+"processing remaining {} unassigned".format(len(vep_unassigned.index))
    for inner in ["ENST","NM","XM"]:
        if len(vep_unassigned.index)==0: break
        group = vep_unassigned[vep_unassigned.Transcript.str.startswith(inner)]\
                                .groupby(["Varcode","Uniprot"])
        current_vars = group.apply(lambda x: x.sort_values(["quickpos","Transcript"],
                                             ascending=False).head(1))
        vep_final = pd.concat([vep_final,current_vars])
        vep_unassigned = vep_unassigned[~vep_unassigned.Varcode.isin(vep_final.Varcode)]
        if debug:
            print debug_head+"Selected {} rows; remaining unassigned {} rows".format(len(vep_final.index),len(vep_unassigned.index))      
    
    # Step 8: select any remaining variants that could not
    # be assigned a uniprot entry based on the read in uniprot data.
    # Only consider instances where VEP provided a uniprot entry.
    # In some cases, the VEP uniprot may not match the current uniprot
    # which is why only those in the sequence datasets provided to this
    # program were considered first and anything left over the uniprot is 
    # taken from VEP
    vep_nullunp = vep_nullunp[~vep_nullunp.Varcode.isin(vep_final.Varcode)]
    if debug:
        print debug_head+"Processing {} rows without uniprot".format(len(vep_nullunp.index))
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
            if debug:
                print debug_head+"Selected {} rows: remaining nounp {} rows".format(len(vep_final.index),len(vep_nullunp.index))
    # Finally, replace all "Unassigned" isoforms with the uniprot entry
    vep_final["Isoform"] = npwhere(vep_final.Isoform=="Unassigned",
                                   vep_final.Uniprot,
                                   vep_final.Isoform)
    print "{} unique variants with uniprot assignment selected".format(
                                            vep_final.Varcode.nunique())
    # report any variants that have multiple assignments
    if vep_final.Varcode.nunique() < vep_final.shape[0]:
        print "The following variant codes have multiple uniprot assignments and will be treated as individual variants for all pairs:"
        dups = vep_final[vep_final.duplicated(subset='Varcode', keep=False)]
        print vep_final[vep_final.Varcode.isin(dups.Varcode)][['Varcode','Uniprot']].to_string(index=False)
    vep = pd.concat([vep,vep_nullunp,vep_unassigned])
    if len(vep.index)>0:
        print "The following {} unique coding variants were unable to be assigned:".format(
            vep.Varcode.nunique())
        print ",".join(set(vep.Varcode))
    vep_final.to_csv("outtmp.tab",sep="\t",header=True,index=False)
    return vep_final.drop_duplicates()    
      
##### DATASET LOADING #####

def open_seqfile(filename):
    filename = CONFIG['SEQ_PATH']+filename
    try:
        infile = open(filename)
    except IOError:
        sys.exit("Failed to open {}".format(filename))
    return infile
            
# Load the transcript sequences
def load_transcripts():
    infile = open_seqfile(CONFIG['TRANS_SEQ'])
    trans = pickle.load(infile)
    print "loaded {} transcript sequences".format(len(trans))
    return trans
    
# Load the sequences required for PDB casting
def load_uniprot():
    infile = open_seqfile(CONFIG['UNP_SEQ'])
    unp = pickle.load(infile)
    print "loaded {} uniprots".format(len(unp))
    return unp

# Load sifts uniprot-PDB alignments
def load_sifts():
    infile = open_seqfile(CONFIG['SIFTS_SEQ'])
    sifts = pickle.load(infile)
    print "loaded {} sifts uniprots".format(len(sifts))
    return sifts

# Load model mappings
def load_models(source):
    if source=="swissmodel":
        infile = open_seqfile(CONFIG['SWISS_SEQ'])
        models = pickle.load(infile)
        print "loaded {} swissmodel uniprots".format(len(models))
    else:
        sys.exit("Critical: unrecognized model-type {}".format(source))
    return models

# Load table of canonical uniprot isoforms
def load_canonical(debug):
    if debug:
        print "DEBUG: IO: load_canonical {}".format(CONFIG['SEQ_PATH']+CONFIG['UNP_CANONICAL'])
    df = pd.read_csv(open(CONFIG['SEQ_PATH']+CONFIG['UNP_CANONICAL']),sep="\t")
    df.drop_duplicates(inplace=True)
    return df

# Load secondary to primary AC mapping
def load_sec2prime(debug):
    if debug:
        print "DEBUG: IO: load_sec2prime {}".format(CONFIG['SEQ_PATH']+CONFIG['UNP_MAP'])
    with open(CONFIG['SEQ_PATH']+CONFIG['UNP_MAP']) as infile:
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
    pdbfile = CONFIG['SIFTS_PATH']+structure+".xml.gz"
    try:
        infile = gzip.open(pdbfile)
    except:
        raise AlignException("parse_sifts", "failed to open {}".format(pdbfile))
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

def remove_previous_files():
    for f in outfiles:
        if path.isfile(outfiles[f]):
            remove(outfiles[f])
            
##### MODELS #####

parser = PDB.PDBParser(PERMISSIVE=1)

def load_model(modelfile,debug):
    '''Loads the model file and generates a fasta
       with dashes for any skipped residue #'s
       Takes a filename and returns a dict with
       filename: model filename
       fasta: fastaseq (of chain, with - for skipped res #)
       chain: chainID (only uses first chain if multiple)
       icodes: sequence of icodes (' ' for none)
       resnums: list of residue numbers'''
    debug_head = "DEBUG: IO: load_model: "
    if debug:
        print debug_head+"Loading model {}".format(modelfile)
    structure = parser.get_structure("Model",modelfile)
    complex_state = model_information(modelfile)
    if debug:
        print debug_head+"model has complex_state: {}".format(complex_state)
    # I'm assuming that all chains in a model are of the target protein
    # I can't find any models that are hetero even if they come from a heteroolig template
    # So, for now I need to stick with this assumption because I can't see any way to determine
    # which chain is the target protein in any cases of a heterooligomer model
    # which is why I think that never happens
    for cc in structure[0]: # For now, only the first chain in the model is taken
        chain = cc
        break
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
    if debug:
        print debug_head+"{} has resnums".format(modelfile, resnums)        
    return {'filename':modelfile,
            'fasta':fasta,
            'chain':chain.get_id(),
            'icodes':icodes,
            'resnums':resnums,
            'complex_state':complex_state}

def model_information(modelfile):
    '''Gets specific model information.
       Regarding complex status and current chain.
       Relevant for non-monomers'''
#    chain = re.search("[.][\w][_][\w]+\.pdb$",modelfile).group(0)[1:2]
    #Get the complex state of the model - assumes model complex state is first
    with open(modelfile) as infile:
        complex_state = 'monomer'
        for line in infile:
            if not line.strip().startswith("REMARK   3  OSTAT"): continue
            complex_state = line.strip().split()[-1]
            break
    return complex_state

def load_artifacts():
    '''
    Loads the artifacts from the artifact file
    Puts into dictionary with keys ligand and peptide
    ligands is a list of 3-character ligand ID
    peptides is a dict of keys pdbid
    each pdbid is a list of chain1_chain2
    '''
    global ARTIFACTS
    x = 0
    ARTIFACTS = {'ligand':[],'peptide':[]}
    if not path.isfile(CONFIG['ARTIFACTS_FILE']):
        raise ParseWarning("Parse artifacts","{} file not found".format(CONFIG['ARTIFACTS_FILE']))
    with open(CONFIG['ARTIFACTS_FILE']) as infile:
        for line in infile:
            line = line.strip()
            if line[0]=="#" or len(line)<3: continue                               
            line = line.split("#")[0]
            if len(line)==3:
                ARTIFACTS['ligand'].append(line)
            elif len(line.split()==3):
                x += 1
                line = line.split()
                if line[0] in ARTIFACTS['peptide']:
                    ARTIFACTS['peptide'][line[0]].append("_".join(line[1:]))
                else:
                    ARTIFACTS['peptide'][line[0]]=["_".join(line[1:])]
            else:
                continue
    print "Loaded artifacts: {} ligands, and {} protein-protein interfaces".format(len(ARTIFACTS['ligand']),x)
                                                                                                
def get_artifacts(pdbid=None):
    if ARTIFACTS is None:
        try:
            load_artifacts()
        except (ParseWarning,IOError) as e:
            print "Warning: Failed to open artifacts file {}, unable to filter artifacts".format(ARTIFACTS_FILE)
            try:
                print "Error message: {}".format(e.fullmsg)
            except:
                print "Error message: {}".format(e)
    return ARTIFACTS
