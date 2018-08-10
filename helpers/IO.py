from os import path
from exceptions import *
import sys
import pickle

##### INITIALIZATION FUNCTIONS #####

# Define what the sequence files are
SEQ_PATH = "./sequences/" # Where pickles are
SWISS_SEQ = "swissmodel.pickle" # Model mapping
UNP_SEQ = "unp.pickle" # Uniprot sequences
TRANS_SEQ = "trans.pickle" # Transcript sequences
SIFTS_SEQ = "sifts.pickle" # uniprot - pdb mapping
SIFTS_PATH = "/dors/capra_lab/data_clean/sifts/xml/" # Where sifts alignments are

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
        'conversions': basefile+".conversions",
        'failures': basefile+".failures",
        'completed': basefile+".completed"}

# These are the header and column formats for output
# For models, the unp columns will be NA and therefore
# filtering out models will require uniprot_position=="NA" filter
CONV_HEADER = ["transcript","uniprot","isoform",
               "transcript_identity", "transcript_position",
               "uniprot_position", "structure_position", "icode",
               "transcript_aa", "uniprot_aa", "structure_aa",
               "ref_aa", "alt_aa", "structure_identity",
               "structure", "chain"]

ALN_HEADER = ["transcript","uniprot","isoform",
              "transcript_identity", "transcript_position",
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
                pass              
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

##### PDB LOADING #####


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
