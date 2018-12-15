from exceptions import *
from IO import CONFIG,get_artifacts
from Bio import PDB
from collections import OrderedDict
import os
import uuid
import pandas as pd
import urllib2
from numpy import where as npwhere
from numpy.linalg import norm as norm
import sys
from numpy import vectorize

INSULIN = "P01308"
BETA2MG = "P61769"
PDB_URL_FORMAT = "https://files.rcsb.org/view/{}.pdb" #uppercase 4char

# Headers are the defined columns for each descriptor
HEADERS = {'dssp': ['SS','SASA_complex','SASA_self'],
           'ligand': ['closest_ligand_distance','ligands_within_5A'],
           'nucleotide': ['closest_nucleotide_distance','closest_nucleotide_type'],
           'peptide': ['closest_chain_distance','chains_within_5A'],
           'unp': ['cleaved','sorting','membrane','binding region',
                   'motif','binding site','ptm','covalent']}
           
# Holders are the placeholder values for each of the descriptor columns
# Typically ? for character descriptors and -999 for numeric
HOLDERS = {'SS':'?','SASA_complex':-999,'SASA_self':-999,
           'closest_ligand_distance':-999,'ligands_within_5A':'?',
           'closest_nucleotide_distance':-999,'closest_nucleotide_type':'?',
           'closest_chain_distance':-999,'chains_within_5A':'?',
           'cleaved':False,'sorting':False,'membrane':False,'binding region':False,
           'motif':False,'binding site':False,'ptm':False,'covalent':False}

# Current set of residues considered DNA or RNA
RNA = ["U","C","A","G"]
DNA = ["DT","DC","DA","DG"]

# Preprocessed dataframe containing all uniprot category features
UNP_DF = None

def load_uniprot_df():
    print "Loading preprocessed uniprot features"
    df = pd.read_pickle(CONFIG['SEQ_PATH']+CONFIG['UNP_FEATURES'])
    print "Loaded {} total residues with uniprot features".format(len(df.index))
    # Previously, I had bad column so if it's there, fix it
    df.rename({'residue_number':'uniprot_position'},axis="columns",inplace=True)
    return df
    
class DictQue(OrderedDict):
    '''
    A dict version of deque that allows
    keyed retrieval and lookup
    '''
    def __init__(self, *args, **kwargs):
        self._max = kwargs.pop('maxlen',0)
        super(DictQue,self).__init__(*args, **kwargs)

    def __setitem__(self, key, value):
        OrderedDict.__setitem__(self, key, value)
        if self._max >0:
            if len(self) > self._max:
                self.popitem(False)
                
class PDBList(object):
    '''
    Collection of the last 100 PDB files parsed.
    To prevent constant reparsing and using up memory,
    the last 200 PDBs parsed are retained in case the next variant
    falls on them.
    The majority of proteins should have less than 200 PDBs.
    In the cases of insulin and beta-2-mg that have up to 1000 PDBs
    the current collections are flushed and they are kept in their own
    collection that can hold up to 1500 structures. Once the standard collection
    reaches 10 items, the insulin and beta2 collections are flushed under
    the assumption that they are no longer needed.
    I'm sure there are other highly represented proteins that may want to add.

    This class is unnecessary right now and isn't used
    It's unlikely it will ever become necessary
    '''
    def __init__(self, artifacts=True):
        self.collection = DictQue(maxlen=200)
        self.dssp = CONFIG['DSSP']
        self.parser = PDB.PDBParser(PERMISSIVE=1,QUIET=True)
        self.insulin = DictQue(maxlen=1500)
        self.beta2mg = DictQue(maxlen=1500)
        self.flushes = 0 # In case want to warn user or something else
        self.artifacts = artifacts
        
    def add_pdb(self, pdbid, target=None):
        '''
        Add a PDB object to the que using parser
        Takes PDBID
        Takes target que
        If at any point getting annotations fails.
        It warns user and adds a holder structure that returns
        -1 for all numeric annotations and ? for all string annotations
        '''
        filename = pdbid
        if target is None:
            target = self.collection
        struct_object = Structure(filename,self.artifacts)
            
        target[filename] = struct_object
        return None
           
    def flush(self):
        '''
        Flushes all que's and updates flush count.
        Prints warning if flushcount>warnthresh
        '''
        self.flushes += 1
        self.collection.clear()
        self.insulin.clear()
        self.beta2mg.clear()
        return None
    
    def add_insulin(self, pdbid):
        '''
        Adds insulin to insulin que if not there.
        If insulin que is empty when this is called,
        flush is called
        '''
        if len(self.insulin)==0:
            self.flush()
        if pdbid not in self.insulin:
            self.add_pdb(pdbid,self.insulin)                                
        return None
    
    def add_beta2mg(self, pdbid):
        '''
        Adds beta2mg to beta2mg que if not there.
        If que is empty when called, it calls flush first
        '''
        if len(self.beta2mg)==0:
            self.flush()
        if pdbid not in self.beta2mg:
            self.add_pdb(pdbid,self.beta2mg)                
        return None
    
    def get_descriptors(self, pdbid, uniprot, d):
        if uniprot==INSULIN:
            if pdbid not in self.insulin:
                self.add_insulin(pdbid)
            return self.insulin[pdbid].descriptors(d)

        elif uniprot==BETA2MG:
            if pdbid not in self.beta2mg:
                self.add_beta2mg(pdbid)
            return self.beta2mg[pdbid].descriptors(d)
        else:          
            if pdbid not in self.collection:
                self.add_pdb(pdbid,self.collection)
            return self.collection[pdbid].descriptors(d)

    def __str__(self):
        return "collection: {}, insulin: {}, beta2mg: {}, flushes {}".format(
            len(self.collection), len(self.insulin), len(self.beta2mg), self.flushes)

def distance(a,b):
    return norm(a-b)
    
class Structure(object):
    '''
    Generate a structure object which uses biopython to load
    the pdb file.
    Descriptors can be attached to it as needed.
    Descriptors are stored in the form of a dataframe with each residue
    having its own row and each descriptor as a column
    Takes a PDB file name during initialization and 
    artifacts = T/F (filter out binding partner artifacts (not implemented)
    Can be initialized without a PDB filename in which case all descriptors
    will be set to filler values.
    Filler descriptors are -999 for all numeric values and ? for all character
    '''
    
    def __init__(self,url=None,filename=None,model=None,artifacts=True):
        '''
        Takes a PDB, if url then uses urlstream as parser input
        if filename, then uses filename as parser input
        if model, then uses parsed model object directly without parsing
        Takes PDB filename and
        artifacts = T/F whether to filter out artifact complexes
        '''
        self.parser = PDB.PDBParser(PERMISSIVE=1,QUIET=True)
        self.filename = filename
        self.dssp = None
        self.ligand = None
        self.nucleotide = None
        self.peptide = None       
        self.artifacts = {'ligand':[],'peptide':[]} if not artifacts else get_artifacts()
        self.dssp_path = CONFIG['DSSP']
        self.id = "Holder"
        self.res_header = ['chain','structure_position','icode']
        self.header = list()
        self.descriptors = None
        self.struct = None
        self.debug = False
        self.debug_head = ""
        if url:
            try:
                self.filename = 'url'
                self.id = url
                pdburl = urllib2.urlopen(PDB_URL_FORMAT.format(url.upper()))
                self.struct = self.parser.get_structure(self.id,pdburl)
                pdburl.close()
                self.header += self.res_header
            except IOError as e:
                print e
                self.id = "Holder"
        elif filename:
            try:
                self.filename = filename
                self.id = os.path.splitext(os.path.basename(filename))[0]            
                self.struct = self.parser.get_structure(self.id,filename)
                self.header += self.res_header
            except IOError as e:
                print e
                self.id = "Holder"
        elif model:
            self.filename = 'model'
            self.struct = model
            self.id = model.get_id()
            self.header += self.res_header

    def flag_debug(self):
        self.debug = True
        self.debug_head = "DEBUG: descriptors: Structure({}): ".format(self.id)

    def add_dssp(self,selected_chains=list()):
        '''
        Adds DSSP features.
        DSSP ignores hetatoms but treats oligomers as single units
        This means that interface residues will have lower than expected SASA
        Therefore, each chain is split into individual chains and DSSP is
        calculated over each chain separately and over the whole oligomer
        DSSP takes files directly so need to create a temporary PDB file for
        each chain
        '''
        if self.debug:
            print self.debug_head+"Adding DSSP"
            print self.debug_head+"Current header: {}".format(self.header)
        if self.id == "Holder":
            c = HOLDERS.keys()
            h = ['?','?','?']+[HOLDERS[x] for x in c]
            c = self.res_header+c
            self.dssp = pd.DataFrame([h],columns=c)[self.res_header+HEADERS['dssp']]
            if self.debug:
                print self.debug_head+"added holder: {}".format(h)
        else:
            genempty = False
            dssp = list()
            try:
                #Calculate DSSP for Oligomer
                #DSSP is pain in that it only takes files. Therefore, if this was created
                # using a model or url, create a temporary file
                if self.debug:
                    print self.debug_head+"Oligomer DSSP"
                if self.filename == "url" or self.filename == "model":
                    ofile = "_".join([uuid.uuid4().hex,self.id])
                    io = PDB.PDBIO()
                    io.set_structure(self.struct)
                    io.save(ofile)
                    if self.debug:
                        print self.debug_head+"wrote file for dssp"
                    try:                        
                        oligomer = dict(PDB.DSSP(self.struct[0],ofile,self.dssp_path))
                    except Exception as e: # DSSP failures generate unnamed exceptions
                        raise DescriptorException("dssp calculation",e)
                    os.remove(ofile)
                else:
                    oligomer = dict(PDB.DSSP(self.struct[0],self.filename,self.dssp_path))
                if self.debug:
                    print self.debug_head+"finished oligomer DSSP"
                class chain_select(PDB.Select): #Needed for extracting each chain
                    def accept_chain(self,c):
                        if c.get_id() == chainid:
                            return True
                        else:
                            return False                               
                #Calculate DSSP for isolated chains
                isolated = dict()
                badchains = list()
                if len(self.struct[0].get_list()) == 1:
                    if self.debug:
                        print self.debug_head+"Single chain, no need to run DSSP on isolated chains"
                    isolated[self.struct[0].get_list()[0].get_id()] = oligomer
                else:
                    if self.debug:
                        print self.debug_head+"Running DSSP on isolated chains"
                    io = PDB.PDBIO()            
                    io.set_structure(self.struct)
                    for chain in self.struct[0]:
                        # Make sure there are residues in this chain otherwise unnamed Exception
                        if len([x for x in chain.get_residues() if x.get_id()[0]==' ' and len(x)>2])==0:
                            badchains.append(chain.get_id())
                            if self.debug:
                                print self.debug_head+"skipping no-res chain {}".format(chain.get_id())
                            continue
                        # Generate random filename
                        chainid = chain.get_id()
                        if chainid not in selected_chains:
                            if self.debug:
                                print self.debug_head+"skipping unwanted chain {}".format(chainid)
                                continue
                        cfile = "_".join([uuid.uuid4().hex,self.id,chainid])
                        # Write isolated chain to random filename
                        # Then calculate DSSP from it
                        io.save(cfile, chain_select())
                        tmp = self.parser.get_structure(chainid,cfile)
                        if self.debug:
                            print self.debug_head+"write file for chain {}, calculating DSSP".format(chainid)
                        try:
                            isolated[chainid] = dict(PDB.DSSP(tmp[0],cfile,self.dssp_path))
                        except Exception: # DSSP failures generate unnamed exceptions
                            print "Warning, DSSP failed for {} isolated chain {}, setting equal to oligomer".format(self.id,chainid)
                            isolated[chainid] = oligomer                        
                        os.remove(cfile)
                #Generate DSSP DataFrame
                if self.debug:
                    print self.debug_head+"generating DSSP dataframe"
                    if len(badchains)>0:
                        print self.debug_head+"skipping badchains {}".format(badchains)
                for res in oligomer:
                    if res[1][0]!=" ": continue
                    c = res[0]
                    if c not in selected_chains or c in badchains: continue
#                    try:
#                        r = int(res[1][1])
#                    except ValueError:
#                        r = -999
                    r = res[1][1]
                    i = res[1][2]
                    ss = oligomer[res][2]
                    try:
                        sasa = round(float(oligomer[res][3]),3)                       
                        isosasa = round(float(isolated[c][res][3]),3)
                    except ValueError: #TODO: Does this every happen?
                        sasa = oligomer[res]
                        isosasa = isolated[c][res][3]                                            
                    dssp.append([c,r,i,ss,sasa,isosasa])                                                                              
            except (OSError,DescriptorException,PDB.PDBExceptions.PDBException) as e: # Generate a holder set on the fly
                print "Warning, DSSP calculation failed for {}: {}".format(self.id,e)
                genempty = True
#                raise ParseWarning("DSSP Calculation","Failed DSSP for {}({})".format(self.id,e))
#            except ParseWarning as e:
                if self.debug:
                    print self.debug_head+"DSSP failed with exception {}".format(e)
            if len(dssp)==0: # Generate a holder set on the fly if it's still empty
                genempty = True
                if self.debug:
                    print self.debug_head+"DSSP set is empty, generating holder set"
            if genempty:
                dssp = list()
                for c in self.struct[0]:
                    for r in c:
                        res = r.get_id()
                        if res[0]!=' ': continue
                        dssp.append([c.get_id(),res[1],res[2],'?',-999,-999])                                 
            self.dssp = pd.DataFrame(dssp,columns=self.res_header+HEADERS['dssp'])
            if self.debug:
                print "Finished dssp datatable has {} rows".format(len(self.dssp.index))
        if self.descriptors is None:
            self.descriptors = self.dssp
        else:
            self.descriptors = self.descriptors.merge(self.dssp,
                                                      on=self.res_header,
                                                      how='outer')
        self.header += HEADERS['dssp']                                                   
        if self.debug:
            print self.debug_head+"Finished DSSP, current header: {}".format(self.header)
    
    def add_ligand(self,selected_residues=list()): #TODO: Need to add skipping NCAA's since MSE is for now an artifact since it's modified MET.
        '''
        Attach the ligand-based descriptors including
        distance to closest ligand and list of ligands within
        4 angstroms
        Can take an option list of selected residues chain_resnum_icode
        to avoid unnecessary computations
        '''
        if self.debug:
            print self.debug_head+"Adding ligand-based descriptor"
            print self.debug_head+"There are {} selected residues".format(len(selected_residues))
        ligdict = {x:list() for x in self.res_header+HEADERS['ligand']}
        if self.id == "Holder":
            c = HOLDERS.keys()
            h = ['?','?','?']+[HOLDERS[x] for x in c]
            c = self.res_header+c
            self.ligand = pd.DataFrame([h],columns=c)[self.res_header+HEADERS['ligand']]
            if self.debug:
                print self.debug_head+"added holder: {}".format(h)
        else:
            if len(selected_residues)==0:
                sel = False
            else:            
                sel = True
            resdict = split_residue_codes(selected_residues)       
            # Get the list of ligands
            ligands = list()
            for c in self.struct[0]:
                for r in c:
                    resid = r.get_id()
                    if resid[0][0]=="H":
                        ligid = resid[0].split("_")[1].strip()
                        if ligid=="HOH": continue #Always skip water
                        if ligid in self.artifacts['ligand']: 
                            if self.debug:
                                print self.debug_head+"Skipping artifact ligand {}".format(ligid)
                            continue
                        ligands.append(r)
            if self.debug:
                print self.debug_head+"Found {} relevant ligands in structure".format(len(ligands))
    
            # If there are no relevant ligands in the structure, generate a holder set
            if len(ligands)==0:
                if self.debug:
                    print self.debug_head+"No ligands in structure, using holders"
                useholder = True
            else:
                useholder = False                            
            for c in self.struct[0]:
                if c.get_id() not in resdict and sel: continue
                for r in c:
                    if r.get_id()[0]!=" ": continue
                    rcode = "_".join(str(i) for i in r.get_id()[1:])
                    if rcode not in resdict[c.get_id()] and sel: continue
                    if useholder:
                        ligdict['chain'].append(c.get_id())
                        ligdict['structure_position'].append(r.get_id()[1])
                        ligdict['icode'].append(r.get_id()[2])
                        for h in HEADERS['ligand']:
                            ligdict[h].append(HOLDERS[h])
                    else:
                        if self.debug:
                            print self.debug_head+"Getting ligand distances for {}_{}".format(c.get_id(),rcode)
                        closest = None
                        within = set()
                        for lig in ligands:
                            ligcode = lig.get_id()
                            ligname = ligcode[0].split("_")[1].strip()
                            # I really hope ligands never have icodes
                            ligcode = "_".join([ligname,c.get_id(),str(ligcode[1])])
                            for a1 in r:
                                for a2 in lig:                           
                                    distance = norm(a1-a2)
                                    if closest is None or distance<closest:
                                        if self.debug:
                                            print self.debug_head+"{} is now closest with dist {}".format(ligcode,distance)
                                        closest = distance
                                    if distance<=5.0:
                                        within.add(ligcode)                                                                                                               
                        ligdict['chain'].append(c.get_id())
                        ligdict['structure_position'].append(r.get_id()[1])
                        ligdict['icode'].append(r.get_id()[2])
                        if closest is None:
                            if self.debug:
                                print self.debug_head+"no closest ligand was found, using holder"
                            closest = HOLDERS['ligand']["closest_ligand_distance"]
                        else:
                            closest = round(closest,2)
                        ligdict['closest_ligand_distance'].append(closest)
                        if len(within)==0:
                            within.add("None")
                        if self.debug:
                            print self.debug_head+"Ligands within cutoff: {}".format(within)
                        ligdict['ligands_within_5A'].append(",".join(within))
                self.ligand = pd.DataFrame.from_dict(ligdict)
                if self.debug:
                    print "Finished ligand datatable has {} rows".format(len(self.ligand.index))
        if self.descriptors is None:
            self.descriptors = self.ligand
        else:
            self.descriptors = self.descriptors.merge(self.ligand,
                                                      on=self.res_header,
                                                      how='outer')
        self.header += HEADERS['ligand']                                                   
        if self.debug:
            print self.debug_head+"Finished ligand, current header: {}".format(self.header)
    
    def add_nucleotide(self, selected_residues = list()):
        '''
        Calculate distance to closest DNA or RNA and indicate whether it's DNA or RNA
        '''
        if self.debug:
            print self.debug_head+"Adding nucleotide-based descriptor"
            print self.debug_head+"There are {} selected residues".format(len(selected_residues))
        nucdict = {x:list() for x in self.res_header+HEADERS['nucleotide']}
        # Return holders if it's a Holder structure
        if self.id == "Holder":
            c = HOLDERS.keys()
            h = ['?','?','?']+[HOLDERS[x] for x in c]
            c = self.res_header+c
            self.nucleotide = pd.DataFrame([h],columns=c)[self.res_header+HEADERS['nucleotide']]
            if self.debug:
                print self.debug_head+"added holder: {}".format(h)
        else:
            if len(selected_residues)==0:
                sel = False
            else:            
                sel = True        
            # First get any DNA/RNA residues in the structure
            nucleotides = list()
            for c in self.struct[0]:
                nucleotides += [x for x in c.get_list() if x.get_resname().strip() in RNA+DNA]
            if len(nucleotides)==0:
                if self.debug:
                    print self.debug_head+"No nucleotide residues found, returning holder"
                useholder = True
            else:
                nucleotide_tups = list()            
                for nuc in nucleotides:
                    t = "DNA" if nuc.get_resname().strip() in DNA else "RNA"
                    nucleotide_tups += [(a,t) for a in nuc if a.element!="H"]
                if self.debug:
                    print self.debug_head+"{} nucleotide residues found".format(len(nucleotides))
                useholder = False
            resdict = split_residue_codes(selected_residues)                                                                                               
            for c in self.struct[0]:           
                if c.get_id() not in resdict and sel: continue                
                for r in c:
                    if r.get_id()[0]!=" ": continue
                    rcode = "_".join(str(i) for i in r.get_id()[1:])
                    if rcode not in resdict[c.get_id()] and sel: continue
                    if useholder:
                        nucdict['chain'].append(c.get_id())
                        nucdict['structure_position'].append(r.get_id()[1])
                        nucdict['icode'].append(r.get_id()[2])
                        for h in HEADERS['nucleotide']:
                            nucdict[h].append(HOLDERS[h])
                    else:
                        if self.debug:
                            print self.debug_head+"Getting nucleotide atom distances for {}_{}".format(c.get_id(),rcode)
                        closest = list()
                        vecd = vectorize(distance)                                                                  
                        atoms = [a for a in r if a.element!="H"]
                        n = nucleotide_tups * len(atoms)
                        atoms = atoms * len(nucleotide_tups)
                        na = [x[0] for x in n]
                        nt = [x[1] for x in n]
                        cv = vecd(atoms,na)
                        mi = cv.argmin()
                        closest = [cv[mi],nt[mi]]
                        nucdict['chain'].append(c.get_id())
                        nucdict['structure_position'].append(r.get_id()[1])
                        nucdict['icode'].append(r.get_id()[2])
                        nucdict['closest_nucleotide_distance'].append(round(closest[0],2))
                        nucdict['closest_nucleotide_type'].append(closest[1])
                self.nucleotide = pd.DataFrame.from_dict(nucdict)
                if self.debug:
                    print "Finished nucleotide datatable has {} rows".format(len(self.nucleotide.index))
        if self.descriptors is None:
            self.descriptors = self.nucleotide
        else:
            self.descriptors = self.descriptors.merge(self.nucleotide,
                                                      on=self.res_header,
                                                      how='outer')
        self.header += HEADERS['nucleotide']                                                   
        if self.debug:
            print self.debug_head+"Finished nucleotide, current header: {}".format(self.header)                                                 

    def add_peptide(self, selected_residues=list()):
        '''
        Calculate distance to closest chain
        '''
        if self.debug:
            print self.debug_head+"Adding peptide-based descriptor"
            print self.debug_head+"There are {} selected residues".format(len(selected_residues))       
        pepdict = {x:list() for x in self.res_header+HEADERS['peptide']}
        # Use holders if it's a Holder structure
        if self.id == "Holder":
            c = HOLDERS.keys()
            h = ['?','?','?']+[HOLDERS[x] for x in c]
            c = self.res_header+c
            self.peptide = pd.DataFrame([h],columns=c)[self.res_header+HEADERS['peptide']]
            if self.debug:
                print self.debug_head+"added holder: {}".format(h)
        else:
            if len(selected_residues)==0:
                sel = False
            else:            
                sel = True        
            # Get all the residue atoms that will be used for neighbor arrays
            residue_atoms = dict()
            for c in self.struct[0]:
                ca = [a for r in c for a in r if r.get_id()[0]==' ' and a.element!='H']
                if len(ca)>0:
                    residue_atoms[c.get_id()] = ca
                    
            # Check if it's only 1 chain and if so use holders
            # Need to check after gathering atoms to avoid ligand chains, etc.
            if len(residue_atoms.keys())<2:
                if self.debug:
                    print self.debug_head+"Only one chain with residues found, returning holder"
                useholder = True
            else:
                if self.debug:
                    print self.debug_head+"{} neighbor chains found".format(len(residue_atoms.keys()))
                useholder = False
            resdict = split_residue_codes(selected_residues)                                                                                               
            for c in self.struct[0]:           
                if c.get_id() not in resdict and sel: continue                
                for r in c:
                    if r.get_id()[0]!=" ": continue
                    rcode = "_".join(str(i) for i in r.get_id()[1:])
                    if rcode not in resdict[c.get_id()] and sel: continue
                    if useholder:
                        pepdict['chain'].append(c.get_id())
                        pepdict['structure_position'].append(r.get_id()[1])
                        pepdict['icode'].append(r.get_id()[2])
                        for h in HEADERS['peptide']:
                            pepdict[h].append(HOLDERS[h])
                    else:
                        if self.debug:
                            print self.debug_head+"Getting peptide distances for {}_{}".format(c.get_id(),rcode)
                        vecd = vectorize(distance)                                                                  
                        atoms = [a for a in r if a.element!="H"]
                        # Generate a list that contains all atoms for each chain
                        # duplicated for all atoms in current res
                        # could probably just use broadcasting
                        current_neighbors = list()
                        for nc in residue_atoms:
                            if nc == c.get_id(): continue
                            current_neighbors += residue_atoms[nc]
                        if self.debug:
                            print self.debug_head+"Getting distances for {} atoms".format(len(current_neighbors))
                        cna = current_neighbors * len(atoms)
                        curatoms = atoms * len(current_neighbors)
                        cv = vecd(curatoms,cna)
                        mindist = cv.min()
                        if self.debug:
                            print self.debug_head+"Min dist = {}".format(mindist)                        
                        if mindist < 5:
                            # Put the chains in a temporary dataframe to collect those with min distance less than 5
                            tmpp = pd.DataFrame.from_dict({'c':[cl.get_parent().get_parent().get_id() for cl in cna],
                                                           'd':cv})
                            if self.debug:
                                tmpp['id'] = cna
                            cutoffchains = set(tmpp[tmpp.d<5].c)
                            if self.debug:
                                cutoffatoms = set(tmpp[tmpp.d<5].id)
                                print self.debug_head+"Residues within 5A: {}".format([dx.get_parent().get_full_id() for dx in cutoffatoms])
                        else:
                            cutoffchains = set(['None'])                            
                        pepdict['chain'].append(c.get_id())
                        pepdict['structure_position'].append(r.get_id()[1])
                        pepdict['icode'].append(r.get_id()[2])
                        pepdict['closest_chain_distance'].append(round(mindist,2))
                        pepdict['chains_within_5A'].append(",".join(cutoffchains))
                self.peptide = pd.DataFrame.from_dict(pepdict)
                if self.debug:
                    print self.debug_head+"Finished peptide datatable has {} rows".format(len(self.peptide.index))
        if self.descriptors is None:
            self.descriptors = self.peptide
        else:
            self.descriptors = self.descriptors.merge(self.peptide,
                                                      on=self.res_header,
                                                      how='outer')
        self.header += HEADERS['peptide']                                                   
        if self.debug:
            print self.debug_head+"Finished peptide, current header: {}".format(self.header)                       
                     
    def attach_descriptors(self, partner):
        '''
        Attach descriptors to partner df
        Requires a partner df
        Does this by merge (chain/res#/icode) if not a holder
        If holder, attaches using bind
        returns the merged/binded result
        '''
        if self.debug:
            print self.debug_head+"Attaching descriptors to df with {} rows".format(len(partner.index))
        if self.descriptors is None:
            return partner
        if self.id == "Holder":
            # Generate a holder descriptor table
            desc_dict = dict()
            nrow = len(partner.index)
            for d in self.header:
                desc_dict[d] = [HOLDERS[d]]*nrow
            holder_desc = pd.DataFrame.from_dict(desc_dict)
            holder_desc.set_index(partner.index)
            if self.debug:
                print self.debug_head+"Concatenating holder descriptors"
            return pd.concat([partner,holder_desc[self.header]],axis=1)          
        else:
#            if self.filename=="url":
#                self.descriptors[self.header].to_csv("test1")
#                partner.to_csv("test2")
#                print self.res_header
#                print self.descriptors[self.res_header]
#                print partner[['structure','icode','chain']]
#                with open("tmp","a") as outfile:
#                    outfile.write("{}\ndf:{} = {}\npartner:{} = {}\n".format(self.id,list(self.descriptors[self.header]),self.descriptors[self.header].dtypes,list#(partner),partner.dtypes))
#            print partner.merge(self.descriptors[self.header],how='left',on=self.res_header)
#            print self.descriptors[self.header]
#            print partner
            newdf = partner.merge(self.descriptors[self.header],
                                 how='left',
                                 on=self.res_header)
            # Fill in any residues that were missing descriptors
            for x in self.header:
                if x in self.res_header: continue
                newdf[x] = npwhere(newdf[x].isnull(),HOLDERS[x],newdf[x])
            return newdf

    def get_dssp(self,chain,resnum,icode=" "):
        '''
        Get the DSSP for a particular residue
        Takes chain, resnum, and icode
        Returns a tuple of dicts with keys
        ss (sec structure) and sasa (acc surface area)
        First dict is isolated chain
        Second dict is from oligomer
        '''
        if self.id == "Holder":
            print "Warning, getting dssp from holder, returning holder value"
            return pd.DataFrame([[chain,resnum,icode,"?",-999,-999]],
                                columns=self.res_header+HEADERS['dssp'])
        elif self.dssp is None:
            print "Warning, dssp not run returning holder value"
            return pd.DataFrame([[chain,resnum,icode,"?",-999,-999]],
                                columns=self.res_header+HEADERS['dssp'])
        else:
            return self.dssp[(self.dssp.chain==chain) &
                             (self.dssp.structure_position==resnum) &
                             (self.dssp.icode==icode)]                                

    def __str__(self):       
        return "{} ({}), dssp: {}, ligand: {}, nucleotide: {}, peptide: {}".format(
            self.id, self.filename, self.dssp is not None, self.ligand is not None,
            self.nucleotide is not None, self.peptide is not None)

def add_structure_descriptors(df,descriptors,debug):
    '''
    Add artifacts to dataframe
    Takes dataframe and a list of descriptors
    '''
    debug_head = "DEBUG: descriptors: add_structure_descriptors: "
    if debug:
        print debug_head+"adding descriptors to df with {} rows".format(len(df.index))
    try:
       assert 'structure_position' in list(df.columns.values)
       assert 'icode' in list(df.columns.values)
       assert 'structure' in list(df.columns.values)
       assert 'chain' in list(df.columns.values)
    except AssertionError:
       raise DescriptorException("initialization",
                                 "dataframe passed to add_structure_descriptors missing required column(s)")
    if 'artifacts' in descriptors:
        a = False
    else:
        a = True
    structures = df.structure.drop_duplicates()
    if debug:
        print debug_head+"After dropping duplicates, df has {} rows".format(len(df.index))
    allstructs = None
    for model in structures:
        # Get chains and residues of interest to minimize calculation time
        tmp_df = df[df.structure==model][['chain','structure_position','icode']]
        chains_of_interest = list(tmp_df.chain.drop_duplicates())
        if debug:
            print debug_head+"Restricting to chains of interest: {}".format(chains_of_interest)

        residues_of_interest = list(tmp_df.chain.astype(str)+"_"+tmp_df.structure_position.astype(str)+"_"+tmp_df.icode.astype(str))
        if debug:
            print debug_head+"Restricting to residues of interest: {}".format(residues_of_interest)                                              
        if model=="Uniprot":
            # Uniprot is holder structure
            if debug:
                print debug_head+"adding descriptors to uniprot"
            structure = Structure()
        elif len(model)==4:
            # PDB needs to be retrieved
            if debug:
                print debug_head+"adding descriptors to {}".format(model)
            structure = Structure(url=model,artifacts=a)           
        else:
            # Model is local file
            if debug:
                print debug_head+"adding descriptors to {}".format(model)
            pdbfile = model
            structure = Structure(filename=pdbfile,artifacts=a)                                   
        if debug:
            structure.flag_debug()
        if 'dssp' in descriptors:
            structure.add_dssp(chains_of_interest)
        if 'ligand' in descriptors:
            structure.add_ligand(residues_of_interest)
        if 'peptide' in descriptors:
            structure.add_peptide(residues_of_interest)
        if 'nucleotide' in descriptors:
            structure.add_nucleotide(residues_of_interest)
        current_df = structure.attach_descriptors(df[df.structure==model][
                                                     ['structure',
                                                     'chain',
                                                     'structure_position',
                                                     'icode']])
        if allstructs is None:
            allstructs = current_df
        else:
            allstructs = pd.concat([current_df,allstructs])
        if debug:
            print "Current allstructs df has {} rows".format(len(allstructs.index))                                                                             
    return pd.merge(df,allstructs,how='left',
                    on=['structure','chain','structure_position','icode'])

def split_residue_codes(rescodes):
    resdict = dict()
    for r in rescodes:
        c,rn,i = r.split("_")
        ri = "{}_{}".format(rn,i)
        if c in resdict:
            resdict[c].append(ri)
        else:
            resdict[c] = [ri]
    return resdict

## Non-structure based descriptors
def add_uniprot_descriptors(df, debug):
    '''
    Adds uniprot category flags from preprocessed uniprot feature set
    Can handle holder since it just relies on residue having a uniprot position
    '''
    current_unp = df.iloc[0].uniprot
    global UNP_DF
    if UNP_DF is None:
        UNP_DF = load_uniprot_df() #TODO: move this to do at the beginning in main
            
    debug_head = "DEBUG: descriptors: add_uniprot_descriptors: "
    if debug:
        print debug_head+"Adding uniprot category descriptors"
    unpdf = UNP_DF[UNP_DF.uniprot==current_unp]
    if unpdf.empty: # Probably don't need to do this as nulls are filled in at the end
        #No features assigned to this uniprot, need to use holders
        if debug:
            print debug_head+"No uniprot features found for uniprot {}, using holders".format(current_unp)
        holderdict = {'uniprot_position': list(df.uniprot_position.unique())}
        totalres = len(holderdict['uniprot_position'])
        for h in HEADERS['unp']:
            holderdict[h] = [HOLDERS[h]]*totalres
        unpdf = pd.DataFrame.from_dict(holderdict)
    elif debug:
        print debug_head+"Uniprot features found for {} residues in {}".format(len(unpdf.index),current_unp)
    ndf = pd.merge(df,unpdf,how="left",on='uniprot_position')
    for h in HEADERS['unp']:
        ndf[h] = npwhere(ndf[h].isnull(),HOLDERS[h],ndf[h])
    return ndf        
