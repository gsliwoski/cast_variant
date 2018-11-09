from exceptions import *
from IO import DSSP
from Bio import PDB
from collections import OrderedDict
import os
import uuid

INSULIN = "P01308"
BETA2MG = "P61769"

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
    '''
    def __init__(self, artifacts=True):
        self.collection = DictQue(maxlen=200)
        self.dssp = DSSP
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
        pdbid = os.path.splitext(os.path.basename(filename))[0]
        try:
            struct = self.parser.get_structure(pdbid,filename)
            struct_object = Structure(struct,self.artifacts,filename)
            struct_object.add_dssp()
        except ParseWarning as e:
            print e.fullmsg
            struct_object = Structure(None,self.artifacts)
            
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
    
    def get_descriptors(self, pdbid, uniprot):
        if uniprot==INSULIN:
            if pdbid not in self.insulin:
                self.add_insulin(pdbid)
            return self.insulin[pdbid]

        elif uniprot==BETA2MG:
            if pdbid not in self.beta2mg:
                self.add_beta2mg(pdbid)
            return self.beta2mg[pdbid]
        else:          
            if pdbid not in self.collection:
                self.add_pdb(pdbid,self.collection)
            return self.collection[pdbid]

    def __str__(self):
        return "collection: {}, insulin: {}, beta2mg: {}, flushes {}".format(
            len(self.collection), len(self.insulin), len(self.beta2mg), self.flushes)

class Structure(object):
    '''
    Container for structure-based descriptors.
    Takes a parsed PDB object or nothing.
    If provided nothing, generates a holder structure
    that returns -1 for all numeric descriptors
    and ? for all string descriptors
    '''
    
    def __init__(self, struct=None,artifacts=True,filename=None):
        '''
        Takes PDB Structure object and
        artifact flag (filter artifact complexes)
        '''
        self.struct = struct
        self.filename = filename
        self.dssp = None
        self.ligand = None
        self.nucleotide = None
        self.peptide = None       
        self.artifacts = artifacts
        self.dssp_path = DSSP
        print struct.get_id()
    def add_dssp(self):
        '''
        Adds DSSP features.
        DSSP ignores hetatoms but treats oligomers as single units
        This means that interface residues will have lower than expected SASA
        Therefore, each chain is split into individual chains and DSSP is
        calculated over each chain separately and over the whole oligomer
        DSSP takes files directly so need to create a temporary PDB file for
        each chain
        '''
        if self.struct is None:
            print "Can't add DSSP to holder Structure"
            return None
        try:
            oligomer = PDB.DSSP(self.struct[0],self.filename,self.dssp_path)
        except OSError as e:
            raise ParseWarning("DSSP Calculation","Failed DSSP for {}({})".format(self.filename,e))
#        chainid = "A" #Used for isolating chains
        class chain_select(PDB.Select): #Needed for extracting each chain
            def accept_chain(self,c):
                if c.get_id() == chainid:
                    return True
                else:
                    return False                    
        #TODO: add exception handling to this part
        cparser = PDB.PDBParser(PERMISSIVE=1,QUIET=True)
        self.dssp = {'ALL':dict(oligomer)}
        if len(self.struct[0].get_list()) == 1:
            self.dssp[self.struct[0].get_list()[0].get_id()] = self.dssp['ALL']
        else:
            for chain in self.struct[0]:
                io = PDB.PDBIO()
                io.set_structure(self.struct)
                # Generate random filename
                cfile = uuid.uuid4().hex
                chainid = chain.get_id()
                io.save(cfile, chain_select())
                tmp = cparser.get_structure(chain,cfile)
                self.dssp[chainid] = dict(PDB.DSSP(tmp[0],cfile,self.dssp_path))
                os.remove(cfile)

    def get_holder(self,source=None):
        '''
        Placeholder features whenever a retrieval fails
        1 argument, source: Sources: DSSP
        If no source given, returns None
        '''
        if source=="DSSP":
            return (-1,'?','?',-1,-1,-1,-1,-1,-1,-1,-1,-1,-1)
        else:
            return None            
        
    def get_dssp(self,chain,resnum,icode=" "):
        '''
        Get the DSSP for a particular residue
        Takes chain, resnum, and icode
        Returns a tuple of dicts with keys
        ss (sec structure) and sasa (acc surface area)
        First dict is isolated chain
        Second dict is from oligomer
        '''
        res_tup = (chain,(" ", resnum, icode))
        try:
            if self.struct is None:
                raise KeyError("no loaded pdb")
            dssp_tup_chain = self.dssp[chain][res_tup]
            dssp_tup_olig = self.dssp['ALL'][res_tup]
        except KeyError as e:
            print "Warning, unable to calculate DSSP for {}:{}{}({})".format(chain,resnum,icode,e)
            print "Returning Placeholder values"
            dssp_tup_chain = self.get_holder("DSSP")
            dssp_tup_olig = self.get_holder("DSSP")
        return ({'ss': dssp_tup_chain[2], 'sasa': dssp_tup_chain[3]},
                {'ss': dssp_tup_olig[2], 'sasa': dssp_tup_olig[3]})

    def __str__(self):
        if self.struct is None:
            sid = "None"
        else:
            sid = self.struct.get_id()
        return "{} ({}), dssp: {}, ligand: {}, nucleotide: {}, peptide: {}".format(
            sid, self.filename, self.dssp is not None, self.ligand is not None,
            self.nucleotide is not None, self.peptide is not None)
