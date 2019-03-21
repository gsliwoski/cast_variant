# Tutorial
Please see tutorial.pdf for a walkthrough and description of the datasets used and an example run of the cast_variant. Example input and output files are also precreated in main directory as example_VEP.

# cast_variant
Cast variant onto protein structure (PDB or model)

Aligns with SWISSMODEL models and PDBs via uniprot-sifts alignments

Returns complete alignments and variant conversions

Right now, skips any NMD, start changes, stop loss, and multi-AA missense changes

## Usage

### Required input

Input style 1: List of pre-formatted variants with all necessary columns:

Takes list of variants that have already been processed with effect predictor (hopefully VEP38)

Tab delimited with one variant per line (no header)

Anything after # is ignored

Incomplete lines will be reported with warning and put into .skipped

Format:

Consequence ENST#   Position    Amino\_acids(R/A)    ENSP#   Uniprot Uniprot-isoform Variant\_Code

**Example:**

>\#Ignore this line
>
>this line incomplete #ignore this part makes it long enough tho
>
>missense_variant        ENST00000377038 235     V/I     ENSP00000366237 O00273  O00273-1        9

>frameshift_variant      ENST00000377038 187-188 QS/X    ENSP00000366237 O00273  O00273-1        10

>missense_variant,splice\_region\_variant  ENST00000358465 989     S/L     ENSP00000351250 Q9UPN9  Q9UPN9-1        23

>missense_variant        ENST00000358465 974     Y/H     ENSP00000351250 Q9UPN9  Q9UPN9-1        24

>missense_variant        ENST00000358465 961     V/M     ENSP00000351250 Q9UPN9  Q9UPN9-1        25

>missense_variant        ENST00000358465 938     E/Q     ENSP00000351250 Q9UPN9  Q9UPN9-1        26

>missense_variant        ENST00000358465 898     D/V     ENSP00000351250 Q9UPN9  Q9UPN9-1        27

>stop_gained,NMD     ENST00000358465 865     R/*     ENSP00000351250 Q9UPN9  Q9UPN9-1        28

>inframe_deletion        ENST00000358465 720     S/-     ENSP00000351250 Q9UPN9  Q9UPN9-1        36

**With this input the following variants will be skipped without expand variant flag (untested)**

* frameshift
* inframe_deletion
* stop_gained = skipped regardless since NMD

Input style 2: Raw VEP output

https://useast.ensembl.org/Tools/VEP

When giving variants from VEP, it will detect this by checking for #Uploaded_variation as the first cell, so if it fails to detect this, make sure you have a header in the VEP file.

When running VEP, the following settings beyond default are recommended:

Transcript database to use:

Ensembl/GENCODE and RefSeq transcripts

Identifiers: check all of them (REQUIRED identifiers: Protein, UniProt)

Variants and frequency data: Leave as default

Additional annotation: Leave as default

Predictions: Leave as default

Filtering options: Return results for variants in coding regions only (Leave restrict results as show all results)

Example VEP command line that uses all these settings:

./vep --af --appris --biotype --ccds --check_existing --coding_only --distance 5000 --hgvs --merged --polyphen b --protein --pubmed --regulatory --sift b --species homo_sapiens --symbol --tsl --uniprot --cache --input_file [input_data]

Must use flag -v or --vep to indicate it's vep input.

First, a single representation is selected for each unique variant based on priorities below.

Selected variants are written to input_filename.vep_selected for use directly on any following runs.

Then, selected variants are processed as a normal set of variants.

**Priority selection**

Selection relies on uniprot canonical transcript identification and therefore uniprot mapping files must be premade.

Uniprots are first assigned based on transcript library and, in cases where no uniprot can be assigned, vep assignments are used

Anything without an assigned uniprot is filtered.

Canonical sequences are always defined as the sequence presented on the webpage for the uniprot entry.

Each unique variant is represented by the first occurence in the following list of priorities:

(Any cases that return multiple identical effects on different transcripts, select highest transcript #)

1. Canonical ENST as defined by uniprot canonical library
2. Canonical NM as defined by uniprot canonical library
3. ENST with unassigned canonical (in many cases where there is only one isoform, uniprot will not assign an isoform number to canonical sequence and so these represent the only sequence for this uniprot)
4. NM with unassigned canonical
5. XM transcript that has an assigned uniprot

**Result:**

* All missense will be cast
* frameshift and inframe_deletion will only be cast if expand is set
* stop\_gained,NMD will be skipped either because of stop\_gained or NMD
* Incomplete line will be skipped with warning

### Required datasets

**See /datasets/downloading_datasets.txt for complete downloading instructions**

1. A set of swissmodel structures in the directory hierarchy as created through direct retrieval of SWISSMODEL database
2. A set of all SIFTS alignments downloaded from SIFTS database. (split_xml)
3. A set of all PDB structures downloaded from rcsb (.pdb format)
4. Homo_sapiens.GRCh38.pep.all.fa.gz = FASTA of all human protein sequences from ENSEMBL (in datasets folder
5. pdb_chain_uniprot.tsv.gz = Summary file from SIFTS that lists all alignments in SIFTS (in datasets folder)
6. uniprot_sprot.fasta.gz = FASTA of all canonical human protein sequences from Uniprot (in datasets folder)
7. refseq_mammalian.fasta.gz = fasta of all refseq sequences (in datasets folder)
8. uniprot_sprot_human.dat (in datasets folder)
9. uniprot_sprot.xml (in datasets folder)
10. sifts.pickle = generated in processing step (datasets folder)
8. swissmodel.pickle = generated in processing step (datasets folder)
9. trans.pickle = generated in processing step (datasets folder)
10. uniprot_canonical_isoforms.tab = generated in processing step, used to filter secondary uniprot IDs (datasets folder)
11. uniprot_sec2prim_ac.txt = generated in processing step, used to match transcripts to uniprot IDs (datasets folder)
12. unp_features.pickle = generated in processing step (datasets folder)
13. unp.pickle = generated in processing step (datasets folder)

**config.sys** = contains all paths for necessary sequences and applications. Each setting is space separated ID and value.


### Command line

```
./cast_variant.py variant_list.tab
```

**Optional Flags:**

**--expand, -e** = expands non-missense variants to encompass all residues affected

Note: in-frame deletions have alternate amino acid '-' for deleted residues; frameshift and stop_gain have X as alternate AA from the start position to the end of the transcript.

**--nopdb, -p** = skip casting to PDB structures

**--nomodel, -m** = skip casting to models

**--nouniprot,-u** = only relevant for casting to models since PDB casting uses SIFTS which relies on uniprot sequence. This is so you still get the uniprot positions if you cast to a  model which is done directly.

**--list,-l FILENAME** = include a set of local models not part of SWISSMODEL. Filename contains list of all models (one per line) with transcript[tab]full_path_file.pdb. This accepts any number of unique transcripts. Note, this hasn't been thoroughly tested but appears to work.

**--completed,-c** = suppress continue function. By default, a .completed file is created which tracks which variants have already been cast so, in case of interruption, can continue from where it left off.

**--num_procs \#, -n \#** = number of multiprocesses to run. Default = 1

**--noalign, -a** = suppress creation of full alignment file. By default, writes all full alignments to separate .alignments file. Set if only care about variants.

**--vep, -v** = indicates input file is raw VEP output and should be pre-processed into usable variants.

**--descriptors, -d []** = attach desired descriptors to each variant, see Descriptors section below

**--debug, -x** = Very verbose output for debugging

### Descriptors

Descriptors are provided as a comma-separated list.

Any unrecognized descriptors are ignored

Using 'all' will attach all possible descriptors

By default, for complex-based descriptors including protein-protein interfaces and protein-ligand interfaces, known or suspected artifacts are filtered. These are defined in the ARTIFACTS_FILE and can be changed as desired.

Alternatively, including 'artifacts' in the descriptor list will not filter for any artifacts

**Available descriptors**

Note: All structure-based descriptors are drawn from the first model in a structure

Note: descriptor columns are attached to the .variant output file only (not the .alignment)

**dssp** = Runs DSSP on the whole structure, and then on the variant's isolated chain. returns SS (secondary structure), SASA_complex (relative solvent accessibility in entire structure), SASA_self (relative solvent accessibility for the isolated chain)

**nucleotide** = returns closest_nucleotide_distance (minimum atom-pair distance from variant residue to any RNA or DNA residues in the structure), closest_nucleotide_type (either DNA or RNA or None)

**ligand** = May be subject to artifact filtering. Returns closest_ligand_distance (minimum atom-pair distance from variant residue distance to any ligands in structure), ligands_within_5A (comma separated list of all ligands as ID_chain_resnum of all ligands in structure with atom-pair distance less than 5 with variant residue).

**peptide** = May be subject to artifact filtering. Returns closest_chain_distance (minimum atom-pair distance from variant residue to any other chains in structure), chains_within_5A (comma separated list of chains with minimum atom-pair distance from variant residue less than 5)

**unp** = Uniprot-based annotations. Downloading uniprot XML and running ./prepare_uniprot_features.py is required to use this. True/False features are added for pre-specified uniprot groups. These groups include

1. cleaved = initiator methionine, propeptide [regions of the protein that are not present in its mature state]
2. sorting = signal, transit peptide, signal peptide [regions that are used to deliver protein to proper organelle/place in cell
3. membrane = transmembrane, intramembrane, transmembrane region, intramembrane region [regions found within a membrane]
4. binding region = calcium binding, DNA binding, DNA-binding region, nucleotide binding, calcium-binding region, phosphate-binding region, nucleotide phosphate-binding region, lipid moiety-binding region [regions that bind various ligands. distinguished from sites since they are typically much larger and less precise]
5. motif = motif, short sequence motif [predefined motifs less than 20 residues in length]
6. binding site = active site, metal binding, binding site, site, metal ion-binding site [precise binding sites for various ligands, typically smaller than binding regions]
7. ptm = modified residue, lipidation, glycosylation, glycosylation site [residues that undergo post-transcriptional modification either to anchor the protein to a specific site (ie: lipidation/glycosylation) or regulation (phosphorylation)]
8. covalent = disulfide bond, cross-link [residues involved in covalent interactions with other residues typically important for conformational stability] 

The following features are known and not used for any category either because they are too imprecise or are not currently considered useful information:

chain, peptide, domain, repeat, coiled coil, coiled-coil region, compositional bias, non-standard residue, alternative sequence, natural variant, mutagenesis, mutagenesis site, sequence uncertainty, sequence conflict, non-adjacent residues, non-terminal residues, helix, turn, beta strand, strand, sequence variant, splice variant, compositionally biased region, zinc finger region, non-standard amino acid, non-terminal residue, topological domain, region, region of interest, non-consecutive residues

Note: for any features that are defined by start and end residues, all residues in between are flagged. Disulfide bonds are flagged for both cysteines involved.

### Custom local models

If user wants to map to any local models not part of SWISSMODEL, must list all models (one per line) in a TAB delimited file in format

Transcript1  Full_model_path1.pdb
Transcript1  Full_model_path2.pdb
Transcript2  Full_model_path3.pdb

Note this hasn't been thoroughly tested but has worked so far.

### Results

Assuming noalign and completed is not set, will create 2 files and potentially 3 more depending on variant list.

Inputfilename.variants = all variant positions aligned with uniprot, PDBs, and models along with any selected descriptors

Inputfilename.alignments = full alignments for transcript and all partners

**Potential Files**

Inputfilename.skipped = any variants that were not cast and, if possible, a short reason why

Inputfilename.failures = any variants that failed during casting with a reason why when possible otherwise a traceback for general exceptions

Inputfilename.completed = list of completed transcripts

