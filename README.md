# cast_variant
Cast variant onto protein structure (PDB or model)

Version 1:

Takes list of variants that have already been processed with effect predictor (hopefully VEP38)

Aligns with SWISSMODEL models and PDBs via uniprot-sifts alignments

Returns complete alignments and variant conversions

Right now, skips any NMD, start changes, stop loss, and multi-AA missense changes

## Usage

### Required input

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

**Result:**

* All missense will be cast
* frameshift and inframe_deletion will only be cast if expand is set
* stop\_gained,NMD will be skipped either because of stop\_gained or NMD
* Incomplete line will be skipped with warning

### Required datasets

1. A set of swissmodel structures in the directory hierarchy as created through direct retrieval of SWISSMODEL database
2. A set of all SIFTS alignments in hierarchy as directly downloaded from SIFTS database. (split_xml)
3. Homo_sapiens.GRCh38.pep.all.fa.gz = FASTA of all human protein sequences from ENSEMBL (put in /sequences/)
4. pdb_chain_uniprot.tsv.gz = Summary file from SIFTS that lists all alignments in SIFTS (put in /sequences/)
5. uniprot_sprot.fasta.gz = FASTA of all canonical human protein sequences from Uniprot (put in /sequences/)
6. refseq_mammalian.fasta.gz = fasta of all refseq sequences

Note: See /sequences/downloading_datasets.txt for more information

**Run:**
```
/sequences/pickle_sequences.py
```

to prepare files.

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

**--completed,-c** = suppress continue function. By default, a .completed file is created which tracks which variants have already been cast so, in case of interruption, can continue from where it left off.

**--num_procs \#, -n \#** = number of multiprocesses to run. Default = 1

**--noalign, -a** = suppress creation of full alignment file. By default, writes all full alignments to separate .alignments file. Set if only care about variants.

### Results

Assuming noalign and completed is not set, will create 2 files and potentially 3 more depending on variant list.

Inputfilename.variants = all variant positions aligned with uniprot, PDBs, and models

Inputfilename.alignments = full alignments for transcript and all partners

**Potential Files**

Inputfilename.skipped = any variants that were not cast and, if possible, a short reason why

Inputfilename.failures = any variants that failed during casting with a reason why when possible otherwise a traceback for general exceptions

Inputfilename.completed = list of completed transcripts


