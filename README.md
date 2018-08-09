# cast_variant
Cast variant onto protein structure (PDB or model)

Version 1:

Takes list of variants that have already been processed with effect predictor (hopefully VEP38)

Aligns with SWISSMODEL models and PDBs via uniprot-sifts alignments

Returns complete alignments and variant conversions

Right now, skips any NMD, start changes, stop loss, and multi-AA missense changes
