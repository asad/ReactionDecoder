****************************
Reaction Decoder Tool (RDT)
****************************

1) Atom Atom Mapping (AAM) Tool
2) Reaction Annotator (Extract Bond Changes, Identify & Mark Reaction Centres) and
3) Reaction Comparator (Reaction Similarity based on the Bond Changes, Reaction Centres or Substructures)
-----------
Perform AAM
-----------

AAM using SMILES

  java -jar ReactionDecoder.jar -Q SMI -q "CC(O)CC(=O)OC(C)CC(O)=O.O[H]>>[H]OC(=O)CC(C)O.CC(O)CC(O)=O" -g -j AAM -f TEXT
Annotate Reaction using SMILES

  java -jar ReactionDecoder.jar -Q SMI -q "CC(O)CC(=O)OC(C)CC(O)=O.O[H]>>[H]OC(=O)CC(C)O.CC(O)CC(O)=O" -g -j ANNOTATE -f XML
-----------------
Compare Reactions
-----------------
Compare Reactions using SMILES with precomputed AAM mappings

  java -jar ReactionDecoder.jar -Q SMI -q -Q RXN -q example/ReactionDecoder_mapped.rxn  -T RXN -t example/ReactionDecoder_mapped.rxn -j COMPARE -f BOTH -u
Compare Reactions using RXN files

  java -jar ReactionDecoder.jar -Q RXN -q example/ReactionDecoder_mapped.rxn  -T RXN -t example/ReactionDecoder_mapped.rxn -j COMPARE -f BOTH
License

--------------------------
Reaction Decoder Tool (RDT)
--------------------------

Author: Syed Asad Rahman
e-mail: asad@ebi.ac.uk
c/o EMBL-European BioInformatics Institute (EBI)
WTGC, CB10 1SD Hinxton
UK

Note: The copyright of this software belongs to the author
and EMBL-European BioInformatics Institute (EBI).

--------
Citation
--------

References

Rahman, S.A. et.al.(2015) Reaction Decoder Tool (RTD): 
Extracting Features from Chemical Reactions (submitted)
