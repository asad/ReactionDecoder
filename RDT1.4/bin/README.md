Introduction
============

`Reaction Decoder Tool (RDT)`
-----------------------------

`1. Atom Atom Mapping (AAM) Tool`

`2. Reaction Annotator (Extract Bond Changes, Identify & Mark Reaction Centres) and `

`3. Reaction Comparator (Reaction Similarity based on the Bond Changes, Reaction Centres or Substructures)`


Installation
============

`a)` You could [download the latest RDT] (https://github.com/asad/ReactionDecoder/releases) release version from the github.

`b)` Compile the code using `ant`?:

`Ant Build` commands

```
CLEAN:
  ant clean
BUILD:
  ant compile
DIST:
  ant jar
DOC:
  ant javadoc
TEST:
  ant test
HELP:
  ant run
Fat Jar:
 ant package-for-store
```


Licence
=======

`RDT` is released under the [GNU General Public License version 3](http://www.gnu.org/licenses/gpl.html).

```
Author: Syed Asad Rahman
e-mail: asad@ebi.ac.uk
c/o EMBL-European BioInformatics Institute (EBI)
WTGC, CB10 1SD Hinxton
UK

Note: The copyright of this software belongs to the author
and EMBL-European BioInformatics Institute (EBI).
```

How to Cite RDT?
================

`SA Rahman, G Torrance, L Baldacci, SM Cuesta, F Fenninger, N Gopal, S Choudhary, JW May, GL Holliday, C Steinbeck and JM Thornton: Reaction Decoder Tool (RDT): Extracting Features from Chemical Reactions, Bioinformatics (2016), doi: 10.1093/bioinformatics/btw096`


Subcommands
===========


`Perform AAM`
-------------

`AAM using SMILES`
  
  ```
  java -jar RDT1.4.1.jar -Q SMI -q "CC(O)CC(=O)OC(C)CC(O)=O.O[H]>>[H]OC(=O)CC(C)O.CC(O)CC(O)=O" -g -j AAM -f TEXT
  ```

  ```
  java -cp dist/*:lib/* aamtool.ReactionDecoder -Q SMI -q "CC(O)CC(=O)OC(C)CC(O)=O.O[H]>>[H]OC(=O)CC(C)O.CC(O)CC(O)=O" -g -j  AAM -f TEXT
  ```

`Annotate Reaction using SMILES`
---------------------------------

  ```
  java -jar RDT1.4.1.jar -Q SMI -q "CC(O)CC(=O)OC(C)CC(O)=O.O[H]>>[H]OC(=O)CC(C)O.CC(O)CC(O)=O" -g -j ANNOTATE -f XML
  ```


`Compare Reactions`
--------------------

`Compare Reactions using SMILES with precomputed AAM mappings`
  
  ```
  java -jar RDT1.4.1.jar -Q RXN -q example/ReactionDecoder_mapped.rxn  -T RXN -t example/ReactionDecoder_mapped.rxn -j COMPARE -f BOTH -u
  ```


`Compare Reactions using RXN files`
  
  ```
  java -jar RDT1.4.1.jar -Q RXN -q example/ReactionDecoder_mapped.rxn  -T RXN -t example/ReactionDecoder_mapped.rxn -j COMPARE -f BOTH
  ```
