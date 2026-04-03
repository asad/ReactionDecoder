Introduction
============

`Reaction Decoder Tool (RDT) v3.8.1`
--------------------------------------

`1. Atom Atom Mapping (AAM) Tool`

`2. Reaction Annotator (Extract Bond Changes, Identify & Mark Reaction Centres)`

`3. Reaction Comparator (Reaction Similarity based on Bond Changes, Reaction Centres or Substructures)`

Contact
============
Author: Dr. Syed Asad Rahman
e-mail: asad.rahman@bioinceptionlabs.com
Organisation: BioInception PVT LTD

Installation
============

`a)` [Download the latest RDT](https://github.com/asad/ReactionDecoder/releases) release from GitHub.

`b)` Compile using `maven`:

```
use pom.xml and mvn commands to build your project
1) mvn clean compile                                  (compile only)
2) mvn clean test                                     (compile and run tests)
3) mvn clean install -DskipTests=true                 (install, skip tests)
4) mvn clean install                                  (install with tests)
5) mvn -P local clean install -DskipTests=true        (fat jar, skip tests)
6) mvn -P local clean install                         (fat jar with tests)
```

Atom Atom Mapping — Simple Java API
=====================================

```java
import com.bioinceptionlabs.reactionblast.api.RDT;
import com.bioinceptionlabs.reactionblast.api.ReactionResult;

ReactionResult result = RDT.map("CC(=O)O.OCC>>CC(=O)OCC.O");
System.out.println("Mapped: " + result.getMappedSmiles());
System.out.println("Bond changes: " + result.getTotalBondChanges());
```

Atom Atom Mapping — Advanced CDK API
======================================

```java
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import com.bioinceptionlabs.reactionblast.mechanism.ReactionMechanismTool;
import com.bioinceptionlabs.reactionblast.tools.StandardizeReaction;

SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
IReaction rxn = sp.parseReactionSmiles("CC(=O)C=C.CC=CC=C>>CC1CC(CC=C1)C(C)=O");
rxn.setID("DielsAlder");

ReactionMechanismTool rmt = new ReactionMechanismTool(
        rxn, true, true, false, true, false, new StandardizeReaction());
System.out.println("Algorithm: " + rmt.getSelectedSolution().getAlgorithmID());
```

License
=======

`RDT` is released under the [GNU Lesser General Public License (LGPL) version 3.0](https://www.gnu.org/licenses/lgpl-3.0.en.html).

```
Author: Syed Asad Rahman
e-mail: asad.rahman@bioinceptionlabs.com
BioInception PVT LTD

Note: The copyright of this software belongs to the author and BioInception PVT LTD.
```

How to Cite RDT?
================

`SA Rahman, G Torrance, L Baldacci, SM Cuesta, F Fenninger, N Gopal, S Choudhary, JW May, GL Holliday, C Steinbeck and JM Thornton: Reaction Decoder Tool (RDT): Extracting Features from Chemical Reactions, Bioinformatics (2016)`

[doi: 10.1093/bioinformatics/btw096](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4920114/)

Sub-commands
============

`Perform AAM`
-------------

`AAM using SMILES`

```
java -jar rdt-3.8.1-jar-with-dependencies.jar -Q SMI -q "CC(O)CC(=O)OC(C)CC(O)=O.O[H]>>[H]OC(=O)CC(C)O.CC(O)CC(O)=O" -g -c -j AAM -f TEXT
```

`Perform AAM for Transporters` (accept mapping with no bond changes: `-b`)

```
java -jar rdt-3.8.1-jar-with-dependencies.jar -Q SMI -q "O=C(O)C(N)CC(=O)N.O=C(O)C(N)CS>>C(N)(CC(=O)N)C(=O)O.O=C(O)C(N)CS" -b -g -c -j AAM -f TEXT
```

`Annotate Reaction using SMILES`
---------------------------------

```
java -jar rdt-3.8.1-jar-with-dependencies.jar -Q SMI -q "CC(O)CC(=O)OC(C)CC(O)=O.O[H]>>[H]OC(=O)CC(C)O.CC(O)CC(O)=O" -g -c -j ANNOTATE -f XML
```

`Compare Reactions`
--------------------

`Compare using precomputed AAM mappings`

```
java -jar rdt-3.8.1-jar-with-dependencies.jar -Q RXN -q example/ReactionDecoder_mapped.rxn -T RXN -t example/ReactionDecoder_mapped.rxn -j COMPARE -f BOTH -u
```

`Compare using RXN files`

```
java -jar rdt-3.8.1-jar-with-dependencies.jar -Q RXN -q example/ReactionDecoder_mapped.rxn -T RXN -t example/ReactionDecoder_mapped.rxn -j COMPARE -f BOTH
```
