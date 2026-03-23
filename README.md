Introduction
============

`Reaction Decoder Tool (RDT) v3.0.0`
--------------------------------------

`1. Atom Atom Mapping (AAM) Tool`

`2. Reaction Annotator (Extract Bond Changes, Identify & Mark Reaction Centres)`

`3. Reaction Comparator (Reaction Similarity based on the Bond Changes, Reaction Centres or Substructures)`

Contact
============
Author: Dr. Syed Asad Rahman
e-mail: asad.rahman@bioinceptionlabs.com

Installation
============

`a)` You could [download the latest RDT](https://github.com/asad/ReactionDecoder/releases) release version from the github.

`b)` Compile the core code using `maven`:

```
use pom.xml and mvn commands to build your project
1) mvn clean compile                                  (compile only)
2) mvn clean test                                     (compile and run tests)
3) mvn clean install -DskipTests=true                 (install, skip tests)
4) mvn clean install                                  (install with tests)
5) mvn -P local clean install -DskipTests=true        (fat jar, skip tests)
6) mvn -P local clean install                         (fat jar with tests)
```

Atom Atom Mapping using Java API
=================================

```java
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import com.bioinceptionlabs.reactionblast.mechanism.MappingSolution;
import com.bioinceptionlabs.reactionblast.mechanism.ReactionMechanismTool;
import com.bioinceptionlabs.reactionblast.tools.StandardizeReaction;

public class Example {
    public static void main(String[] args) throws Exception {
        final SmilesGenerator sg = new SmilesGenerator(SmiFlavor.AtomAtomMap);
        final SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());

        String reactionSM = "CC(=O)C=C.CC=CC=C>>CC1CC(CC=C1)C(C)=O";
        String reactionName = "DielsAlder";

        IReaction cdkReaction = smilesParser.parseReactionSmiles(reactionSM);

        IReaction performAtomAtomMapping = performAtomAtomMapping(cdkReaction, reactionName);
        System.out.println("AAM sm: " + sg.create(performAtomAtomMapping));
    }

    public static IReaction performAtomAtomMapping(IReaction cdkReaction, String reactionName) throws Exception {
        cdkReaction.setID(reactionName);
        boolean forceMapping = true;
        boolean generate2D = true;
        boolean generate3D = false;
        boolean complexMapping = true;
        boolean acceptNoChange = false;
        StandardizeReaction standardizeReaction = new StandardizeReaction();
        ReactionMechanismTool rmt = new ReactionMechanismTool(
                cdkReaction, forceMapping, generate2D, generate3D,
                complexMapping, acceptNoChange, standardizeReaction);
        MappingSolution s = rmt.getSelectedSolution();
        return s.getReaction();
    }
}
```


Migrating from v2.x
====================

The package namespace has changed from `uk.ac.ebi` to `com.bioinceptionlabs` in v3.0.0.

**Maven dependency**

```xml
<!-- Old (v2.x) -->
<groupId>uk.ac.ebi.rdt</groupId>

<!-- New (v3.0.0) -->
<groupId>com.bioinceptionlabs.rdt</groupId>
```

**Import changes**

Replace imports in your code:

| Old (v2.x) | New (v3.0.0) |
|-------------|--------------|
| `uk.ac.ebi.aamtool.*` | `com.bioinceptionlabs.aamtool.*` |
| `uk.ac.ebi.reactionblast.*` | `com.bioinceptionlabs.reactionblast.*` |
| `uk.ac.ebi.centres.*` | `com.bioinceptionlabs.centres.*` |

A simple find-and-replace of `uk.ac.ebi` with `com.bioinceptionlabs` in your import statements is sufficient. The API itself is unchanged.


License
=======

`RDT` is released under the [GNU Lesser General Public License (LGPL) version 3.0](https://www.gnu.org/licenses/lgpl-3.0.en.html).

```
Author: Syed Asad Rahman
e-mail: asad.rahman@bioinceptionlabs.com
BioInception

Note: The copyright of this software belongs to the author
and BioInception.
```

How to Cite RDT?
================

`SA Rahman, G Torrance, L Baldacci, SM Cuesta, F Fenninger, N Gopal, S Choudhary, JW May, GL Holliday, C Steinbeck and JM Thornton: Reaction Decoder Tool (RDT): Extracting Features from Chemical Reactions, Bioinformatics (2016)`

[doi: 10.1093/bioinformatics/btw096](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4920114/)


Sub-commands
===========


`Perform AAM`
-------------

`AAM using SMILES`

  ```
  java -jar rdt-3.0.0-jar-with-dependencies.jar -Q SMI -q "CC(O)CC(=O)OC(C)CC(O)=O.O[H]>>[H]OC(=O)CC(C)O.CC(O)CC(O)=O" -g -c -j AAM -f TEXT
  ```

`Perform AAM` for Transporters
-------------

`AAM using SMILES` (accept mapping with no bond changes -b)

  ```
  java -jar rdt-3.0.0-jar-with-dependencies.jar -Q SMI -q "O=C(O)C(N)CC(=O)N.O=C(O)C(N)CS>>C(N)(CC(=O)N)C(=O)O.O=C(O)C(N)CS" -b -g -c -j AAM -f TEXT
  ```

`Annotate Reaction using SMILES`
---------------------------------

  ```
  java -jar rdt-3.0.0-jar-with-dependencies.jar -Q SMI -q "CC(O)CC(=O)OC(C)CC(O)=O.O[H]>>[H]OC(=O)CC(C)O.CC(O)CC(O)=O" -g -c -j ANNOTATE -f XML
  ```


`Compare Reactions`
--------------------

`Compare Reactions using SMILES with precomputed AAM mappings`

  ```
  java -jar rdt-3.0.0-jar-with-dependencies.jar -Q RXN -q example/ReactionDecoder_mapped.rxn  -T RXN -t example/ReactionDecoder_mapped.rxn -j COMPARE -f BOTH -u
  ```


`Compare Reactions using RXN files`

  ```
  java -jar rdt-3.0.0-jar-with-dependencies.jar -Q RXN -q example/ReactionDecoder_mapped.rxn  -T RXN -t example/ReactionDecoder_mapped.rxn -j COMPARE -f BOTH
  ```
