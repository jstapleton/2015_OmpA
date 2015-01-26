##Scripts used in the redesign of the OmpA surface

###OmpA_redesign123.cc
This is a C++ script that designs an OmpA surface using a Metropolis Monte Carlo search. This script was used to create Redesigns #1-3. Energies are computed with EzBeta and modified with a sequence entropy term. Without this term EzBeta will choose leucine, alanine, or phenylalanine at all positions. Note that the sequence entropy term discourages overrepresentation of any individual amino acid. As a result, the redesigned surfaces do not have as many aromatic residues as wild-type OmpA. Future versions will use a different method of encouraging surface diversity.

The version of this script used in the creation of Redesigns #1-3 contained a typo in one of the EzBeta energy calculations. This typo has been corrected in this version, but is pointed out in a comment in the script. 

###OmpA_redesign4.cc
This is a C++ script that designs an OmpA surface using a Metropolis Monte Carlo search. This script was used to create Redesign #4. This script differs from OmpA_redesign123.cc in the addition of a third energy term (in addition to EzBeta and sequence entropy) that promotes a surface that is energetically favorable when centered in the membrane but unfavorable when it is moved up or down in the membrane. 

The version of this script used in the creation of Redesign #4 contained a typo in one of the EzBeta energy calculations. This typo has been corrected in this version, but is pointed out in a comment in the script. 

###EzB_energy_calculator.cc
This is a C++ script that calculates the EzBeta energy of a given OmpA surface. The surface is entered as a hard-coded string of 43 surface residues.

#PDB files of redesigns

 - OR1.pdb
 - OR2.pdb
 - OR3.pdb
 - OR4.pdb

These were created by threading the redesigned sequences onto PBD 1bxw using Rosetta remodel with the supplied blueprint files and the RosettaRemodel config file (modified to take the appropriate blueprint file).

###Files to reproduce PDB threading:
- blueprint_R1.txt
- blueprint_R2.txt
- blueprint_R3.txt
- blueprint_R4.txt
- blueprint_R4cons.txt
- 1bxw.pdb
- RosettaRemodel.sh
To run, ensure the path to the Rosetta database is correct in RosettaRemodel.sh for your system, then run ./RosettaRemodel.sh on the command line.

# OmpA variant sequences
DNA and amino acid sequences of all variants described in the manuscript are provided in fasta format in all_sequences.fa.