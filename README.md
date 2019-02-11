# Contact Divergence
Software related with the paper: "**Quantifying the evolutionary divergence of protein structures: the role of function change and function conservation**, Pascual-Garci­a et al _Proteins_ 2010 78:181-96"

* Program: Contact_divergence
* Author: Ugo Bastolla, (Centro de Biologia Molecular Severo Ochoa, CSIC-UAM, Spain)
* Email: <ubastolla@cbm.csic.es>
* Short description: Given a multiple sequence alignment (MSA) of proteins with known structures, for all aligned pairs it computes and prints sequence and structure similarity measures.
* License: Please see license file, and note that it includes the needlemanwunsch aligner developed by Dr. Andrew C. R. Martin in the Profit suite of programs, (c) SciTech Software 1993-2007.

USAGE
------

* COMPILE: make -f make_Contact_divergence
* RUN:     ./Contact_divergence <alignment file>
* EXAMPLE: ./Contact_divergence Input_Cont_Div_50044.aln

(you have to modify the names of the pdb files and directory in
Input_Cont_Div_50044.aln)

INPUT: 
------
Multiple sequence alignment of proteins in FASTA format.
The protein name is the name of a PDB file, optionally followed
by the chain index (Ex: >1opd.pdb A)
The first line may be PDBDIR=<directory of PDB files>
(default: current directory)

OUTPUT: 
-------
For each pair proteins, structural scores are printed:
* Contact divergence (Pascual-García et al Proteins 2010 78:181-96)
* Relatedness: 
 * 1 is CD is computed using the formula -log((q-q0)/(1-q0))
 * 0 if q is too small to be informative, and CD=(CD)_0-Z(q)
* Contact overlap (optional)
* TM score (Zhang & Skolnick Proteins 2004 57:702-10, currently not available)
* Sequence Indentity (optional)

OPTIONS:
--------
-cont Print contact overlap
-seq  Print sequence identity
