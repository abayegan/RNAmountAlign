Amir Bayegan, Peter Clote --Boston College

/******************************************************************************
 *   Copyright (C) 2018  Amir Bayegan, Peter Clote                            *
 *                                                                            *
 *  This program is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU General Public License as published by      *
 *  the Free Software Foundation, either version 3 of the License, or         *
 *  (at your option) any later version.                                       *
 *                                                                            *
 *  This program is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU General Public License for more details.                              *
 *                                                                            *
 *  You should have received a copy of the GNU General Public License         *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.     *
 ******************************************************************************/

RNAmountAlign performs a local/semi-global/global sequence structure alignment. 

RNAmountAlignScan executable is used for searching a query in a given target sequence.
 
Statistics are reported based on the alignment type either from Karlin-Altschul
or parameter fitting(see the paper).


Running the program:

Executables:
	RNAmountAlign
	RNAmountAlignScan

WARNING: The current libRNA.a is compiled for Linux systems. 
		 If there is a problem with linking the library, please download
		 and compile Vienna RNA Package on your system and replace 
		 the current libRNA.a with the one obtained from Vienna.

Usage: ./RNAmountAlign [options]
REQUIRED:
	RNAmountAlign:
		-f <string>	the input fasta file containing two sequences 
		OR
		-s <string> <string> provide sequence1 and sequence2


OPTIONS:

ALIGNMENT:
	-gi <float>	Gap initiation penalty for sequence alignment(Default:-3).
			A negative value should be provided

	-ge <float>	Gap extension penalty for sequence alignment(Default:-1).
			A negative value should be provided

	-gamma <float>	Weight of the structural homology. Must be in [0,1](Default:0.5)
			similarity = gamma*str_sim + (1-gamma)*seq_sim

	-m <string>	Similarity matrix file in RIBOSUM format.(Default:RIBOSUM85-60.mat)
			All RIBOSUM files are included in ./matrices directory

	-semi 		Perform semi-global alignment.
			Both gap ends of the firts sequence will be free of penalty.

	-local 		Perform local alignment

	-global 		Perform global alignment.(Default alignment type)
	
	-alifold		Output the consensus structure of the alignment from RNAalifold


STATISTICS:
		KA=Karlin-Altschul; EVD=extreme value dirstribution, ND=normal distribution

	-stat		Report statistics based on the alignment type(Default: off). 
			For local alignments: E-value from Karlin-Altschul(default) or EVD fitting
			For global alignments: p-value from ND fitting
			For semi-global alignments: p-value from ND fitting

	-evd		EVD fitting will be used for local alignments instead of KA

	-num <int>	Number of random sequences generated for fitting.(Default:500)

	-gc <int>	Size of GC bins (an integer between [0-100]). (Default: 5)
			This is used only with with parameter fitting and not KA.


OUTPUT:
	-o <string>	Write the output to a file.
			If not used the output will be printed to stdout.

	-format <clustal|fasta>	Format of the alignment output. (Default: clustal)

	-v		Verbose output. Prints MFE structures, the ensmeble expected and incremental heights.
	
	-h		Print help

/********************************************************************************
	 
Usage: ./RNAmountAlignScan [options]
REQUIRED:
	RNAmountAlignScan:
		-qf <string>	fasta file containing query sequence
		AND
		-tf <string>	fasta file containing target sequence


OPTIONS:

SLIDING
	 -window <int>	Size of the sliding window(Default:300).

	 -step <int>	Step size for incrementing the window start(Default:200).


ALIGNMENT:
	-gi <float>	Gap initiation penalty for sequence alignment(Default:-3).
			A negative value should be provided

	-ge <float>	Gap extension penalty for sequence alignment(Default:-1).
			A negative value should be provided

	-gamma <float>	Weight of the structural homology. Must be in [0,1](Default:0.5)
			similarity = gamma*str_sim + (1-gamma)*seq_sim

	-m <string>	Similarity matrix file in RIBOSUM format.(Default:RIBOSUM85-60.mat)
			All RIBOSUM files are included in ./matrices directory

	-semi 		Perform semi-global alignment.
			Both gap ends of the firts sequence will be free of penalty.

	-local 		Perform local alignment

	-global 		Perform global alignment.(Default alignment type)


STATISTICS:
		KA=Karlin-Altschul; EVD=extreme value dirstribution, ND=normal distribution

	-stat		Report statistics based on the alignment type(Default: off). 
			For local alignments: E-value from Karlin-Altschul or EVD fitting
			For global alignments: p-value from ND fitting
			For semi-global alignments: p-value from ND fitting

	-evd		EVD fitting will be used for local alignments instead of KA

	-num <int>	Number of random sequences generated for fitting.(Default:500)

	-gc <int>	Size of GC bins (an integer between [0-100]). (Default: 5)
			This is used only with with parameter fitting and not KA.


OUTPUT:
	-o <string>	Write the output to a file.
			If not used the output will be printed to stdout.

	-format <clustal|fasta>	Format of the alignment output. (Default: clustal)

	-v		Verbose output. Prints MFE structures, the ensmeble expected and incremental heights.
	-h		Print help
