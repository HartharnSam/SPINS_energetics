This set of scripts is to help the user interrogate and analyse the energetics-related outputs from SPINS. 

There are two toolboxes:
-Offline\_Sorting : This is to calculate the BPE, PE, KE, APE, etc offline using the output data. It implements a Winters et al (1995)-like sorting algorithm for a mapped case in SPINS, and helps the user analyse these
-Analyse\_Diagnos : This is an addition to the existing plot\_diagnos (part of SPINSmatlab). Instead, it reduces the temporal frequency consistently onto a regularly spaced grid across all outputs, and implements a different filter onto BPE which works for any size array. 

Author: Sam Hartharn-Evans, 2022

