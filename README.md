What is pes_nbo?
================

pes_nbo is Python script version of the java app NboScan (marpat/NboScan). It processes Gaussian output files containing results of the Potential Energy Scan (PES) jobs. NboScan extracts essential parts of the PES output file and creates formatted Gaussian input files (.gjf) to perform a single point NBO analysis. The number of files corresponds to the number of incremental steps in the coordinate scan.  Newly created .gjf files include NBO keywords and instructions to perform NBO analysis at each geometry on the PES.

*pes_nbo:*
-	Reads Gaussian output files containing geometries and energies of the PES scan (output of coordinate driving).
-	Generates multiple Gaussian input files (.gjf) containing custom route card, comments, geometries, and NBO keylists with the corresponding keywords.
-	Generates a plot view of energy profile along the scanned coordinate.
- Generates a *.dat file of energies and coordinates.
-	Creates Gaussian batch list file (.bcf) for processing multiple input (.gjf) files.

*How to run:*
Run from the terminal window (cmd) or shell as:
>> python pes_nbo3.py output_file.out