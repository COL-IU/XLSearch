XLSearch, Version 1.0
Copyright of School of Informatics and Computing, Indiana University
Contact: jic@indiana.edu, sujunli@indiana.edu

This software is intended to perform database sequence search for inter-peptide cross-links on MS/MS spectra.
Usage of this software is free of charge for academic purposes.

The user can perform the database search in the searching mode which uses pre-trained classification models. 
Alternatively, the user can perform in-sample training (i.e. training mode) and update the coefficients of the
classification models.

Searching mode:

Input:	1)	parameter.txt	-	Contains parameters for database searching
	2)	protein sequence database	-	Text file containing amino acid sequences in fasta format. Specified in 'parameter.txt'
	3)	MS/MS spectra	-	All MS/MS spectra files should be moved to a single directory. Currently only mzXML files are supported. Specified in 'paramter.txst'
Output:	1)	output.txt	-	Text file containing top hit of each query spectrum

Usage:	'python xlsearch.py paramter.txt output.txt'

Training mode:

Usage:	'python train_mode.py parameters.txt output.txt'

Note: The coefficients of logistic regression models will be shown in the standard ouput. They can be copy and
pasted to the paramter.txt file, and database search can be performed using the updated paramters.
