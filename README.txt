XLSearch, Version 1.0
Copyright of School of Informatics and Computing, Indiana University
Contact: jic@indiana.edu, sujli@indiana.edu

I. INTRODUCTION
This software is intended to perform database sequence search for identifying 
chemically cross-linked peptide pairs from tandem mass spectra. Usage of this 
 software is free of charge for academic purposes.

II. PREREQUISITES

	i. Software packages

		This software can be run on Unix/Linux operating systems.
		1. Python version 2.6 or higher is required.
		2. To perform the in-sample training (i.e. 'training mode'), additional
		 python modules (Numpy 1.6.1 or higher, Scipy 0.9 or higher, Scikit-learn 
		 0.15 or higher) are required.

	ii. Data
		1. mzXML files containing tandem mass spectra converted using msconvert
		(http://proteowizard.sourceforge.net/tools.shtml) from RAW files.
		NOTE: Currently only mzXML format is supported.
		2. Fasta file containing the desired protein sequences to be searched
		against.

III. USAGE

	XLSearch can be run in 'searching mode' and 'training mode'. Searching
	mode is intended to perform the database sequence search where the peptide
	spectrum matches (PSM) are assigned a score based on the computed features 
	that describe the maching quality between spectrum and each individual peptide, 
	as well as weights of pre-trained logisitic models. Training mode is intended 
	to re-train the logistic models using authentic cross-link PSMs obtained from
	the new data.

	i. Searching mode
	Input:	1) PARAM.txt		Contains parameters for performing the database searching.
			2) database.fasta	Fasta format text file containing amino acid sequences 
								in fasta format. Specified in 'PARAM.txt'.
		
		Steps:
		1. Preparation:
			a. Unzip the .zip file to a directory (i.e. '/xlsearch_install_dir/'). It 
			should contain the python modules in '/xlsearch_install_dir/lib/', as well 
 			as the pipline script for searching and training model ('xlsearch_search.py' 
			and 'xlsearch_train.py').

			b. Create directory where search is to be performed (i.e. '/xlsearch_search_dir/').
			c. Copy the file 'xlsearch_search.py' and 'PARAM.txt' to this directory.
			d. Copy the fasta sequence file (i.e. 'database.fasta') to this directory.
			e. Create directory where the mzXML files are located (i.e. '/xlsearch_search_dir/mzxml/').
			f. Edit the parameter file 'PARAM.txt' as needed.

		2. Perform datbase search
			Under directory '/xlearch_search_dir/'

				$ python xlsearch_search.py -l /xlsearch_install_dir/
											-p PARAM.txt
											-o output.txt

			where '-l', '-p' and '-o' indicates the path to the library, parameter file
			and the output file name. All three arguments are required.

		3. Output file
			A tab-delimited text file contains top-scoring PSM for each query spectrum.
			Sorted by the joint probability score assigned to each PSM.
			The first line contains the headers of the columns:
				a. Rank of PSM
				b. Sequence of alpha peptide
				c. Sequence of beta peptide
				d. Index of cross-linking site on alpha
				e. Index of cross-linking site on beta
				f. Protein header of alpha peptide
				g. Protein header of beta peptide
				h. Charge state
				i. Joint probability score P(alpha = T, beta = T)
				j. Margianl probability P(alpha = T)
				k. Marginal probability P(beta = T)
				l. The title of the query spectrum

		4. Evalutating identified PSMs 
			The output file contains the top-scoring PSMs for each query spectrum sorted in descending
			order of the joint probability score. The percentage of false positive identification at a
			given score cutoff $S$ is estimated by counting the numbers of true-true, true-false, and
			false-false PSMs whose scores are greater than $S$. Specifically,

										FDR = (#(TF) - #(FF)) / #(TT)

			To filter the output PSMs at a given score cutoff, provide the value of 'cutoff' and 
			'is_unique' in the parameter file, where 'cutoff' indicates the desired fdr cutoff, 
			and 'is_unique' ('True' or 'False') indicates whether the unique cross-linked peptides 
			(i.e. the combination of cross-linked peptides and charge) or the redundant PSMs are counted 
			in the FDR calculation. For example, to filter for the results at 1% FDR cutoff where the 
			redundant PSMs are counted, set 'cutoff' to 0.01 and 'is_unique' to False.


			The filtered result will be written to file 'intra0.01.txt' and 'inter0.01txt' for intra-protein
			and inter-protein cross-links.

	ii. Training mode
	Input:	1) PARAM.txt        Contains parameters for performing the database searching.
			2) target_database.fasta	Contains only the TARGET sequences from which true-true
					PSMs can be identified.
			3) uniprot_database.fasta	Contains the pool of protein sequences from which the
				true-false and false-false PSMs can be generated based on the true-true PSMs.
			4) true_true.psm (Optional)	Contains the authentic true-true PSMs from which 
				the true-false and false-false PSMs can be genearted. Check the sample file for
				format.
			
		Steps:
		1. Preparation:
			a. Same as in searching mode.
			b. Create directory where training is to be performed (i.e. '/xlsearch_train_dir/')
			c. Copy 'xlsearch_train.py' to the current directory
			d. Copy the fasta sequence file ('target_database.fasta', 'uniprot_database.fasta')
				to the current directory
			e. Create directory where the mzXML files are located (i.e. '/xlsearch_search_dir/mzxml/').
			f. Edit the parameter file 'PARAM.txt' as needed.

		2. Perform training

			Under directory '/xlearch_train_dir/'
				$ python xlsearch_search.py -l /xlsearch_install_dir/
											-p PARAM.txt
											-o output.txt

			where '-l', '-p' indicates the path to the library, parameter file, and the output
			 file name. All three arguments are required.
			
		3. Output file
			The output will be in the following format:

			CI00	...	weight 0 of classifier I
			...
			CI15	...	weight 15 of classifier I

			CII00	... weight 0 of classifier II
			...
			CII15	... weight 15 of classifier II
			
			nTT		... number of true-true PSMs 
			nTF		... number of true-false PSMs
			nFF		...	number of false-false PSMs
			
			These lines correspond to the logistic regression parameters for classfier I and II 
			('CI' and 'CII'), and the numbers of true-true, true-false and false-false PSMs used 
			 to train them ('nTT', 'nTF', 'nFF'). The parameters in the 'PARAM.txt' can be 
			overwritten by these lines to use the updated model.

