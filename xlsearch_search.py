import sys
import pickle
import getopt
from time import ctime

usage = '''
USAGE: python xlsearch_search.py -l [path to xlsearch library]
                                 -p [parameter file]
                                 -o [output file]'''

(pairs, args) = getopt.getopt(sys.argv[1:], 'l:p:o:')
cmd_arg = dict()
for i in range(len(pairs)):
	cmd_arg[pairs[i][0]] = pairs[i][1]

if len(cmd_arg) != 3:
	print usage
	sys.exit(1)

lib_path = cmd_arg['-l']
param_file = cmd_arg['-p']
output_file = cmd_arg['-o']

sys.path.append(lib_path)

from utility import *
from index import EnumIndexBuilder
from match import Match

print 'XLSearch, version 1.1'
print 'Copyright of School of Informatics and Computing, Indiana University'
print 'Current time: %s' % ctime() 

print 'Performing database search for identifying cross-linked peptides...'

print '\nReading paramters from: %s...' % param_file
[param, mass] = read_param(param_file)
print 'Reading parameters done!'

print '\nReading MSMS spectra files from directory: %s...' % param['ms_data']
spec_dict = read_spec(param['ms_data'], param, mass)
print 'Total number of spectra: %d' % len(spec_dict)
print 'Reading MSMS spectra files done!'

if param['deisotope'] == True:
	print '\nDeisotoping MSMS spectra...'
	deisotoped = dict()
	titles = spec_dict.keys()

	ms2tol = 0.02 # 0.02Da for deisotoping

	for i in range(len(titles)):
		title = titles[i]
		(one, align) = spec_dict[title].deisotope(mass, param['ndeisotope'], ms2tol)
		deisotoped[title] = one

	spec_dict = deisotoped

	print 'Deisotoping MSMS spectra done!'
	print 'Current time: %s' % ctime()

print '\nBuilding index for all possible inter-peptide cross-links...'
index = EnumIndexBuilder(param['database'], spec_dict, mass, param)
print 'Building index done!'
print 'Current time: %s' % ctime()

print '\nComputing features for candidate PSMs for query spectra...'
results = []
titles = spec_dict.keys()

length = len(titles)
print 'Total number of spectra to be searched: %d' % length
for i in range(0, length):
	print '%d / %d' % (i, length)
	sys.stdout.flush()
	title = titles[i]
	result = get_matches_per_spec(mass, param, index, title)
	if len(result) > 0:	
		result = [title, result]
		results.append(result)
print 'Computing features done!\n'
print 'Current time: %s' % ctime()

print '\nScoring and ranking PSMs corresponding to the same query spectrum...'
print 'Ranking all top-hit PSMs...'
tophits = get_tophits(index, results)
print 'Scoring and ranking PSMs done!'
print 'Current time: %s' % ctime()

print '\nOutputting results...'
pickle.dump(tophits, file('tophits.pickle', 'w'))
print 'Writing result to file: %s' % output_file
write_results(output_file, tophits)
print 'Outputtting results done!'
print 'Current time %s' % ctime()

if param['annotation'] == True:
	for i in range(10): #range(len(tophits)):
		[match, filename] = get_matches_from_tophits(tophits[i], param, mass, spec_dict)
		match.annotate(filename, mass)

filter_by_fdr(tophits, param['decoy_string'], param['cutoff'], param['is_unique'])

print 'XLSearch finished running!'
