import sys
import pickle
import os
import getopt
from time import ctime
import numpy as np

usage = '''
USAGE: python xlsearch_train.py -l [path to xlsearch library]
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
from fastareader import FastaReader

print 'XLSearch, version 1.0'
print 'Copyright of School of Informatics and Computing, Indiana University'
print 'Current time %s' % ctime()

print 'Training logistic regression models using authetic true-true PSMs...'

print '\nReading paramters from: %s...' % param_file
[param, mass] = read_param(param_file)
param['ntermxlink'] = False
param['neutral_loss']['h2o_loss']['aa'] = set('DEST')
param['neutral_loss']['nh3_loss']['aa'] = set('KNQR')
param['neutral_loss']['h2o_gain']['aa'] = set()
mass['C'] = 103.009184
print 'Reading parameters done!'

print '\nReading MSMS spectra files from directory: %s...' % param['ms_data']
spec_dict = read_spec(param['ms_data'], param, mass)
pickle.dump(spec_dict, file('spectra.pickle', 'w'))
print 'Total number of spectra: %d' % len(spec_dict)
print 'Reading MSMS spectra files done!'

print '\nDeisotoping MSMS spectra...'
spec_dict = pickle.load(file('spectra.pickle'))
deisotoped = dict()
titles = spec_dict.keys()
for i in range(len(titles)):
        title = titles[i]
        (one, align) = spec_dict[title].deisotope(mass, 4, 0.02)
        deisotoped[title] = one
pickle.dump(deisotoped, file('deisotoped.pickle', 'w'))
deisotoped = pickle.load(file('deisotoped.pickle'))
spec_dict = deisotoped
print 'Deisotoping MSMS spectra done!'
print 'Current time %s' % ctime()

print '\nBuilding index for all possible inter-peptide cross-links...'
index = EnumIndexBuilder(param['target_database'], spec_dict, mass, param)
pickle.dump(index, file('index.pickle', 'w'))
index = pickle.load(file('index.pickle'))
print 'Building index done!'
print 'Current time %s' % ctime()

print '\nComputing features for candidate PSMs for query spectra...'
results = []
titles = []
for title in index.search_index.keys():
	if len(index.search_index[title]) != 0:
		titles.append(title)
length = len(titles)
for i in range(0, length):
	print '%d / %d' % (i, length)
	sys.stdout.flush()
	title = titles[i]
	result = get_matches_per_spec(mass, param, index, title)
	result = [title, result]
	results.append(result)
print 'Computing features done!\n'
print 'Current time: %s' % ctime()

pickle.dump(results, file('results.pickle', 'w'))
results = pickle.load(file('results.pickle'))

print 'Extracting authentic true-true PSMs...'
true_true = get_true_true(results, index, param, mass)
pickle.dump(true_true, file('TT.pickle', 'w'))
print 'Extracting authentic true-true PSMs done!'

print 'Extracting true-false PSMs based on true-true PSMs as seeds...'
true_false = get_true_false(true_true, param, mass)
pickle.dump(true_false, file('TF.pickle', 'w'))
print 'Extracting true-false PSMs done!'

print 'Extracting false-false PSMs based on true-true PSMs as seeds...'
false_false = get_false_false(true_true, param, mass)
pickle.dump(false_false, file('FF.pickle', 'w'))
print 'Extracting false-false PSMs done!'

print 'Computing feature matrix for true-true, true-false, false-false PSMs...'
X_true_true = get_feature_matrix(true_true)
X_true_false = get_feature_matrix(true_false)
X_false_false = get_feature_matrix(false_false)

X_TT_TF = np.concatenate((X_true_true, X_true_false), axis = 0)
y_TT_TF = []
y_TT_TF.extend([1.0] * len(true_true))
y_TT_TF.extend([0.0] * len(true_false))
y_TT_TF = np.asarray(y_TT_TF)
y_TT_TF = y_TT_TF.T

X_TF_FF = np.concatenate((X_true_false, X_false_false), axis = 0)
y_TF_FF = []
y_TF_FF.extend([1.0] * len(true_false))
y_TF_FF.extend([0.0] * len(false_false))
y_TF_FF = np.asarray(y_TF_FF)
y_TF_FF = y_TF_FF.T
print 'Computing features done!'

from sklearn import linear_model
log_reg = linear_model.LogisticRegression()
log_reg.fit(X_TT_TF, y_TT_TF)
model_TT_TF = []
model_TT_TF.extend(log_reg.intercept_.tolist())
model_TT_TF.extend(log_reg.coef_.tolist())

log_reg = linear_model.LogisticRegression()
log_reg.fit(X_TF_FF, y_TF_FF)
model_TF_FF = []
model_TF_FF.extend(log_reg.intercept_.tolist())
model_TF_FF.extend(log_reg.coef_.tolist())

f = open(output_file, 'w')
f.write('# Classifier I (TT-TF) coefficients')
for i in range(len(model_TT_TF)):
	f.write('CI%02d\t')
	f.write('%.60f\n' % model_TT_TF[i])

f.write('# Classifier II (TF-FF) coefficients')
for i in range(len(model_TF_FF)):
	f.write('CII%02d\t')
	f.write('%.60f\n' % model_TF_FF[i])

f.write('nTT\t%d\n' % len(true_true))
f.write('nTF\t%d\n' % len(true_false))
f.write('nFF\t%d\n' % len(false_false))
f.close()
print 'XLSearch train mode finished!'
