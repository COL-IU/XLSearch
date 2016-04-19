import sys
import pickle
import os
sys.path.append('./library')


from Utility import *
from EnumIndexBuilder import *
from FastaReader import *

print 'XLSearch, version 1.0'
print 'Copyright of School of Informatics and Computing, Indiana University'

if len(sys.argv) != 3:
	print 'Usage: python train_model.py paramter.txt output.txt'
else:
	para_file = sys.argv[1]
	output_file = sys.argv[2]

print '\nReading paramters from: %s...' % para_file
[param, mass] = readParam(para_file)
print 'Reading parameters done!'

print '\nReading MSMS spectra files from directory: %s...' % param['msdata']
specdict = readSpectra(param['msdata'], param, mass)
pickle.dump(specdict, file('spectra.pickle', 'w'))
print 'Total number of spectra: %d' % len(specdict)
print 'Reading MSMS spectra files done!'

print '\nDeisotoping MSMS spectra...'
specdict = pickle.load(file('spectra.pickle'))
deisotoped = dict()
titles = specdict.keys()
for i in range(len(titles)):
        title = titles[i]
        (one, align) = specdict[title].deisotope(mass, 4, 0.02)
        deisotoped[title] = one
pickle.dump(deisotoped, file('deisotoped.pickle', 'w'))
deisotoped = pickle.load(file('deisotoped.pickle'))
specdict = deisotoped
print 'Deisotoping MSMS spectra done!'

print '\nBuilding index for all possible inter-peptide cross-links...'
param['neutralloss']['h2oLoss']['aa'] = set('DEST')
param['neutralloss']['nh3Loss']['aa'] = set('KNQR')
param['neutralloss']['h2oGain']['aa'] = set()
mass['C'] = 103.009184

index = EnumIndexBuilder(param['target_database'], specdict, mass, param)
pickle.dump(index, file('trainindex.pickle', 'w'))
index = pickle.load(file('trainindex.pickle'))

results = []
titles = []
for title in index.searchIndex.keys():
	if len(index.searchIndex[title]) != 0:
		titles.append(title)

length = len(titles)
for i in range(0, length):
	print i
	sys.stdout.flush()
	title = titles[i]
	result = getMatchesPerSpectrum(mass, param, index, title)
	result = (title, result)
	results.append(result)

#pickle.dump(results, file('trainresults.pickle', 'w'))

TT = get_true_true(results, index, param, mass)
AB = get_true_false(TT, param, mass)
FF = get_false_false(TT, param, mass)

[xTT, xTF, xFF] = calculate_feature(TT, AB, FF)

f = open('tt', 'w')
for i in range(len(xTT)):
	for j in range(16):
		f.write('%f\t' % xTT[i][j])

	f.write('%s\n' % '')
f.close()

f = open('tf', 'w')
for i in range(len(xTF)):
	for j in range(16):
		f.write('%f\t' % xTF[i][j])

	f.write('%s\n' % '')
f.close()

f = open('ff', 'w')
for i in range(len(xFF)):
	for j in range(16):
		f.write('%f\t' % xFF[i][j])

        f.write('%s\n' % '')
f.close()
pickle.dump([xTT, xTF, xFF], file('X.pickle', 'w'))

print 'Classifier I coefficients:'
os.system('python lg_r.py tt tf')
print 'Classifier II coefficients:'
os.system('python lg_r.py tf ff')

print 'nTT = %d' % len(xTT)
print 'nTF = %d' % len(xTF)
print 'nFF = %d' % len(xFF)
