import sys
import pickle
sys.path.append('./library')

from Utility import *
from EnumIndexBuilder import *
from Match import *
from time import ctime

print 'XLSearch, version 1.0'
print 'Copyright of School of Informatics and Computing, Indiana University'
print 'Current time: %s' % ctime() 

if len(sys.argv) != 3:
	print 'Usage: python xlsearch.py paramter.txt output.txt'
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
print 'Current time: %s' % ctime()

print '\nBuilding index for all possible inter-peptide cross-links...'
index = EnumIndexBuilder(param['database'], specdict, mass, param)
del index.param['database']
del index.param['msdata']
del index.param['fix_mod_mass']
del index.param['fix_mod_res']
pickle.dump(index, file('index.pickle', 'w'))
index = pickle.load(file('index.pickle')) 
print 'Building index done!'
print 'Current time: %s' % ctime()

print '\nComputing features for candidate PSMs for query spectra...'
results = []
titles = []
for title in index.searchIndex.keys():
        if len(index.searchIndex[title]) != 0:
                titles.append(title)
length = len(titles)
print 'Total number of spectra to be searched: %d' % length 
for i in range(0, length):
        print '%d / %d' % (i, length)
        sys.stdout.flush()
        title = titles[i]
        result = getMatchesPerSpectrum(mass, param, index, title)
        result = (title, result)
        results.append(result)

print 'Computeing features done!\n'
print 'Current time: %s' % ctime()

pickle.dump(results, file('results.pickle', 'w'))
results = pickle.load(file('results.pickle'))

print '\nScoring and ranking PSMs...'
tophits = getTophits(index, results)
print 'Scoring and ranking PSMs done!'

print '\nOutputting results...'
pickle.dump(tophits, file('tophits.pickle', 'w'))
write_results(output_file, tophits)

if param['annotation'] == True:
	for i in range(len(tophits)):
		pep1 = tophits[i][0][0]
		pep2 = tophits[i][0][1]
		pro1 = tophits[i][2][0]
		pro2 = tophits[i][2][1]
		pos1 = tophits[i][1][0]
		pos2 = tophits[i][1][1]
		ch = tophits[i][3]
		title = tophits[i][5]
		filename = pep1 + '_' + pep2 + '_' + str(pos1) + '_' + str(pos2) + '_' + str(ch) + '_' + title + '.annotation'

		modification = dict(position=[], deltaMass=[])
		for j in range(len(pep1)):
			if pep1[j].islower():
				modification['position'].append(j)
				modification['deltaMass'].append(param['modMass'])
		pep1 = Peptide(pep1, ', '.join(pro1), modification, mass)

		modification = dict(position=[], deltaMass=[])
		for j in range(len(pep2)):
			if pep2[j].islower():
				modification['position'].append(j)
				modification['deltaMass'].append(param['modMass'])
		pep2 = Peptide(pep2, ', '.join(pro2), modification, mass)
		xl = XLink(pep1, pep2, [pos1, pos2], ch, mass, param)
		match = Match(specdict[title], xl, mass, param)
		match.annotate(filename, mass)

print 'Outputting results done!'
print 'XLSearch finished running!'
