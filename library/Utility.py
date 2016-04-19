import re
import math
import glob
import pickle
import sys

from Peptide import *
from XLink import *
from Match import *
from MZXMLReader import *
from FastaReader import *
from random import shuffle

def logisticEval(b, x):
	tmp = x
	x = [1]
	x.extend(tmp)
#	b = [1].extend(b)
	val = 0
#	print '%d\t%d\n' % (len(b), len(x))

	for i in range(len(b)):
		val += b[i] * x[i]
	
	return float(1) / (1 + math.exp(-val))

def adjustPrior(prior_t, p, iter):

	posterior_t = []
	for unit in p:
		posterior_t.extend(unit)

	N = len(posterior_t)
	prior = prior_t
	posterior = [0.0] * N 

	for i in range(iter):
		for k in range(N):
			denom = (float(prior) / prior_t) * posterior_t[k] + (float(1 - prior) / (1 - prior_t)) * (1 - posterior_t[k])

			if denom == 0:
				print 'denom == 0'
			posterior[k] = float(prior * posterior_t[k]) / (prior_t * denom)

		prior_prev = prior
		prior = sum(posterior) / len(posterior)
#		print prior
	newp = []

	b = 0
	for i in range(len(p)):
		e = b + len(p[i])
		newp.append(posterior[b : e])
		b = e

	return newp
def adjustPriorMarginal(prior_t, p, marginal, iter):

	posterior_t = []
	for unit in p:
		posterior_t.extend(unit)

	N = len(posterior_t)
	prior = prior_t
	posterior = [0.0] * N
#	marginal = [1.0] * N
	for i in range(iter):
		for k in range(N):
			denom = (float(prior) / prior_t) * posterior_t[k] + (float(1 - prior) / (1 - prior_t)) * (1 - posterior_t[k])

			if denom == 0:
				print 'denom == 0'
			posterior[k] = float(prior * posterior_t[k]) / (prior_t * denom)

		prior_prev = prior
#		prior = sum(posterior) / len(posterior)
		prior = 0.0
		for j in range(N):
			prior = prior + posterior[j] * marginal[j]
		prior = prior / sum(marginal)

		print 'iteration = %d, %.10f' % (i, prior)
		if i > 0 and abs(float(prior_prev - prior) / prior_prev) <= 0.01:
			print '%f' % abs(float(prior_prev - prior) / prior_prev)
			break
	newp = []

	b = 0
	for i in range(len(p)):
		e = b + len(p[i])
		newp.append(posterior[b : e])
		b = e

	return newp

def getMarginal(p21, p11, p12, p22):
	alphaT = []
	betaT = []
	alphaF =  []
	betaF = []

	for i in range(len(p21)):
		for j in range(len(p21[i])):
			denom = 1.0 - (p11[i][j] - p12[i][j]) * (p21[i][j] - p22[i][j])

			at = (p12[i][j] + p22[i][j] * (p11[i][j] - p12[i][j])) / denom
			bt = (p22[i][j] + p12[i][j] * (p21[i][j] - p22[i][j])) / denom
			af = 1.0 - at
			bf = 1.0 - bt
			alphaT.append(at)
			betaT.append(bt)
			alphaF.append(af)
			betaF.append(bf)

	return [alphaT, betaT, alphaF, betaF]

def getMatchesPerSpectrum(mass, param, index, title):
	spectraDict = index.spectraDict
	uniquePepObjs = index.uniquePepObjs[0]
	precMassPepIndexTuple = index.precMassPepIndexTuple
	searchIndex = index.searchIndex
	xresidue = param['xresidue']

	indexList = searchIndex[title]
	spectrum = spectraDict[title]

	matches = []
	XL = []
	for i in indexList:
		index1 = precMassPepIndexTuple[1][i]
		index2 = precMassPepIndexTuple[2][i]
		
		pepObj1 = uniquePepObjs[index1]
		pepObj2 = uniquePepObjs[index2]

		pepSorted = sorted([pepObj1, pepObj2], key = lambda x : x.sequence)
		pepObj1 = pepSorted[0]
		pepObj2 = pepSorted[1]

		charge = spectraDict[title].ch

		mz = spectraDict[title].mz
		it = spectraDict[title].it

		kPos1 = []
		kPos2 = []
		if param['ntermxlink'] == True:
			if pepObj1.isNterm == True:
				kPos1.append(0)
			if pepObj2.isNterm == True:
				kPos2.append(0)

		pepseq1 = pepObj1.sequence
		kPos1.extend(list(zip(*filter(lambda x : x[1] == xresidue, enumerate(pepseq1[:-1])))[0]))
		pepseq2 = pepObj2.sequence
		kPos2.extend(list(zip(*filter(lambda x : x[1] == xresidue, enumerate(pepseq2[:-1])))[0]))

		for p1 in kPos1:
			for p2 in kPos2:
				positions = [p1, p2]
				xl = XLink(pepObj1, pepObj2, positions, charge, mass, param);	XL.append(xl)
			
				match = Match(spectrum, xl, mass, param)
				match.match(mass)
				matches.append(match.getMatchInfo(index))

	return matches




def getPepObjsFromProtein(header, sequence, patternString, mass, param):
	missedSites = param['missedsites']
	minLen = param['minlength']
	maxLen = param['maxlength']
	modRes = param['modRes']

	pattern = re.compile(patternString)

	sites = [0]
	for i in range(len(sequence)):
		if i == len(sequence) - 1:
			sites.append(i + 1)
		elif (sequence[i] == 'K' or sequence[i] == 'R') and sequence[i + 1] != 'P':
			sites.append(i + 1)

	peptides = []
	for i in range(len(sites)):
		if i < len(sites) - missedSites - 1:
			for j in range(missedSites + 1):
				seq = sequence[sites[i] : sites[i + j + 1]]
				if len(seq) >= minLen and len(seq) <= maxLen and pattern.match(seq):
					peptides.append(seq)

		else:
			for j in range(i + 1, len(sites)):
				seq = sequence[sites[i] : sites[j]]
				if len(seq) >= minLen and len(seq) <= maxLen and pattern.match(seq):
					peptides.append(seq)

	peptides = list(set(peptides))

	pepObjs = []
	for pep in peptides:
		modification = dict(position=[], deltaMass=[])

		isNterm = False
		if pep == sequence[:len(pep)]:
			isNterm = True
		pepObjs.append(Peptide(pep, header, modification, mass, isNterm))

		if len(modRes) != 0:
			modMass = param['modMass']
			index = [i for i, ltr in enumerate(pep) if ltr == modRes]

			if len(index) != 0:
				pep = list(pep)
				for i in index:
					modification = dict(position=[], deltaMass=[])
					modification['position'].append(i)
					modification['deltaMass'].append(modMass)

					pep[i] = pep[i].lower()
				pep = ''.join(pep)

				pepObjs.append(Peptide(pep, header, modification, mass, isNterm))

	return pepObjs

def readParam(filename):
	param = dict(
		useAIon=True,
		verbose=False,
		chargePreXlinkIons=[1, 3],
		chargePostXlinkIons=[2, 5],
		basepeakint = 100.0,
		dynamicrange = 0.001,
		missedsites = 2,
		minlength = 4,
		maxlength = 51,
		modRes = '',
		modMass = 0.0,
		linkermass = 136.10005,
		ms1tol = dict(measure='ppm', val=5),
		ms2tol = dict(measure='da', val=0.01),
		minmz = 200,
		maxmz = 2000,
		mode = 'conservative',
		xresidue = 'K',
#		patternstring = '^[ACDEFGHIKLMNPQRSTVWY]*K[ACDEFGHIKLMNPQRSTVWY]+$',
		aa = 'ACDEFGHIKLMNPQRSTVWY',
		neutralloss=dict(
			h2oLoss=dict(
				mass=-18.010565,
				aa=set('ACDEFGHIKLMNPQRSTVWY')),
			nh3Loss=dict(
				mass=-17.026549,
				aa=set('ACDEFGHIKLMNPQRSTVWY')),
			h2oGain=dict(
				mass=18.010565,
				aa=set('ACDEFGHIKLMNPQRSTVWY'))),
		model_TT_TF = [0.0] * 17,
		model_TF_FF = [0.0] * 17,
		nTT = 169,
		nTF = 8568,
		nFF = 91242)
		
	mass = dict(
		A=71.037114,
		R=156.101111,
		N=114.042927,
		D=115.026943,
		C=103.009184,
		E=129.042593,
		Q=128.058578,
		G=57.021464,
		H=137.058912,
		I=113.084064,
		L=113.084064,
		K=128.094963,
		M=131.040485,
		F=147.068414,
		P=97.052764,
		S=87.032028,
		T=101.047678,
		W=186.079313,
		Y=163.063329,
		V=99.068414,
		Hatom=1.007825032,
		Oatom=15.99491462,
		neutronmass = 1.008701,
		BIonRes=1.0078246,
		AIonRes=-26.9870904,
		YIonRes=19.0183888,
		isotopeInc = [1.008701/4, 1.008701/3, 1.008701/2, 1.008701/1])

	f = open(filename)	
	lines = f.readlines()
	for l in lines:
		l = l[:-1]
		columns= l.split('\t')
		if len(l) == 0 or l[0] == '#' or len(columns) < 2:
			continue

		name = columns[0]
		val = columns[1]
		if name == 'database':
			param['database'] = val
		elif name == 'MS_data_directory':
			param['msdata'] = val
		elif name == 'XLresidue':
			param['xresidue'] = val
		elif name == 'ms1tol_unit':
			param['ms1tol']['measure'] = val
		elif name == 'ms1tol_val':
			param['ms1tol']['val'] = int(val)
		elif name == 'ms2tol_unit':
			param['ms2tol']['measure'] = val
		elif name == 'ms2tol_val':
			param['ms2tol']['val'] = float(val)
		elif name == 'linker_mass':
			param['linkermass'] = float(val)
		elif name == 'miss_cleave':
			param['missedsites'] = int(val)
		elif name == 'include_a_ions':
			param['useAIon'] = True if val == 'True' else False 
		elif name == 'min_peplen':
			param['minlength'] = int(val)
		elif name == 'max_peplen':
			param['maxlength'] = int(val)
		elif name == 'fix_mod_res':
			param['fix_mod_res'] = val
		elif name == 'fix_mod_mass':
			param['fix_mod_mass'] = float(val)
		elif name == 'var_mod_res':
			param['modRes'] = val
		elif name == 'var_mod_mass':
			param['modMass'] = float(val)
		elif name == 'min_preXL_ions_ch':
			param['chargePreXlinkIons'][0] = int(val)
		elif name == 'max_preXL_ions_ch':
			param['chargePreXlinkIons'][1] = int(val)
		elif name == 'min_postXL_ions_ch':
			param['chargePostXlinkIons'][0] = int(val)
		elif name == 'max_postXL_ions_ch':
			param['chargePostXlinkIons'][1] = int(val)
		elif name == 'target_database':
			param['target_database'] = val
		elif name =='uniprot_database':
			param['uniprot_database'] = val
		elif name == 'max_iterations':
			param['MAX_ITERATIONS'] = int(val)
		elif name == 'annotate_spec':
			param['annotation'] = True if val == 'True' else False
		elif name == 'ntermxlink':
			param['ntermxlink'] = True if val == 'True' else False
		elif len(name) >= 4 and name[:2] == 'CI':
			if len(name) == 4:
				s = int(name[2:])
				param['model_TT_TF'][s] = float(val)
			elif len(name) == 5:
				s = int(name[3:])
				param['model_TF_FF'][s] = float(val)
		elif name == 'nTT':
			param['nTT'] = int(val)
		elif name == 'nTF':
			param['nTF'] = int(val)
		elif name == 'nFF':
			param['nFF'] = int(val)

	param['patternstring'] = '^[' + param['aa'] + ']*' + param['xresidue'] + '[' + param['aa'] + ']+$'
	param['prior_t_TT_TF'] = float(param['nTT']) / (param['nTT'] + param['nTF'])
	param['prior_t_TF_FF'] = float(param['nTF']) / (param['nTF'] + param['nFF'])

	f.close()
	mass[param['fix_mod_res']] += param['fix_mod_mass']

	return [param, mass]

def readSpectra(directory, param, mass):
	files = glob.glob(directory + '*.mzXML')

	specdict = dict()
	total = []
	for filename in files:
		reader = MZXMLReader(filename)
		spectra = reader.getSpectraList(mass, param)
		total.append(spectra)

	ss = []
	for i in range(len(total)):
		ss.append(set())

		tmp = []
		for j in range(len(total[i])):
			if total[i][j].rt >= 0 and total[i][j].rt <= 110*60 and total[i][j].ch >= 2 and total[i][j].ch <= 7:
				tmp.append(total[i][j])

		tmp = sorted(tmp, key = lambda s : s.mr)

		tolerance = 0.01
		lower_ratio = 0.3
		upper_ratio = 1 / float(lower_ratio)

		for j in range(len(tmp) - 1):
			MZ = []
			IT = []
			mz = tmp[j].mz
			it = tmp[j].it
			lastindex = 0
			ik = 0
			jk = 0
			for ik in range(len(mz)):
				if lastindex == 0:
					jk = 0
				else:
					jk = lastindex
				while not (lastindex > len(mz) - 1 or jk > len(mz) - 1 or ik > len(mz) - 1 or mz[jk] > mz[ik] + tolerance):
					if mz[jk] <= mz[ik] - tolerance:
						lastindex = jk

					ratio = float(it[ik]) / float(it[jk])
					if abs(mz[ik] - mz[jk]) <= tolerance and ratio >= lower_ratio and ratio <= upper_ratio:
						MZ.append(mz[ik])
						IT.append(it[ik])
					jk = jk + 1	
			if len(MZ) >= 25:
				specdict[tmp[j].title] = tmp[j]
	return specdict
	# deisotope spectra
#	deisotoped = dict()
#	titles = specdict.keys()
#	for i in range(len(titles)):
#		title = titles[i]
#		(one, align) = specdict[title].deisotope(mass, 4, 0.02)
#		deisotoped[title] = one
#	return deisotoped


def getTophits(index, result):

	model_TT_TF = index.param['model_TT_TF']
	model_TF_FF = index.param['model_TF_FF']

	prior_t_TT_TF = index.param['prior_t_TT_TF']
	prior_t_TF_FF = index.param['prior_t_TF_FF']
	it = index.param['MAX_ITERATIONS']

	p21 = []
	p11 = []
	p12 = []
	p22 = []

	for i in range(len(result)):
		print i

#		title = result[i][0]
#		ch = str(title.split('.')[-1])

		p21.append([])
		p11.append([])
		p12.append([])
		p22.append([])
#		S = []

		for j in range(len(result[i][1])):

#			pep1 = index.uniquePepObjs[0][result[i][1][j][0][0]]
#			pep2 = index.uniquePepObjs[0][result[i][1][j][0][1]]

#			positions = result[i][1][j][1]
			feature = result[i][1][j][2]

			x = list(feature[0])
			x.extend(feature[1])

			xflip = list(feature[1])
			xflip.extend(feature[0])

#			s = [pep1.sequence, pep2.sequence, str(positions[0] + 1), str(positions[1] + 1)]
#			s = '_'.join(s)
#			s = s + '_' + ch + '_' + title

#			S.append(s)

			b = model_TT_TF

			p21[-1].append(logisticEval(b, x))
			p11[-1].append(logisticEval(b, xflip))

			b = model_TF_FF

			p12[-1].append(logisticEval(b, x))
			p22[-1].append(logisticEval(b, xflip))

	[alphaT, betaT, alphaF, betaF] = getMarginal(p21, p11, p12, p22)
	p21 = adjustPriorMarginal(prior_t_TT_TF, p21, alphaT, it)
	p11 = adjustPriorMarginal(prior_t_TT_TF, p11, betaT, it)
	p12 = adjustPriorMarginal(prior_t_TF_FF, p12, betaF, it)
	p22 = adjustPriorMarginal(prior_t_TF_FF, p22, alphaF, it)

	for i in range(len(result)):
		print i
		result[i] = list(result[i])

		for j in range(len(result[i][1])):

			pep1 = index.uniquePepObjs[0][result[i][1][j][0][0]]
			pep2 = index.uniquePepObjs[0][result[i][1][j][0][1]]

			ap21 = p21[i][j]
			ap11 = p11[i][j]
			ap12 = p12[i][j]
			ap22 = p22[i][j]

			denom = 1 - (ap11 - ap12) * (ap21 - ap22)

			alaph_T = (ap12 + ap22 * (ap11 - ap12)) / denom
			beta_T = (ap22 + ap12 * (ap21 - ap22)) / denom
			prob1 = ap11 * beta_T
			prob2 = ap21 * alaph_T
			score = (prob1 + prob2) / float(2)

			info = {'alpha' : alaph_T, 'beta' : beta_T, 'prob' : [prob1, prob2], 'score' : score}

			result[i][1][j] = list(result[i][1][j])
			result[i][1][j].append(info)

	for r in result:
		r[1] = sorted(r[1], key = lambda x : x[3]['score'], reverse = True)

	result = sorted(result, key = lambda x : x[1][0][3]['score'], reverse = True)

	tophits = []

	for r in result:
		scan = r[0]
		pep = [index.uniquePepObjs[0][r[1][0][0][0]].sequence, index.uniquePepObjs[0][r[1][0][0][1]].sequence]
		pos = [int(r[1][0][1][0]), int(r[1][0][1][1])]
		pro = [index.uniquePepObjs[0][r[1][0][0][0]].proteinID, index.uniquePepObjs[0][r[1][0][0][1]].proteinID]
		ch = int(scan.split('.')[-1])
		score = r[1][0][3]['score']
		tophits.append([pep, pos, pro, ch, score, scan])

	return tophits

def write_results(output_file, tophits):
	f = open(output_file, 'w')
	f.write('Rank\tPep_alpha\tPep_beta\tSite_alpha\tSite_beta\tPro_alpha\tPro_beta\tCharge\tpr(alpha=T,beta=T)\tSpectrum\n')
	for i in range(len(tophits)):
		f.write('%d\t' % (i + 1))
		f.write('%s\t%s\t' % (tophits[i][0][0], tophits[i][0][1]))
		f.write('%d\t%d\t' % (tophits[i][1][0], tophits[i][1][1]))
		f.write('%s\t%s\t' % (','.join(tophits[i][2][0]), ','.join(tophits[i][2][1])))
		f.write('%d\t' % tophits[i][3])
		f.write('%E\t' % tophits[i][4])
		f.write('%s\n' % tophits[i][5])
	f.close()

def get_true_true(result, index, param, mass):
        TT = []
        for i in range(len(result)):
                scan = result[i][0]
                spectrum = index.spectraDict[scan]
                charge = spectrum.ch

                cand = []
                sumint = []
                for j in range(len(result[i][1])):
                        pepObj1 = index.uniquePepObjs[0][result[i][1][j][0][0]]
                        pepObj2 = index.uniquePepObjs[0][result[i][1][j][0][1]]
                        SL = [set(), set()]
                        positions = result[i][1][j][1]

                        for pro in pepObj1.proteinID:
                                cols = pro.split('|R')
                                SL[0].add(cols[1][0])

                        for pro in pepObj2.proteinID:
                                cols = pro.split('|R')
                                SL[1].add(cols[1][0])


                        feature = list(result[i][1][j][2][0])
                        feature.extend(result[i][1][j][2][1])
                        if feature[0] / float(feature[7]) >= 0.20 and feature[1] / float(feature[7]) >= 0.20 and feature[8] / float(feature[15]) >= 0.20 and feature[9] / float(feature[15]) >= 0.20 and feature[2] >= 0.1 and feature[10] >= 0.1 and len(SL[0].intersection(SL[1])) > 0:

                                xl = XLink(pepObj1, pepObj2, positions, charge, mass, param)
                                match = Match(spectrum, xl, mass, param)
                                match.match(mass)

                                cand.append(match)
                                sumint.append(feature[2] + feature[10])

                if len(cand) == 0:
                        continue
                comb = zip(cand, sumint)
                cand = list(zip(*sorted(comb, key = lambda x : x[1], reverse = True))[0])
                sumint = list(zip(*sorted(comb, key = lambda x : x[1], reverse = True))[1])

                TT.append(cand[0])


        for i in range(len(TT)):

                pepObj1 = TT[i].xlink.pepObjs[0]
                pepObj2 = TT[i].xlink.pepObjs[1]

                s = pepObj1.sequence + '\t' + pepObj2.sequence + '\t' + ','.join(pepObj1.proteinID) + '\t' + ','.join(pepObj2.proteinID)
                print s

        return TT

def get_true_false(TrueTrue, param, mass):
        linkermass = param['linkermass']

        ttseq = set()
        m = []
        for match in TrueTrue:
                ttseq.add(match.xlink.pepObjs[0].sequence)
                ttseq.add(match.xlink.pepObjs[1].sequence)
                m.append(match.xlink.mr - linkermass)

        m = max(m)

        fasta = FastaReader(param['uniprot_database']).readFasta()
        patternString = param['patternstring']
        peptides = dict()
        for header, sequence in fasta:
                if 'MOUSE' in header:
                        pepObjInPro = getPepObjsFromProtein(header, sequence, patternString, mass, param)
                        for pep in pepObjInPro:
                                if pep.sequence not in peptides and pep.sequence not in ttseq and pep.pm < m:
                                        peptides[pep.sequence] = pep

        peptides = peptides.values()

        TrueFalse = []
        A = []
        B = []
        for i in range(len(TrueTrue)):
                print i
                match = TrueTrue[i]
                ch = match.spectrum.ch
                tol = match.xlink.mr * 5 * (10 ** (-6))
                A.append([])
                B.append([])

                pepObjs = match.xlink.pepObjs

                for j in range(len(peptides)):
                        if abs(pepObjs[0].pm + peptides[j].pm + linkermass - match.xlink.mr) <= tol:

                                pepseq1 = pepObjs[0].sequence
                                pepseq2 = peptides[j].sequence
                                kPos1 = list(zip(*filter(lambda x : x[1] == param['xresidue'], enumerate(pepseq1[:-1])))[0])
                                kPos2 = list(zip(*filter(lambda x : x[1] == param['xresidue'], enumerate(pepseq2[:-1])))[0])
                                positions = [kPos1[len(kPos1)/2], kPos2[len(kPos2)/2]]
                                xl = XLink(pepObjs[0], peptides[j], positions, ch, mass, param)

                                tf = Match(match.spectrum, xl, mass, param)
                                tf.match(mass)
                                feature = tf.feature
                                if feature[1][0] / float(feature[1][7]) + feature[1][1] / float(feature[1][7]) >= 0.2:
                                        A[-1].append(tf)

                for j in range(len(peptides)):
                        if abs(pepObjs[1].pm + peptides[j].pm + linkermass - match.xlink.mr) <= tol:

                                pepseq1 = pepObjs[1].sequence
                                pepseq2 = peptides[j].sequence
                                kPos1 = list(zip(*filter(lambda x : x[1] == param['xresidue'], enumerate(pepseq1[:-1])))[0])
                                kPos2 = list(zip(*filter(lambda x : x[1] == param['xresidue'], enumerate(pepseq2[:-1])))[0])
                                positions = [kPos1[len(kPos1)/2], kPos2[len(kPos2)/2]]
                                xl = XLink(pepObjs[1], peptides[j], positions, ch, mass, param)

                                tf = Match(match.spectrum, xl, mass, param)
                                tf.match(mass)
                                feature = tf.feature
                                if feature[1][0] / float(feature[1][7]) + feature[1][1] / float(feature[1][7])>= 0.2:
                                        B[-1].append(tf)

        AB = [A, B]
        return AB

def get_false_false(TrueTrue, param, mass):

        linkermass = param['linkermass']

        ttseq = set()
        m = []
        for match in TrueTrue:
                ttseq.add(match.xlink.pepObjs[0].sequence)
                ttseq.add(match.xlink.pepObjs[1].sequence)
                m.append(match.xlink.mr - linkermass)

        mi = int(min(m) - 0.2)
        ma = int(max(m) + 0.2)

        fasta = FastaReader(param['uniprot_database']).readFasta()
        patternString = param['patternstring']
        peptides = dict()
        for header, sequence in fasta:
                if 'YEAST' in header:
                        pepObjInPro = getPepObjsFromProtein(header, sequence, patternString, mass, param)
                        for pep in pepObjInPro:
                                if pep.sequence not in peptides and pep.sequence not in ttseq:
                                        peptides[pep.sequence] = pep

        peptides = peptides.values()

        intdict = dict()
        for peptide in peptides:
                num = int(peptide.pm)

                if num > ma:
                        continue

                if num not in intdict:
                        intdict[num] = [peptide]
                else:
                        intdict[num].append(peptide)

        FalseFalse = []
        for k in range(len(TrueTrue)):
                match = TrueTrue[k]
                print k
                sys.stdout.flush()

                FalseFalse.append([])
                pm = match.xlink.mr - linkermass
                ch = match.spectrum.ch
                tol = match.xlink.mr * 3 * (10 ** (-6))

                mlist = range(500, ma - 500)
                shuffle(mlist)
                mlist = mlist[:25]

                for m in mlist:
                        if m not in intdict:
                                continue
                        shuffle(intdict[m])
                        intdict[m] = intdict[m][:50]

                        for i in range(len(intdict[m])):
                                num = int(pm - intdict[m][i].pm)
                                if num not in intdict:
                                        continue
                                shuffle(intdict[num])
                                intdict[num] = intdict[num][:50]

                                for j in range(len(intdict[num])):
                                        pepseq1 = intdict[m][i].sequence
                                        pepseq2 = intdict[num][j].sequence
                                        kPos1 = list(zip(*filter(lambda x : x[1] == param['xresidue'], enumerate(pepseq1[:-1])))[0])
                                        kPos2 = list(zip(*filter(lambda x : x[1] == param['xresidue'], enumerate(pepseq2[:-1])))[0])
                                        positions = [kPos1[len(kPos1)/2], kPos2[len(kPos2)/2]]
                                        xl = XLink(intdict[m][i], intdict[num][j], positions, ch, mass, param)

                                        if abs(match.xlink.mr - xl.mr) <= tol:
                                                ff = Match(match.spectrum, xl, mass, param)
                                                ff.match(mass)

                                                feature = ff.feature
                                                if feature[0][0] / float(feature[0][7]) + feature[0][1] / float(feature[0][7]) >= 0.15 and feature[1][0] / float(feature[1][7]) + feature[1][1] / float(feature[1][7]) >= 0.15:
                                                        FalseFalse[-1].append(ff)

        return FalseFalse

def calculate_feature(tt, AB, FF):

        # tt
        Xtt = []
        for i in range(len(tt)):
                feature = list(tt[i].feature[0])
                feature.extend(tt[i].feature[1])
                Xtt.append(feature)

        # tf
        A = AB[0]
        B = AB[1]

        AB = []
        for a in A:
                AB.extend(a)
        for b in B:
                AB.extend(b)

        Xtf = []
        for i in range(len(AB)):
                feature = list(AB[i].feature[0])
                feature.extend(AB[i].feature[1])
                Xtf.append(feature)

        # ff
        ff = []
        for f in FF:
                ff.extend(f)

        Xff = []
        for i in range(len(ff)):
                feature = list(ff[i].feature[0])
                feature.extend(ff[i].feature[1])
                Xff.append(feature)

        return [Xtt, Xtf, Xff]

