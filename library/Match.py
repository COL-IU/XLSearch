from XLink import *
import numpy
import bisect
import copy
import itertools

class Match:
	def __init__(self, spectrum, xlink, mass, param):
		self.spectrum = copy.deepcopy(spectrum)

		self.xlink = xlink
		self.param = param
		self.mr = self.xlink.mr + self.spectrum.ch * mass['Hatom']

                len1 = self.xlink.pepObjs[0].length
                len2 = self.xlink.pepObjs[1].length

                indVectorPref1 = [0] * (len1 - 1)
                indVectorSuff1 = [0] * (len1 - 1)
                indVectorPref2 = [0] * (len2 - 1)
                indVectorSuff2 = [0] * (len2 - 1)
		self.indVector = [indVectorPref1, indVectorSuff1, indVectorPref2, indVectorSuff2]

	def flip(self):
		self.xlink.positions = [self.xlink.positions[1], self.xlink.positions[0]]
		self.xlink.pepObjs = [self.xlink.pepObjs[1], self.xlink.pepObjs[0]]
		self.indVector = [self.indVector[2], self.indVector[3], self.indVector[0], self.indVector[1]]

	def matchIonPerCleaveSite(self, mz, mzList, param):
		ms2tol = self.param['ms2tol']['val']
		if self.param['ms2tol']['measure'] == 'ppm':
			ms2tol = ms2tol * 10**(-6) * mz

		lower = mz - ms2tol
		upper = mz + ms2tol

		li = bisect.bisect_left(mzList, lower)
		ui = bisect.bisect(mzList, upper) - 1	

		return True if li <= ui else False
	def matchIons(self, mz, fragmentIonList):
		param = self.param
		matchedSites = []
		length = len(fragmentIonList)

		for i in range(length):
			mzList = list(zip(*fragmentIonList[i])[0])
			if self.matchIonPerCleaveSite(mz, mzList, param):
				matchedSites.append(i)

		return matchedSites
	def annotate(self, filename, mass):
		f = open(filename, 'w')

		ms2tol = self.param['ms2tol']['val']
		if self.param['ms2tol']['measure'] == 'ppm':
			ms2tol = ms2tol * 10**(-6) * mz

		precprecIonList = self.getPrecPrecIonList(mass)
		(prefList1, suffList1, prefList2, suffList2) = self.xlink.getIonListByCleaveSites(mass, self.param)

		PLIons1 = self.xlink.getPrecursorLinkerIonsPerPeptide(0, mass, self.param).values()
		PLIons1 = list(itertools.chain(*PLIons1))
		PLIons1mzlist = list(zip(*PLIons1)[0])
		PLIons1annlist = list(zip(*PLIons1)[2])

		PLIons2 = self.xlink.getPrecursorLinkerIonsPerPeptide(1, mass, self.param).values()
		PLIons2 = list(itertools.chain(*PLIons2))
		PLIons2mzlist = list(zip(*PLIons2)[0])
		PLIons2annlist = list(zip(*PLIons2)[2])

		mz = self.spectrum.mz
		it = self.spectrum.it

		annotation = ['?'] * len(mz)

		for i in range(len(mz)):
			mi = ms2tol

			for j in range(len(precprecIonList[0])):
				if abs(mz[i] - precprecIonList[0][j]) <= ms2tol:
					if abs(mz[i] - precprecIonList[0][j]) < mi:
						mi = abs(mz[i] - precprecIonList[0][j])
						annotation[i] = precprecIonList[1][j]

		for i in range(len(mz)):
			if annotation[i] != '?':
				continue

			prefSites1 = self.matchIons(mz[i], prefList1)
			suffSites1 = self.matchIons(mz[i], suffList1)
			prefSites2 = self.matchIons(mz[i], prefList2)
			suffSites2 = self.matchIons(mz[i], suffList2)

			mi = ms2tol
			annstr1 = []

			for j in range(len(prefSites1)):
				mzlist = list(zip(*prefList1[prefSites1[j]])[0])
				annlist = list(zip(*prefList1[prefSites1[j]])[2])
			
				for k in range(len(mzlist)):
					if abs(mz[i] - mzlist[k]) <= ms2tol:
						if abs(mz[i] - mzlist[k]) < mi:
							mi = abs(mz[i] - mzlist[k])
							annstr1 = [annlist[k] + ', Alpha']
						elif abs(mz[i] - mzlist[k]) == mi:
							annstr1.append(annlist[k] + ', Alpha')

			for j in range(len(suffSites1)):
				mzlist = list(zip(*suffList1[suffSites1[j]])[0])
				annlist = list(zip(*suffList1[suffSites1[j]])[2])

				for k in range(len(mzlist)):
					if abs(mz[i] - mzlist[k]) <= ms2tol:
						if abs(mz[i] - mzlist[k]) < mi:
							mi = abs(mz[i] - mzlist[k])
							annstr1 = [annlist[k] + ', Alpha']
						elif abs(mz[i] - mzlist[k]) == mi:
							annstr1.append(annlist[k] + ', Alpha')
			mi = ms2tol		
			annstr2 = []
			for j in range(len(prefSites2)):
				mzlist = list(zip(*prefList2[prefSites2[j]])[0])
				annlist = list(zip(*prefList2[prefSites2[j]])[2])

				for k in range(len(mzlist)):
					if abs(mz[i] - mzlist[k]) <= ms2tol:
						if abs(mz[i] - mzlist[k]) < mi:
							mi = abs(mz[i] - mzlist[k])
							annstr2 = [annlist[k] + ', Beta']
						elif abs(mz[i] - mzlist[k]) == mi:
							annstr2.append(annlist[k] + ', Beta')

			for j in range(len(suffSites2)):
				mzlist = list(zip(*suffList2[suffSites2[j]])[0])
				annlist = list(zip(*suffList2[suffSites2[j]])[2])

				for k in range(len(mzlist)):
					if abs(mz[i] - mzlist[k]) <= ms2tol:
						if abs(mz[i] - mzlist[k]) < mi:
							mi = abs(mz[i] - mzlist[k])
							annstr2 = [annlist[k] + ', Beta']
						elif abs(mz[i] - mzlist[k]) == mi:
							annstr2.append(annlist[k] + ', Beta')

			annstr1 = '; '.join(annstr1)
			annstr2 = '; '.join(annstr2)

			if annstr1 != '' and annstr2 != '':
				annotation[i] = ' OR '.join([annstr1, annstr2])
			elif annstr1 != '':
				annotation[i] = annstr1
			elif annstr2 != '':
				annotation[i] = annstr2

			if annstr1 != '' or annstr2 != '':
				continue

			annstr = []
			for j in range(len(PLIons1mzlist)):
				if abs(mz[i] - PLIons1mzlist[j]) <= ms2tol:
					annstr.append(PLIons1annlist[j] + ', Alpha')
			
			for j in range(len(PLIons2mzlist)):
				if abs(mz[i] - PLIons2mzlist[j]) <= ms2tol:
					annstr.append(PLIons2annlist[j] + ', Beta')

			annstr = ' OR '.join(annstr)
			if len(annstr) > 0:
				annotation[i] = annstr

		for i in range(len(mz)):
			f.write('%.4f\t%.1f\t%s\n' % (mz[i], it[i], annotation[i]))

		f.close()
	def match(self, mass):
		ms2tol = self.param['ms2tol']['val']
		if self.param['ms2tol']['measure'] == 'ppm':
			ms2tol = ms2tol * 10**(-6) * mz

		self.filterPrecPrecIons(mass)
		mode = self.param['mode']

		len1 = self.xlink.pepObjs[0].length
		len2 = self.xlink.pepObjs[1].length
		peaks = zip(self.spectrum.mz, self.spectrum.it)
		(prefList1, suffList1, prefList2, suffList2) = self.xlink.getIonListByCleaveSites(mass, self.param)	
		(PLIons1, PLIons2) = self.xlink.getPrecursorLinkerIons(mass, self.param)

		indVectorPref1 = self.indVector[0]
		indVectorSuff1 = self.indVector[1]
		indVectorPref2 = self.indVector[2]
		indVectorSuff2 = self.indVector[3]
		bigInt1 = []
		bigInt2 = []
		PLIons = [0, 0]

		intCutoff = numpy.median(self.spectrum.it)
		doubleMatching = []

		for mz, it in peaks:	
			prefSites1 = self.matchIons(mz, prefList1)
			suffSites1 = self.matchIons(mz, suffList1)
			prefSites2 = self.matchIons(mz, prefList2)
			suffSites2 = self.matchIons(mz, suffList2)

			matched1 = len(prefSites1) != 0 or len(suffSites1) != 0
			matched2 = len(prefSites2) != 0 or len(suffSites2) != 0

			if matched1 and matched2:
				if mode != 'neutral':
					doubleMatching.append((mz, it, prefSites1, suffSites1, prefSites2, suffSites2))
			elif matched1:
				self.assignIndVectorPrefSuff(indVectorPref1, indVectorSuff1, prefSites1, suffSites1)

				if it >= intCutoff:
					bigInt1.append(it)
			elif matched2:
				self.assignIndVectorPrefSuff(indVectorPref2, indVectorSuff2, prefSites2, suffSites2)

				if it >= intCutoff:
					bigInt2.append(it)
			else:
				matchedPLIons1 = any(map(lambda x : abs(x - mz) <= ms2tol, PLIons1))
				matchedPLIons2 = any(map(lambda x : abs(x - mz) <= ms2tol, PLIons2))
				if matchedPLIons1:
					if it >= intCutoff:
						bigInt1.append(it)
					PLIons[0] = 1
				if matchedPLIons2:
					if it >= intCutoff:	
						bigInt2.append(it)
					PLIons[1] = 1	
		
		alphaWins = (sum(indVectorPref1) + sum(indVectorSuff1)) / float(len(indVectorPref1)) > (sum(indVectorPref2) + sum(indVectorSuff2)) / float(len(indVectorPref2))
		betaWins = (sum(indVectorPref1) + sum(indVectorSuff1)) / float(len(indVectorPref1)) < (sum(indVectorPref2) + sum(indVectorSuff2)) / float(len(indVectorPref2))

		if mode == 'conservative':
			if alphaWins:
				for (mz, it, prefSites1, suffSites1, prefSites2, suffSites2) in doubleMatching:
					self.assignIndVectorPrefSuff(indVectorPref1, indVectorSuff1, prefSites1, suffSites1)
					if it >= intCutoff:
						bigInt1.append(it)
			else:
				for (mz, it, prefSites1, suffSites1, prefSites2, suffSites2) in doubleMatching:
					self.assignIndVectorPrefSuff(indVectorPref2, indVectorSuff2, prefSites2, suffSites2)
					if it >= intCutoff:
						bigInt2.append(it)
		elif mode == 'liberal':
			if betaWins:
				for (mz, it, prefSites1, suffSites1, prefSites2, suffSites2) in doubleMatching:
					self.assignIndVectorPrefSuff(indVectorPref1, indVectorSuff1, prefSites1, suffSites1)
					if it >= intCutoff:
						bigInt1.append(it)
			else:
				for (mz, it, prefSites1, suffSites1, prefSites2, suffSites2) in doubleMatching:
					self.assignIndVectorPrefSuff(indVectorPref2, indVectorSuff2, prefSites2, suffSites2)
					if it >= intCutoff:
						bigInt2.append(it)
	
		self.indVector = [indVectorPref1, indVectorSuff1, indVectorPref2, indVectorSuff2] 
		(matchPref1, longestPref1) = self.getMatchesAndLongestSeries(indVectorPref1)
		(matchSuff1, longestSuff1) = self.getMatchesAndLongestSeries(indVectorSuff1)
		(matchPref2, longestPref2) = self.getMatchesAndLongestSeries(indVectorPref2)
		(matchSuff2, longestSuff2) = self.getMatchesAndLongestSeries(indVectorSuff2)	

		intensityList = list(zip(*peaks)[1])
		intensityList = list(zip(*filter(lambda x : x[1] >= intCutoff, enumerate(intensityList)))[1])
		totalIntensity = sum(intensityList)

		intensityFraction1 = sum(bigInt1) / float(totalIntensity)
		intensityFraction2 = sum(bigInt2) / float(totalIntensity)

		countFraction1 = len(bigInt1) / float(len(intensityList))
		countFraction2 = len(bigInt2) / float(len(intensityList))

		matchPref1 = matchPref1 /  float(len1)
		matchSuff1 = matchSuff1 / float(len1)
		longestPref1 = longestPref1 / float(len1)
		longestSuff1 = longestSuff1 / float(len1)

		matchPref2 = matchPref2 / float(len2)
		matchSuff2 = matchSuff2 / float(len2)
		longestPref2 = longestPref2 / float(len2)
		longestSuff2 = longestSuff2 / float(len2)
	
		feature1 = [matchPref1 * float(len1), matchSuff1 * float(len1), intensityFraction1, longestPref1 * float(len1), longestSuff1 * float(len1), countFraction1, PLIons[0], len1]
		feature2 = [matchPref2 * float(len2), matchSuff2 * float(len2), intensityFraction2, longestPref2 * float(len2), longestSuff2 * float(len2), countFraction2, PLIons[1], len2]

		self.feature = (feature1, feature2)
	def filterPrecPrecIons(self, mass):
		precprecIonList = self.getPrecPrecIonList(mass)[0]
		ms2tol = self.param['ms2tol']['val']
		if self.param['ms2tol']['measure'] == 'ppm':
			ms2tol = ms2tol * 10**(-6) * mz


		peaks = zip(self.spectrum.mz, self.spectrum.it)

		MZ = []
		IT = []
		for mz, it in peaks:
			if all(map(lambda x : abs(x - mz) > ms2tol, precprecIonList)):
				MZ.append(mz)
				IT.append(it)
		self.spectrum.mz = MZ
		self.spectrum.it = IT

		return len(MZ)
	def getPrecPrecIonList(self, mass):
		mr = self.mr
		ch = self.spectrum.ch
		h2oLoss = self.param['neutralloss']['h2oLoss']['mass'];
		nh3Loss = self.param['neutralloss']['nh3Loss']['mass'];
		neutronmass = mass['neutronmass']

		precprecIonList = [mr, mr + h2oLoss, mr + nh3Loss, mr + h2oLoss + nh3Loss]
		precprecIonList *= 2
		for i in range(4, len(precprecIonList)):
			precprecIonList[i] += neutronmass
		precprecIonList = map(lambda x : x / ch, precprecIonList)

		return (precprecIonList, ['A+L+B', 'A+L+B-H2O', 'A+L+B-NH3', 'A+L+B-H2O-NH3', 'A+L+Bi', 'A+L+B-H2Oi', 'A+L+B-NH3i', 'A+L+B-H2O-NH3i'])
	def assignIndVectorPrefSuff(self, indVectorPref, indVectorSuff, prefSites, suffSites):
		for s in prefSites:
			indVectorPref[s] = 1
		for s in suffSites:
			indVectorSuff[s] = 1
	def getMatchesAndLongestSeries(self, indVector):
		vec = [0]
		vec.extend(indVector)
		vec.append(0)

		zeroIndex = list(zip(*filter(lambda x : x[1] == 0, enumerate(vec)))[0])
		matches = len(indVector) - len(zeroIndex) + 2 
		consecutiveSize = 0

		for i in range(len(zeroIndex) - 1):
			if zeroIndex[i + 1] - zeroIndex[i] - 1>= consecutiveSize:
				consecutiveSize = zeroIndex[i + 1] - zeroIndex[i] - 1 

		return (matches, consecutiveSize)
	def getMatchInfo(self, index):
		pepDict = index.uniquePepObjs[1]
		pepIndex = (pepDict[self.xlink.pepObjs[0].sequence], pepDict[self.xlink.pepObjs[1].sequence])
		positions = self.xlink.positions
		feature = self.feature

		return (pepIndex, positions, feature)
	def getInfoString(self):
		return self.xlink.getInfoString() + '_' + self.spectrum.title
