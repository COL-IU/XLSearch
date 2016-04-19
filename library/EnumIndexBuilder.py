from FastaReader import *
from Peptide import *
from XLink import *
from Utility import *
#from multiprocessing.dummy import Pool as ThreadPool
import bisect
import re

class EnumIndexBuilder:
	def __init__(self, fastaFileName, spectraDict, mass, param):
		self.fastaFileName = fastaFileName
		self.param = param

		self.spectraDict = spectraDict
		self.uniquePepObjs = self.getUniquePepObjs(mass)
		self.precMassPepIndexTuple = self.getPrecMassPepIndexTuple(mass)
		self.searchIndex = self.buildIndex()

	def getUniquePepObjs(self, mass):
		fastaFileName = self.fastaFileName
		param = self.param
		patternString = param['patternstring']
		fr = FastaReader(fastaFileName)
		fasta = fr.readFasta()

		pepDict = dict()

		for header, sequence in fasta:
			pepObjInPro = getPepObjsFromProtein(header, sequence, patternString, mass, param)
			for pepObj in pepObjInPro:
				if pepObj.sequence not in pepDict:
					pepDict[pepObj.sequence] = pepObj
				else:
					pepDict[pepObj.sequence].proteinID.extend(pepObj.proteinID)
					pepDict[pepObj.sequence].proteinID.sort()

		uniquePepObjs = []
		pepKeys = pepDict.keys()

		for i in range(len(pepKeys)):
			uniquePepObjs.append(pepDict[pepKeys[i]])
			pepDict[pepKeys[i]] = i
		
		return (uniquePepObjs, pepDict)
	def getPrecMassPepIndexTuple(self, mass):
		uniquePepObjs = self.uniquePepObjs[0]
	
		linkerMass = self.param['linkermass']

		mass = []
		index1 = []
		index2 = []

		for i in range(len(uniquePepObjs) - 1):
			for j in range(i + 1, len(uniquePepObjs)):
				pm1 = uniquePepObjs[i].pm
				pm2 = uniquePepObjs[j].pm
				
				pmxlink = pm1 + pm2 + linkerMass
				mass.append(pmxlink)
				index1.append(i)
				index2.append(j)
			
		precMassPepIndexTuple = sorted(zip(mass, index1, index2), key = lambda tup : tup[0])
		mass = zip(*precMassPepIndexTuple)[0] 
		index1 = zip(*precMassPepIndexTuple)[1]
		index2 = zip(*precMassPepIndexTuple)[2]
		precMassPepIndexTuple = (mass, index1, index2)

		return precMassPepIndexTuple
	def buildIndex(self):
		titles = self.spectraDict.keys()
		searchIndex = dict()

#		pool = ThreadPool(4)
#		results = pool.map(self.findCandidates, titles)				
#		pool.close()
#		pool.join()

#		keys = zip(*results)[0]
#		values = zip(*results)[1]
#		for i in range(len(keys)):
#			searchIndex[keys[i]] = values[i]
		for title in titles:
			tup = self.findCandidates(title)
			searchIndex[tup[0]] = tup[1]

		return searchIndex

	def findCandidates(self, title):
		spectraDict = self.spectraDict
		precMassPepIndexTuple = self.precMassPepIndexTuple

		spectrum = spectraDict[title]
		mass = precMassPepIndexTuple[0]
#		ms1tol = self.param['ms1tol']
#		upper = spectrum.mr * (1 + ms1tol)
#		lower = spectrum.mr * (1 - ms1tol)

		ms1tol = self.param['ms1tol']['val']
		if self.param['ms1tol']['measure'] == 'ppm':
			ms1tol = ms1tol * 10**(-6) * spectrum.mr
				
		upper = spectrum.mr + ms1tol
		lower = spectrum.mr - ms1tol

		li = bisect.bisect_left(mass, lower)
		ui = bisect.bisect(mass, upper) - 1

		indexCandidates = range(li, ui + 1)

		return (title, indexCandidates)
