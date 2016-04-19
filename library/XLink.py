from Peptide import *
import itertools

class XLink:
	def __init__(self, pepObj1, pepObj2, positions, charge, mass, param):
		self.pepObjs = [pepObj1, pepObj2]
		self.positions = positions
		self.charge = charge
		self.linkerMass = param['linkermass']
		self.mr = self.pepObjs[0].pm + self.pepObjs[1].pm + self.linkerMass

	def adjustMassListPerPeptide(self, indexThis, indexThat, mass, param):
		massList = self.pepObjs[indexThis].getMassList(mass, param)
		linkPos = self.positions[indexThis]
		pm = self.pepObjs[indexThat].pm
		linkerMass = self.linkerMass
		length = self.pepObjs[indexThis].length

		useAIon = param['useAIon']

		for i in range(length - 1):
			if i < linkPos:
				massList[i]['y'] += (linkerMass + pm)
			else:
				massList[i]['b'] += (linkerMass + pm)
				if useAIon:
					massList[i]['a'] += (linkerMass + pm)

		return massList
	def getFragmentIons(self, mass, param):
		fragmentIonList = []
		fragmentIonList.append(self.getFragmentIonsPerPeptide(0, 1, mass, param))
		fragmentIonList.append(self.getFragmentIonsPerPeptide(1, 0, mass, param))

		return fragmentIonList
	def getFragmentIonsPerPeptide(self, indexThis, indexThat, mass, param):
		massList = self.adjustMassListPerPeptide(indexThis, indexThat, mass, param)

		thisPep = self.pepObjs[indexThis].sequence
		thatPep = set(self.pepObjs[indexThat].sequence)

		linkPos = self.positions[indexThis]
		length = self.pepObjs[indexThis].length

		chPreXL = param['chargePreXlinkIons']
		chPostXL = param['chargePostXlinkIons']
		chPreXL = range(chPreXL[0], int(min(chPreXL[1], self.charge)) + 1)
		chPostXL = range(chPostXL[0], int(min(chPostXL[1], self.charge)) + 1)
		h2oLoss = param['neutralloss']['h2oLoss']
		nh3Loss = param['neutralloss']['nh3Loss']
		h2oGain = param['neutralloss']['h2oGain']
		protonmass = mass['Hatom']	
		useAIon = param['useAIon']

		singleCharge = []
		fragmentIonList = []

		for i in range(length - 1):
			fragment = dict()
			fragment['b'] = dict()
			fragment['b']['none'] = massList[i]['b']
			fragment['y'] = dict()
			fragment['y']['none'] = massList[i]['y']

			if useAIon:
				fragment['a'] = dict()
				fragment['a']['none'] = massList[i]['a']
							
			prefAA = set(thisPep[0:(i + 1)])
			suffAA = set(thisPep[(i + 1):])

			if i < linkPos:
				suffAA.update(thatPep)
			else:
				prefAA.update(thatPep)

			if h2oLoss['aa'].intersection(prefAA):
				fragment['b']['-H2O'] = fragment['b']['none'] + h2oLoss['mass']
				if useAIon:
					fragment['a']['-H2O'] = fragment['a']['none'] + h2oLoss['mass']

			if h2oLoss['aa'].intersection(suffAA):
				fragment['y']['-H2O'] = fragment['y']['none'] + h2oLoss['mass']

			if nh3Loss['aa'].intersection(prefAA):
				fragment['b']['-NH3'] = fragment['b']['none'] + nh3Loss['mass']
				if useAIon:
					fragment['a']['-NH3'] = fragment['a']['none'] + nh3Loss['mass']

			if nh3Loss['aa'].intersection(suffAA):
				fragment['y']['-NH3'] = fragment['y']['none'] + nh3Loss['mass']

			if len(h2oGain['aa']) != 0:
				if i == length - 2:
					fragment['b']['+H2O'] = fragment['b']['none'] + h2oGain['mass']
					if useAIon:
						fragment['a']['+H2O'] = fragment['a']['none'] + h2oGain['mass']

			singleCharge.append(fragment)

		for i in range(length - 1):
			key1 = singleCharge[i].keys()
			fragment = dict()

			for j in key1:
				chList = []
				if i < linkPos:
					if j == 'b' or j == 'a':
						chList = chPreXL
					else:
						chList = chPostXL
				else:
					if j == 'b' or j == 'a':
						chList = chPostXL
					else:
						chList = chPreXL

				key2 = singleCharge[i][j].keys()
				fragment[j] = dict()
				for k in key2:
					if j == 'b' or j == 'a':
						ionString = j + ', ' + str(i + 1) + ', ' + k  
					else:
						ionString = j + ', ' + str(length - i - 1) + ', ' + k
					fragment[j][k] = self.getMZListFromMass(singleCharge[i][j][k], chList, True, ionString, mass)
		 		
			fragmentIonList.append(fragment)

		return fragmentIonList
	def getMZListFromMass(self, mz, chList, singlyCharged, ionString, mass):
		mzList = []
		protonmass = mass['Hatom']

		for i in range(len(chList)):
			if singlyCharged:
				imz = (mz + (chList[i] - 1) * protonmass) / chList[i]
			else:
				imz = (mz + chList[i] * protonmass) / chList[i]
			string = str(ionString + ', ' + str(chList[i]) + '+')
			peak = (imz, chList[i], string)
			mzList.append(peak)

		return mzList
	def getIonListByCleaveSites(self, mass, param):
		fragmentIonList = self.getFragmentIons(mass, param)		
		(prefList1, suffList1) = self.getIonListByCleaveSitesPerPeptide(fragmentIonList[0])
		(prefList2, suffList2) = self.getIonListByCleaveSitesPerPeptide(fragmentIonList[1])
		ionList = (prefList1, suffList1, prefList2, suffList2)

		return ionList
	def getIonListByCleaveSitesPerPeptide(self, fragmentIonList):
		ionList = []
		for cleaveSiteIons in fragmentIonList:
			prefList = []
			suffList = []

			for ionType, ions in cleaveSiteIons.items():
				if ionType == 'a' or ionType == 'b':	
					for variants, ionsch in ions.items():
						prefList.extend(ionsch)
				else:
					for variants, ionsch in ions.items():
						suffList.extend(ionsch)
			prefList = sorted(prefList, key = lambda tup : tup[0])
			suffList = sorted(suffList, key = lambda tup : tup[0])

			ionList.append((prefList, suffList))
		ionList = zip(*ionList)
		prefList = list(ionList[0])
		suffList = [i for i in reversed(list(ionList[1]))]

		return (prefList, suffList)
	def getPrecursorLinkerIons(self, mass, param):
		PLIons1 = self.getPrecursorLinkerIonsPerPeptide(0, mass, param).values()
		PLIons1 = list(itertools.chain(*PLIons1))
		PLIons1 = list(zip(*PLIons1)[0])
		PLIons2 = self.getPrecursorLinkerIonsPerPeptide(1, mass, param).values()
		PLIons2 = list(itertools.chain(*PLIons2))
		PLIons2 = list(zip(*PLIons2)[0])

		return (PLIons1, PLIons2)
	def getPrecursorLinkerIonsPerPeptide(self, indexThis, mass, param):
		pm = self.pepObjs[indexThis].pm
		pepstr = 'Alpha' if indexThis == 0 else 'Beta'
		linkerMass = self.linkerMass
		protonmass = mass['Hatom']
		h2oLoss = param['neutralloss']['h2oLoss']['mass'];
		nh3Loss = param['neutralloss']['nh3Loss']['mass'];
		chList = range(1, int(self.charge) + 1)

		precursorLinkerIons = dict()

		precursorLinkerIons[pepstr] = pm
		precursorLinkerIons[pepstr + '-H2O'] = pm + h2oLoss
		precursorLinkerIons[pepstr + '-NH3'] = pm + nh3Loss
		precursorLinkerIons[pepstr + '+L'] = pm + linkerMass
		precursorLinkerIons[pepstr + '+L-H2O'] = pm + linkerMass + h2oLoss
		precursorLinkerIons[pepstr + '+L-NH3'] = pm + linkerMass + nh3Loss
		precursorLinkerIons[pepstr + '+L-2H2O'] = pm + linkerMass + 2 * h2oLoss

		for ionString, mz in precursorLinkerIons.items():
			precursorLinkerIons[ionString] = self.getMZListFromMass(precursorLinkerIons[ionString], chList, False, ionString, mass)

		return precursorLinkerIons
	def getInfoString(self):
		string = ''
		if self.pepObjs[0].sequence < self.pepObjs[1].sequence:
			string += self.pepObjs[0].sequence + '_' + self.pepObjs[1].sequence + '_' + str(self.positions[0] + 1) + '_' + str(self.positions[1] + 1)
		else:
			string += self.pepObjs[1].sequence + '_' + self.pepObjs[0].sequence + '_' + str(self.positions[1] + 1) + '_' + str(self.positions[0] + 1)
		string += '_' + str(self.charge)

		return string
