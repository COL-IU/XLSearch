class Peptide:
	def __init__(self, sequence, proteinID, modification, mass, isNterm):
		self.sequence = sequence
		self.length = len(self.sequence)
		self.proteinID = [proteinID]
		self.modification = modification
		self.massArray = self.getMassArray(mass)
		self.totalResidueMass = self.getTotalResidueMass()
		self.pm = self.getPrecursorMass(mass)
		self.isNterm = isNterm

	def getMassArray(self, mass):
		massArray = []
		sequence = self.sequence
		position = self.modification['position']
		deltaMass = self.modification['deltaMass']

		for i in range(self.length):

			massArray.append(mass[sequence[i].upper()])

		for i in range(len(position)):
			massArray[position[i]] += deltaMass[i]

		return massArray
	def getTotalResidueMass(self):
		totalResidueMass = sum(self.massArray)

		return totalResidueMass
	def getPrecursorMass(self, mass):
		pm = self.totalResidueMass
		pm = pm + mass['Hatom'] * 2 + mass['Oatom']

		return pm
	def getMassList(self, mass, param):
		massList = []
		fwdMass = 0
		massArray = self.massArray
		totalResidueMass = self.totalResidueMass
		useAIon = param['useAIon']

		for i in range(self.length - 1):
			fragment = dict()
	
			fwdMass += massArray[i]
			revMass = totalResidueMass - fwdMass
			
			fragment['b'] = fwdMass + mass['BIonRes']
			fragment['y'] = revMass + mass['YIonRes']
			if useAIon:
				fragment['a'] = fwdMass + mass['AIonRes']
			massList.append(fragment)
		
		return massList	
