import copy
class Spectrum:
	def __init__(self, title, scanNum, precursorMZ, ch, mz, it, rt, mass):
		self.title = title
		self.scanNum = scanNum
		self.precursorMZ = precursorMZ
		self.ch = ch
		self.mz = mz
		self.it = it
		self.size = len(self.mz)
		self.mr = self.precursorMZ * self.ch - mass['Hatom'] * self.ch
		self.rt = rt
	def deisotope(self, mass, maxIso, tol):
		deisotoped = copy.deepcopy(self)
		mz = self.mz
		it = self.it
		isotopeInc = mass['isotopeInc']		

		MZ = []
		IT = []		
		alignment = ['          '] * len(mz)
		i = 0
		while i < len(mz) - 1:
			alignment[i] = '%.6f' % mz[i]
			count = 0
			for j in isotopeInc:
				referenceMZ = map(lambda x: x * j + mz[i], range(1, maxIso + 1))
				observedMZ = mz[i + 1 : i + maxIso + 1]
				for (mz1, mz2) in zip(referenceMZ, observedMZ):	
					if abs(mz1 - mz2) > tol:
						break	
					count += 1
				if count > 0:
					break

			MZ.append(mz[i])
			IT.append(max(it[i : i + count + 1]))
			i += (count + 1)
		
		if count == 0:
			MZ.append(mz[-1])
			IT.append(it[-1])
			alignment[i] = '%.6f' % mz[-1]

		deisotoped.mz = MZ
		deisotoped.it = IT
		deisotoped.size = len(deisotoped.it)

		mz = map(lambda x : '%.6f' % x, mz)
		it = map(lambda x : '%.6f' % x, it)

		return (deisotoped, zip(alignment, mz, it))
