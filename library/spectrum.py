import copy
class Spectrum:
	def __init__(self, title, scan_num, prec_mz, ch, mz, it, ret_time, mass):
		self.title = title
		self.scan_num = scan_num
		self.prec_mz = prec_mz
		self.ch = ch
		self.mz = mz
		self.it = it
		self.size = len(self.mz)
		self.mol_weight = self.prec_mz * self.ch - mass['Hatom'] * self.ch
		self.ret_time = ret_time
	def deisotope(self, mass, max_iso, tol):
		deisotoped = copy.deepcopy(self)
		mz = self.mz
		it = self.it
		isotope_inc = mass['isotope_inc']		
		MZ = []
		IT = []		
		alignment = ['          '] * len(mz)
		i = 0

		while i < len(mz) - 1:
			alignment[i] = '%.6f' % mz[i]
			count = 0
			for j in isotope_inc:
				reference_mz = map(lambda x: x * j + mz[i], range(1, max_iso + 1))
				observed_mz = mz[i + 1 : i + max_iso + 1]
				for (mz1, mz2) in zip(reference_mz, observed_mz):	
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
