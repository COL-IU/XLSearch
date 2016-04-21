from xlink import XLink
import numpy
import bisect
import copy
import itertools

class Match:
	def __init__(self, spec, xlink, mass, param):
		self.spec = copy.deepcopy(spec)

		self.xlink = xlink
		self.param = param
		self.mol_weight = self.xlink.mol_weight + self.spec.ch * mass['Hatom']

		len1 = self.xlink.pep[0].length
		len2 = self.xlink.pep[1].length

		ind_vector_pref1 = [0] * (len1 - 1)
		ind_vector_suff1 = [0] * (len1 - 1)
		ind_vector_pref2 = [0] * (len2 - 1)
		ind_vector_suff2 = [0] * (len2 - 1)
		self.ind_vector = [ind_vector_pref1, ind_vector_suff1, ind_vector_pref2, ind_vector_suff2]

	def flip(self):
		self.xlink.pos = [self.xlink.pos[1], self.xlink.pos[0]]
		self.xlink.pep = [self.xlink.pep[1], self.xlink.pep[0]]
		self.ind_vector = [self.ind_vector[2], self.ind_vector[3], self.ind_vector[0], self.ind_vector[1]]

	def match_ion_per_cleave_site(self, mz, mz_list, param):
		ms2_tol = self.param['ms2_tol']['val']
		if self.param['ms2_tol']['measure'] == 'ppm':
			ms2_tol = ms2_tol * 10**(-6) * mz

		lower = mz - ms2_tol
		upper = mz + ms2_tol

		li = bisect.bisect_left(mz_list, lower)
		ui = bisect.bisect(mz_list, upper) - 1	

		return True if li <= ui else False
	def match_ions(self, mz, frag_ion_list):
		param = self.param
		matched_sites = []
		length = len(frag_ion_list)

		for i in range(length):
			mz_list = list(zip(*frag_ion_list[i])[0])
			if self.match_ion_per_cleave_site(mz, mz_list, param):
				matched_sites.append(i)

		return matched_sites
	def annotate(self, filename, mass):
		f = open(filename, 'w')

		precprec_ion_list = self.get_precprec_ion_list(mass)
		(pref_list1, suff_list1, pref_list2, suff_list2) = self.xlink.get_ion_list_by_cleave_sites(mass, self.param)

		pl_ions1 = self.xlink.get_prec_linker_ions_per_pep(0, mass, self.param).values()
		pl_ions1 = list(itertools.chain(*pl_ions1))
		pl_ions1_mz_list = list(zip(*pl_ions1)[0])
		pl_ions1_ann_list = list(zip(*pl_ions1)[2])

		pl_ions2 = self.xlink.get_prec_linker_ions_per_pep(1, mass, self.param).values()
		pl_ions2 = list(itertools.chain(*pl_ions2))
		pl_ions2_mz_list = list(zip(*pl_ions2)[0])
		pl_ions2_ann_list = list(zip(*pl_ions2)[2])

		mz = self.spec.mz
		it = self.spec.it

		ann = ['?'] * len(mz)

		for i in range(len(mz)):
			ms2_tol = self.param['ms2_tol']['val']
			if self.param['ms2_tol']['measure'] == 'ppm':
				ms2_tol = ms2_tol * 10**(-6) * mz[i]
			mi = ms2_tol

			for j in range(len(precprec_ion_list[0])):
				if abs(mz[i] - precprec_ion_list[0][j]) <= ms2_tol:
					if abs(mz[i] - precprec_ion_list[0][j]) < mi:
						mi = abs(mz[i] - precprec_ion_list[0][j])
						ann[i] = precprec_ion_list[1][j]

		for i in range(len(mz)):
			if ann[i] != '?':
				continue

			pref_sites1 = self.match_ions(mz[i], pref_list1)
			suff_sites1 = self.match_ions(mz[i], suff_list1)
			pref_sites2 = self.match_ions(mz[i], pref_list2)
			suff_sites2 = self.match_ions(mz[i], suff_list2)

			ms2_tol = self.param['ms2_tol']['val']
			if self.param['ms2_tol']['measure'] == 'ppm':
				ms2_tol = ms2_tol * 10**(-6) * mz[i]
			mi = ms2_tol
			ann_str1 = []

			for j in range(len(pref_sites1)):
				mz_list = list(zip(*pref_list1[pref_sites1[j]])[0])
				ann_list = list(zip(*pref_list1[pref_sites1[j]])[2])
			
				for k in range(len(mz_list)):
					if abs(mz[i] - mz_list[k]) <= ms2_tol:
						if abs(mz[i] - mz_list[k]) < mi:
							mi = abs(mz[i] - mz_list[k])
							ann_str1 = [ann_list[k] + ', Alpha']
						elif abs(mz[i] - mz_list[k]) == mi:
							ann_str1.append(ann_list[k] + ', Alpha')

			for j in range(len(suff_sites1)):
				mz_list = list(zip(*suff_list1[suff_sites1[j]])[0])
				ann_list = list(zip(*suff_list1[suff_sites1[j]])[2])

				for k in range(len(mz_list)):
					if abs(mz[i] - mz_list[k]) <= ms2_tol:
						if abs(mz[i] - mz_list[k]) < mi:
							mi = abs(mz[i] - mz_list[k])
							ann_str1 = [ann_list[k] + ', Alpha']
						elif abs(mz[i] - mz_list[k]) == mi:
							ann_str1.append(ann_list[k] + ', Alpha')
			mi = ms2_tol
			ann_str2 = []
			for j in range(len(pref_sites2)):
				mz_list = list(zip(*pref_list2[pref_sites2[j]])[0])
				ann_list = list(zip(*pref_list2[pref_sites2[j]])[2])

				for k in range(len(mz_list)):
					if abs(mz[i] - mz_list[k]) <= ms2_tol:
						if abs(mz[i] - mz_list[k]) < mi:
							mi = abs(mz[i] - mz_list[k])
							ann_str2 = [ann_list[k] + ', Beta']
						elif abs(mz[i] - mz_list[k]) == mi:
							ann_str2.append(ann_list[k] + ', Beta')

			for j in range(len(suff_sites2)):
				mz_list = list(zip(*suff_list2[suff_sites2[j]])[0])
				ann_list = list(zip(*suff_list2[suff_sites2[j]])[2])

				for k in range(len(mz_list)):
					if abs(mz[i] - mz_list[k]) <= ms2_tol:
						if abs(mz[i] - mz_list[k]) < mi:
							mi = abs(mz[i] - mz_list[k])
							ann_str2 = [ann_list[k] + ', Beta']
						elif abs(mz[i] - mz_list[k]) == mi:
							ann_str2.append(ann_list[k] + ', Beta')

			ann_str1 = '; '.join(ann_str1)
			ann_str2 = '; '.join(ann_str2)

			if ann_str1 != '' and ann_str2 != '':
				ann[i] = ' OR '.join([ann_str1, ann_str2])
			elif ann_str1 != '':
				ann[i] = ann_str1
			elif ann_str2 != '':
				ann[i] = ann_str2

			if ann_str1 != '' or ann_str2 != '':
				continue

			ann_str = []
			for j in range(len(pl_ions1_mz_list)):
				if abs(mz[i] - pl_ions1_mz_list[j]) <= ms2_tol:
					ann_str.append(pl_ions1_ann_list[j] + ', Alpha')
			
			for j in range(len(pl_ions2_mz_list)):
				if abs(mz[i] - pl_ions2_mz_list[j]) <= ms2_tol:
					ann_str.append(pl_ions2_ann_list[j] + ', Beta')

			ann_str = ' OR '.join(ann_str)
			if len(ann_str) > 0:
				ann[i] = ann_str

		for i in range(len(mz)):
			f.write('%.4f\t%.1f\t%s\n' % (mz[i], it[i], ann[i]))

		f.close()
	def match(self, mass):
		self.filter_precprec_ions(mass)
		mode = self.param['mode']

		len1 = self.xlink.pep[0].length
		len2 = self.xlink.pep[1].length
		peaks = zip(self.spec.mz, self.spec.it)
		(pref_list1, suff_list1, pref_list2, suff_list2) = self.xlink.get_ion_list_by_cleave_sites(mass, self.param)
		(pl_ions1, pl_ions2) = self.xlink.get_prec_linker_ions(mass, self.param)

		ind_vector_pref1 = self.ind_vector[0]
		ind_vector_suff1 = self.ind_vector[1]
		ind_vector_pref2 = self.ind_vector[2]
		ind_vector_suff2 = self.ind_vector[3]
		big_int1 = []
		big_int2 = []
		pl_ions = [0, 0]

		int_cutoff = numpy.median(self.spec.it)
		double_matching = []

		for mz, it in peaks:	
			ms2_tol = self.param['ms2_tol']['val']
			if self.param['ms2_tol']['measure'] == 'ppm':
				ms2_tol = ms2_tol * 10**(-6) * mz			

			pref_sites1 = self.match_ions(mz, pref_list1)
			suff_sites1 = self.match_ions(mz, suff_list1)
			pref_sites2 = self.match_ions(mz, pref_list2)
			suff_sites2 = self.match_ions(mz, suff_list2)

			matched1 = len(pref_sites1) != 0 or len(suff_sites1) != 0
			matched2 = len(pref_sites2) != 0 or len(suff_sites2) != 0

			if matched1 and matched2:
				if mode != 'neutral':
					double_matching.append((mz, it, pref_sites1, suff_sites1, pref_sites2, suff_sites2))
			elif matched1:
				self.assign_ind_vector_prefsuff(ind_vector_pref1, ind_vector_suff1, pref_sites1, suff_sites1)

				if it >= int_cutoff:
					big_int1.append(it)
			elif matched2:
				self.assign_ind_vector_prefsuff(ind_vector_pref2, ind_vector_suff2, pref_sites2, suff_sites2)

				if it >= int_cutoff:
					big_int2.append(it)
			else:
				matched_pl_ions1 = any(map(lambda x : abs(x - mz) <= ms2_tol, pl_ions1))
				matched_pl_ions2 = any(map(lambda x : abs(x - mz) <= ms2_tol, pl_ions2))
				if matched_pl_ions1:
					if it >= int_cutoff:
						big_int1.append(it)
					pl_ions[0] = 1
				if matched_pl_ions2:
					if it >= int_cutoff:	
						big_int2.append(it)
					pl_ions[1] = 1	
		
		alpha_wins = (sum(ind_vector_pref1) + sum(ind_vector_suff1)) / float(len(ind_vector_pref1)) > (sum(ind_vector_pref2) + sum(ind_vector_suff2)) / float(len(ind_vector_pref2))
		beta_wins = (sum(ind_vector_pref1) + sum(ind_vector_suff1)) / float(len(ind_vector_pref1)) < (sum(ind_vector_pref2) + sum(ind_vector_suff2)) / float(len(ind_vector_pref2))

		if mode == 'conservative':
			if alpha_wins:
				for (mz, it, pref_sites1, suff_sites1, pref_sites2, suff_sites2) in double_matching:
					self.assign_ind_vector_prefsuff(ind_vector_pref1, ind_vector_suff1, pref_sites1, suff_sites1)
					if it >= int_cutoff:
						big_int1.append(it)
			else:
				for (mz, it, pref_sites1, suff_sites1, pref_sites2, suff_sites2) in double_matching:
					self.assign_ind_vector_prefsuff(ind_vector_pref2, ind_vector_suff2, pref_sites2, suff_sites2)
					if it >= int_cutoff:
						big_int2.append(it)
		elif mode == 'liberal':
			if beta_wins:
				for (mz, it, pref_sites1, suff_sites1, pref_sites2, suff_sites2) in double_matching:
					self.assign_ind_vector_prefsuff(ind_vector_pref1, ind_vector_suff1, pref_sites1, suff_sites1)
					if it >= int_cutoff:
						big_int1.append(it)
			else:
				for (mz, it, pref_sites1, suff_sites1, pref_sites2, suff_sites2) in double_matching:
					self.assign_ind_vector_prefsuff(ind_vector_pref2, ind_vector_suff2, pref_sites2, suff_sites2)
					if it >= int_cutoff:
						big_int2.append(it)
	
		self.ind_vector = [ind_vector_pref1, ind_vector_suff1, ind_vector_pref2, ind_vector_suff2] 
		(match_pref1, longest_pref1) = self.get_matches_and_longest_series(ind_vector_pref1)
		(match_suff1, longest_suff1) = self.get_matches_and_longest_series(ind_vector_suff1)
		(match_pref2, longest_pref2) = self.get_matches_and_longest_series(ind_vector_pref2)
		(match_suff2, longest_suff2) = self.get_matches_and_longest_series(ind_vector_suff2)	

		int_list = list(zip(*peaks)[1])
		int_list = list(zip(*filter(lambda x : x[1] >= int_cutoff, enumerate(int_list)))[1])
		total_int = sum(int_list)

		int_frac1 = sum(big_int1) / float(total_int)
		int_frac2 = sum(big_int2) / float(total_int)

		count_frac1 = len(big_int1) / float(len(int_list))
		count_frac2 = len(big_int2) / float(len(int_list))

		match_pref1 = match_pref1 /  float(len1)
		match_suff1 = match_suff1 / float(len1)
		longest_pref1 = longest_pref1 / float(len1)
		longest_suff1 = longest_suff1 / float(len1)

		match_pref2 = match_pref2 / float(len2)
		match_suff2 = match_suff2 / float(len2)
		longest_pref2 = longest_pref2 / float(len2)
		longest_suff2 = longest_suff2 / float(len2)
	
		feature1 = [match_pref1 * float(len1), match_suff1 * float(len1), int_frac1, longest_pref1 * float(len1), longest_suff1 * float(len1), count_frac1, pl_ions[0], len1]
		feature2 = [match_pref2 * float(len2), match_suff2 * float(len2), int_frac2, longest_pref2 * float(len2), longest_suff2 * float(len2), count_frac2, pl_ions[1], len2]

		self.feature = (feature1, feature2)
	def filter_precprec_ions(self, mass):
		precprec_ion_list = self.get_precprec_ion_list(mass)[0]
		peaks = zip(self.spec.mz, self.spec.it)

		MZ = []
		IT = []
		for mz, it in peaks:
			ms2_tol = self.param['ms2_tol']['val']
			if self.param['ms2_tol']['measure'] == 'ppm':
				ms2_tol = ms2_tol * 10**(-6) * mz

			if all(map(lambda x : abs(x - mz) > ms2_tol, precprec_ion_list)):
				MZ.append(mz)
				IT.append(it)
		self.spec.mz = MZ
		self.spec.it = IT

		return len(MZ)
	def get_precprec_ion_list(self, mass):
		mol_weight = self.mol_weight
		ch = self.spec.ch
		h2o_loss = self.param['neutral_loss']['h2o_loss']['mass'];
		nh3_loss = self.param['neutral_loss']['nh3_loss']['mass'];
		neutron_mass = mass['neutron_mass']

		precprec_ion_list = [mol_weight, mol_weight + h2o_loss, mol_weight + nh3_loss, mol_weight + h2o_loss + nh3_loss]
		precprec_ion_list *= 2
		for i in range(4, len(precprec_ion_list)):
			precprec_ion_list[i] += neutron_mass
		precprec_ion_list = map(lambda x : x / ch, precprec_ion_list)

		return (precprec_ion_list, ['A+L+B', 'A+L+B-H2O', 'A+L+B-NH3', 'A+L+B-H2O-NH3', 'A+L+Bi', 'A+L+B-H2Oi', 'A+L+B-NH3i', 'A+L+B-H2O-NH3i'])
	def assign_ind_vector_prefsuff(self, ind_vector_pref, ind_vector_suff, pref_sites, suff_sites):
		for s in pref_sites:
			ind_vector_pref[s] = 1
		for s in suff_sites:
			ind_vector_suff[s] = 1
	def get_matches_and_longest_series(self, ind_vector):
		vec = [0]
		vec.extend(ind_vector)
		vec.append(0)

		zero_index = list(zip(*filter(lambda x : x[1] == 0, enumerate(vec)))[0])
		matches = len(ind_vector) - len(zero_index) + 2 
		consecutive_size = 0

		for i in range(len(zero_index) - 1):
			if zero_index[i + 1] - zero_index[i] - 1 >= consecutive_size:
				consecutive_size = zero_index[i + 1] - zero_index[i] - 1

		return (matches, consecutive_size)
	def get_match_info(self, index):
		pep_dict = index.unique_pep[1]
		pep_index = (pep_dict[self.xlink.pep[0].seq], pep_dict[self.xlink.pep[1].seq])
		pos = self.xlink.pos
		feature = self.feature

		return [pep_index, pos, feature]
	def get_info_string(self):
		return self.xlink.get_info_string() + '_' + self.spec.title
