from fastareader import FastaReader
from peptide import Peptide
from xlink import XLink
from utility import *

import bisect
import re

class EnumIndexBuilder:
	def __init__(self, fasta_filename, spec_dict, mass, param):
		self.fasta_filename = fasta_filename
		self.param = param

		self.spec_dict = spec_dict
		self.unique_pep = self.get_unique_pep(mass)
		self.mass_bin = self.get_mass_bin()
#		self.search_index = self.build_index()

	def get_unique_pep(self, mass):
		fasta_filename = self.fasta_filename
		param = self.param
		pattern_string = param['pattern_string']
		fr = FastaReader(fasta_filename)
		fasta = fr.read_fasta()

		pep_dict = dict()

		for header, seq in fasta:
			pep_in_pro = get_pep_from_pro(header, seq, pattern_string, mass, param)
			for pep in pep_in_pro:
				if pep.seq not in pep_dict:
					pep_dict[pep.seq] = pep
				else:
					pep_dict[pep.seq].pro_id.extend(pep.pro_id)
					pep_dict[pep.seq].pro_id = list(set(pep_dict[pep.seq].pro_id))
					pep_dict[pep.seq].pro_id.sort()
		unique_pep = []
		pep_keys = pep_dict.keys()

		for i in range(len(pep_keys)):
			unique_pep.append(pep_dict[pep_keys[i]])
			pep_dict[pep_keys[i]] = i
		
		return (unique_pep, pep_dict)
	def build_index(self):
		titles = self.spec_dict.keys()
		search_index = dict()

		for title in titles:
			tup = self.get_candidates(title)
			search_index[tup[0]] = tup[1]

		return search_index

	def get_mass_bin(self):
		unique_pep = self.unique_pep[0]

		mass_bin = dict()
		for i in range(len(unique_pep)):
			int_mass = round(unique_pep[i].prec_mass) 
			if int_mass not in mass_bin:
				mass_bin[int_mass] = [i]
			else:
				mass_bin[int_mass].append(i)

		return mass_bin

	def get_candidates(self, title):
		unique_pep = self.unique_pep[0]
		spec_dict = self.spec_dict
		spec = spec_dict[title]
		mass_bin = self.mass_bin

		ms1_tol = self.param['ms1_tol']['val']
		if self.param['ms1_tol']['measure'] == 'ppm':
			ms1_tol = ms1_tol * 10**(-6) * spec.mol_weight

		upper = spec.mol_weight + ms1_tol
		lower = spec.mol_weight - ms1_tol

		total_mass = spec.mol_weight
		linker_mass = self.param['linker_mass']

		mass_range = mass_bin.keys()
		index_candidates = set()

		mass = []
		index1 = []
		index2 = []

		for i in range(len(mass_range)):
			peps_index1 = mass_bin[mass_range[i]]
		
			int_mass_c = round((total_mass - linker_mass) - mass_range[i])
			peps_index2 = []
			if int_mass_c in mass_bin:
				peps_index2.extend(mass_bin[int_mass_c])
			if int_mass_c - 1 in mass_bin:
				peps_index2.extend(mass_bin[int_mass_c - 1])
			if int_mass_c + 1 in mass_bin:
				peps_index2.extend(mass_bin[int_mass_c + 1])
			if len(peps_index2) == 0:
				continue

			for i1 in range(len(peps_index1)):
				for i2 in range(len(peps_index2)):
					prec_mass1 = unique_pep[peps_index1[i1]].prec_mass
					prec_mass2 = unique_pep[peps_index2[i2]].prec_mass

					obs_mass = prec_mass1 + prec_mass2 + self.param['linker_mass']
					if obs_mass <= upper and obs_mass >= lower:
						if peps_index1[i1] < peps_index2[i2]:
							index_candidates.add(str(peps_index1[i1]) + '_' + str(peps_index2[i2]))
						elif peps_index2[i2] < peps_index1[i1]:
							index_candidates.add(str(peps_index2[i2]) + '_' + str(peps_index1[i1]))

#					if obs_mass <= upper and obs_mass >= lower and peps_index1[i1] < peps_index2[i2]:
#						mass.append(obs_mass)
#						index1.append(peps_index1[i1])
#						index2.append(peps_index2[i2])

#		if len(mass) > 0:
#			prec_mass_pep_index_tuple = sorted(zip(mass, index1, index2), key = lambda tup : tup[0])
#			mass = zip(*prec_mass_pep_index_tuple)[0]
#			index1 = zip(*prec_mass_pep_index_tuple)[1]
#			index2 = zip(*prec_mass_pep_index_tuple)[2]

#			li = bisect.bisect_left(mass, lower)
#			ui = bisect.bisect(mass, upper) - 1
#			mass_index = range(li, ui + 1)
			
#			for j in range(len(mass_index)):				
#				index_candidates.append([index1[mass_index[j]], index2[mass_index[j]]])
#				if index1[mass_index[j]] < index2[mass_index[j]]:
#					index_candidates.add(str(index1[mass_index[j]]) + '_' + str(index2[mass_index[j]]))
#				elif index2[mass_index[j]] < index1[mass_index[j]]:
#					index_candidates.add(str(index2[mass_index[j]]) + '_' + str(index1[mass_index[j]]))

#					if obs_mass >= lower and obs_mass <= upper:
#						if peps_index1[i1] < peps_index2[i2]:
#							index_candidates.add(str(peps_index1[i1]) + '_' + str(peps_index2[i2]))
#						elif peps_index1[i1] > peps_index2[i2]:
#							index_candidates.add(str(peps_index2[i2]) + '_' + str(peps_index1[i1]))

		index_candidates = list(index_candidates)
		for i in range(len(index_candidates)):
			index_candidates[i] = index_candidates[i].split('_')
			index_candidates[i][0] = int(index_candidates[i][0])
			index_candidates[i][1] = int(index_candidates[i][1])
		return index_candidates
	#	return [title, index_candidates]
