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
		self.prec_mass_pep_index_tuple = self.get_prec_mass_pep_index_tuple(mass)
		self.search_index = self.build_index()

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
	def get_prec_mass_pep_index_tuple(self, mass):
		unique_pep = self.unique_pep[0]
	
		linker_mass = self.param['linker_mass']

		mass = []
		index1 = []
		index2 = []

		for i in range(len(unique_pep) - 1):
			for j in range(i + 1, len(unique_pep)):
				prec_mass1 = unique_pep[i].prec_mass
				prec_mass2 = unique_pep[j].prec_mass
				
				prec_mass_xlink = prec_mass1 + prec_mass2 + linker_mass
				mass.append(prec_mass_xlink)
				index1.append(i)
				index2.append(j)
			
		prec_mass_pep_index_tuple = sorted(zip(mass, index1, index2), key = lambda tup : tup[0])
		mass = zip(*prec_mass_pep_index_tuple)[0] 
		index1 = zip(*prec_mass_pep_index_tuple)[1]
		index2 = zip(*prec_mass_pep_index_tuple)[2]
		prec_mass_pep_index_tuple = (mass, index1, index2)

		return prec_mass_pep_index_tuple
	def build_index(self):
		titles = self.spec_dict.keys()
		search_index = dict()

		for title in titles:
			tup = self.find_candidates(title)
			search_index[tup[0]] = tup[1]

		return search_index

	def find_candidates(self, title):
		spec_dict = self.spec_dict
		prec_mass_pep_index_tuple = self.prec_mass_pep_index_tuple

		spec = spec_dict[title]
		mass = prec_mass_pep_index_tuple[0]

		ms1_tol = self.param['ms1_tol']['val']
		if self.param['ms1_tol']['measure'] == 'ppm':
			ms1_tol = ms1_tol * 10**(-6) * spec.mol_weight
				
		upper = spec.mol_weight + ms1_tol
		lower = spec.mol_weight - ms1_tol

		li = bisect.bisect_left(mass, lower)
		ui = bisect.bisect(mass, upper) - 1

		index_candidates = range(li, ui + 1)

		return (title, index_candidates)
