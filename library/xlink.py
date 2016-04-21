from peptide import Peptide
import itertools
import sys

class XLink:
	def __init__(self, pep1, pep2, pos, ch, mass, param):
		self.pep = [pep1, pep2]
		self.pos = pos
		self.ch = ch
		self.linker_mass = param['linker_mass']
		self.mol_weight = self.pep[0].prec_mass + self.pep[1].prec_mass + self.linker_mass

	def adjust_mass_list_per_pep(self, index_this, index_that, mass, param):
		mass_list = self.pep[index_this].get_mass_list(mass, param)
		link_pos = self.pos[index_this]
		prec_mass = self.pep[index_that].prec_mass
		linker_mass = self.linker_mass
		length = self.pep[index_this].length
		use_a_ion = param['use_a_ion']

		for i in range(length - 1):
			if i < link_pos:
				mass_list[i]['y'] += (linker_mass + prec_mass)
			else:
				mass_list[i]['b'] += (linker_mass + prec_mass)
				if use_a_ion:
					mass_list[i]['a'] += (linker_mass + prec_mass)

		return mass_list
	def get_frag_ions(self, mass, param):
		frag_ion_list = []
		frag_ion_list.append(self.get_frag_ions_per_pep(0, 1, mass, param))
		frag_ion_list.append(self.get_frag_ions_per_pep(1, 0, mass, param))
		return frag_ion_list

	def get_frag_ions_per_pep(self, index_this, index_that, mass, param):
		mass_list = self.adjust_mass_list_per_pep(index_this, index_that, mass, param)
		this_pep = self.pep[index_this].seq
		that_pep = set(self.pep[index_that].seq)
		link_pos = self.pos[index_this]
		length = self.pep[index_this].length
		ch_pre_xl = param['ch_pre_xlink_ions']
		ch_post_xl = param['ch_post_xlink_ions']
		ch_pre_xl = range(ch_pre_xl[0], int(min(ch_pre_xl[1], self.ch)) + 1)
		if len(ch_pre_xl) == 0:
			print 'ch_pre_xl: is empty!'
			print 'execution will terminate!'
			sys.exit(1)
		ch_post_xl = range(ch_post_xl[0], int(min(ch_post_xl[1], self.ch)) + 1)
		if len(ch_post_xl) == 0:
			print 'ch_post_xl: is empty!'
			print 'execution will terminate!'
			sys.exit(1)
		h2o_loss = param['neutral_loss']['h2o_loss']
		nh3_loss = param['neutral_loss']['nh3_loss']
		h2o_gain = param['neutral_loss']['h2o_gain']
		proton_mass = mass['Hatom']	
		use_a_ion = param['use_a_ion']
		single_ch = []
		frag_ion_list = []

		for i in range(length - 1):
			frag = dict()
			frag['b'] = dict()
			frag['b']['none'] = mass_list[i]['b']
			frag['y'] = dict()
			frag['y']['none'] = mass_list[i]['y']

			if use_a_ion:
				frag['a'] = dict()
				frag['a']['none'] = mass_list[i]['a']
							
			pref_AA = set(this_pep[0:(i + 1)])
			suff_AA = set(this_pep[(i + 1):])

			if i < link_pos:
				suff_AA.update(that_pep)
			else:
				pref_AA.update(that_pep)

			if h2o_loss['aa'].intersection(pref_AA):
				frag['b']['-H2O'] = frag['b']['none'] + h2o_loss['mass']
				if use_a_ion:
					frag['a']['-H2O'] = frag['a']['none'] + h2o_loss['mass']

			if h2o_loss['aa'].intersection(suff_AA):
				frag['y']['-H2O'] = frag['y']['none'] + h2o_loss['mass']

			if nh3_loss['aa'].intersection(pref_AA):
				frag['b']['-NH3'] = frag['b']['none'] + nh3_loss['mass']
				if use_a_ion:
					frag['a']['-NH3'] = frag['a']['none'] + nh3_loss['mass']

			if nh3_loss['aa'].intersection(suff_AA):
				frag['y']['-NH3'] = frag['y']['none'] + nh3_loss['mass']

			if len(h2o_gain['aa']) != 0:
				if i == length - 2:
					frag['b']['+H2O'] = frag['b']['none'] + h2o_gain['mass']
					if use_a_ion:
						frag['a']['+H2O'] = frag['a']['none'] + h2o_gain['mass']

			single_ch.append(frag)

		for i in range(length - 1):
			key1 = single_ch[i].keys()
			frag = dict()

			for j in key1:
				ch_list = []
				if i < link_pos:
					if j == 'b' or j == 'a':
						ch_list = ch_pre_xl
					else:
						ch_list = ch_post_xl
				else:
					if j == 'b' or j == 'a':
						ch_list = ch_post_xl
					else:
						ch_list = ch_pre_xl

				key2 = single_ch[i][j].keys()
				frag[j] = dict()
				for k in key2:
					if j == 'b' or j == 'a':
						ion_string = j + ', ' + str(i + 1) + ', ' + k  
					else:
						ion_string = j + ', ' + str(length - i - 1) + ', ' + k
					frag[j][k] = self.get_mz_list_from_mass(single_ch[i][j][k], ch_list, True, ion_string, mass)
		 		
			frag_ion_list.append(frag)

		return frag_ion_list

	def get_mz_list_from_mass(self, mz, ch_list, singly_charged, ion_string, mass):
		mz_list = []
		proton_mass = mass['Hatom']

		for i in range(len(ch_list)):
			if singly_charged:
				imz = (mz + (ch_list[i] - 1) * proton_mass) / ch_list[i]
			else:
				imz = (mz + ch_list[i] * proton_mass) / ch_list[i]
			string = str(ion_string + ', ' + str(ch_list[i]) + '+')
			peak = (imz, ch_list[i], string)
			mz_list.append(peak)

		return mz_list

	def get_ion_list_by_cleave_sites(self, mass, param):
		frag_ion_list = self.get_frag_ions(mass, param)
		(pref_list1, suff_list1) = self.get_ion_list_by_cleave_sites_per_pep(frag_ion_list[0])
		(pref_list2, suff_list2) = self.get_ion_list_by_cleave_sites_per_pep(frag_ion_list[1])
		ion_list = (pref_list1, suff_list1, pref_list2, suff_list2)

		return ion_list

	def get_ion_list_by_cleave_sites_per_pep(self, frag_ion_list):
		ion_list = []

		for cleave_site_ions in frag_ion_list:
			pref_list = []
			suff_list = []

			for ion_type, ions in cleave_site_ions.items():
				if ion_type == 'a' or ion_type == 'b':	
					for loss_type, ions_ch in ions.items():
						pref_list.extend(ions_ch)
				else:
					for loss_type, ions_ch in ions.items():
						suff_list.extend(ions_ch)
			pref_list = sorted(pref_list, key = lambda tup : tup[0])
			suff_list = sorted(suff_list, key = lambda tup : tup[0])
			ion_list.append((pref_list, suff_list))

		ion_list = zip(*ion_list)
		pref_list = list(ion_list[0])
		suff_list = [i for i in reversed(list(ion_list[1]))]

		return (pref_list, suff_list)
	def get_prec_linker_ions(self, mass, param):
		pl_ions1 = self.get_prec_linker_ions_per_pep(0, mass, param).values()
		pl_ions1 = list(itertools.chain(*pl_ions1))
		pl_ions1 = list(zip(*pl_ions1)[0])
		pl_ions2 = self.get_prec_linker_ions_per_pep(1, mass, param).values()
		pl_ions2 = list(itertools.chain(*pl_ions2))
		pl_ions2 = list(zip(*pl_ions2)[0])

		return (pl_ions1, pl_ions2)

	def get_prec_linker_ions_per_pep(self, index_this, mass, param):
		prec_mass = self.pep[index_this].prec_mass
		pep_str = 'Alpha' if index_this == 0 else 'Beta'
		linker_mass = self.linker_mass
		proton_mass = mass['Hatom']
		h2o_loss = param['neutral_loss']['h2o_loss']['mass'];
		nh3_loss = param['neutral_loss']['nh3_loss']['mass'];
		ch_list = range(1, int(self.ch) + 1)

		prec_linker_ions = dict()

		prec_linker_ions[pep_str] = prec_mass
		prec_linker_ions[pep_str + '-H2O'] = prec_mass + h2o_loss
		prec_linker_ions[pep_str + '-NH3'] = prec_mass + nh3_loss
		prec_linker_ions[pep_str + '+L'] = prec_mass + linker_mass
		prec_linker_ions[pep_str + '+L-H2O'] = prec_mass + linker_mass + h2o_loss
		prec_linker_ions[pep_str + '+L-NH3'] = prec_mass + linker_mass + nh3_loss
		prec_linker_ions[pep_str + '+L-2H2O'] = prec_mass + linker_mass + 2 * h2o_loss

		for ion_string, mz in prec_linker_ions.items():
			prec_linker_ions[ion_string] = self.get_mz_list_from_mass(prec_linker_ions[ion_string], ch_list, False, ion_string, mass)

		return prec_linker_ions

	def get_info_string(self):
		string = ''
		if self.pep[0].seq < self.pep[1].seq:
			string += self.pep[0].seq + '_' + self.pep[1].seq + '_' + str(self.pos[0] + 1) + '_' + str(self.pos[1] + 1)
		else:
			string += self.pep[1].seq + '_' + self.pep[0].seq + '_' + str(self.pos[1] + 1) + '_' + str(self.pos[0] + 1)
		string += '_' + str(self.ch)

		return string
