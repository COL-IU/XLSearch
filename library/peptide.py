class Peptide:
	def __init__(self, seq, pro_id, modif, mass, is_nterm):
		self.seq = seq
		self.length = len(self.seq)
		self.pro_id = [pro_id]
		self.modif = modif
		self.mass_array = self.get_mass_array(mass)
		self.total_res_mass = self.get_total_res_mass()
		self.prec_mass = self.get_prec_mass(mass)
		self.is_nterm = is_nterm

	def get_mass_array(self, mass):
		mass_array = []
		seq = self.seq
		pos = self.modif['position']
		delta_mass = self.modif['delta_mass']

		for i in range(self.length):
			mass_array.append(mass[seq[i].upper()])

		for i in range(len(pos)):
			mass_array[pos[i]] += delta_mass[i]

		return mass_array

	def get_total_res_mass(self):
		total_res_mass = sum(self.mass_array)

		return total_res_mass

	def get_prec_mass(self, mass):
		prec_mass = self.total_res_mass
		prec_mass = prec_mass + mass['Hatom'] * 2 + mass['Oatom']

		return prec_mass

	def get_mass_list(self, mass, param):
		mass_list = []
		fwd_mass = 0
		mass_array = self.mass_array
		total_res_mass = self.total_res_mass
		use_a_ion = param['use_a_ion']

		for i in range(self.length - 1):
			frag = dict()	
			fwd_mass += mass_array[i]
			rev_mass = total_res_mass - fwd_mass
			frag['b'] = fwd_mass + mass['b_ion_res']
			frag['y'] = rev_mass + mass['y_ion_res']
			if use_a_ion:
				frag['a'] = fwd_mass + mass['a_ion_res']
			mass_list.append(frag)

		return mass_list
