import re
import math
import glob
import pickle
import sys
import numpy as np

from peptide import Peptide
from xlink import XLink
from match import Match
from mzxmlreader import MZXMLReader
from fastareader import FastaReader
from random import shuffle

def logit(coef, x):
	v = x
	x = [1]
	x.extend(v)
	val = 0

	for i in range(len(coef)):
		val += coef[i] * x[i]
	
	return float(1) / (1 + math.exp(-val))

def adjust_prior_marginal(prior_tr, model_output, marginal, iterations, rel_change_perc):
	posterior_tr = []
	for per_spec in model_output:
		posterior_tr.extend(per_spec)

	n = len(posterior_tr)
	prior = prior_tr
	posterior = [0.0] * n

	for i in range(iterations):
		for k in range(n):
			denominator = (float(prior) / prior_tr) * posterior_tr[k] + (float(1 - prior) / (1 - prior_tr)) * (1 - posterior_tr[k])
			if prior_tr * denominator == 0:
				print 'denominator equals 0!'
				print 'execution will terminate!'
				sys.exit(1)
			posterior[k] = float(prior * posterior_tr[k]) / (prior_tr * denominator)
		prior_prev = prior
		prior = 0.0
		for j in range(n):
			prior = prior + posterior[j] * marginal[j]
		prior = prior / sum(marginal)
#		print 'iteration = %d, %.10f' % (i, prior)
		if i > 0 and abs(float(prior_prev - prior) / prior_prev) <= rel_change_perc:
#			print '%f' % abs(float(prior_prev - prior) / prior_prev)
			break

	new_model_output = []
	begin = 0

	for i in range(len(model_output)):
		end = begin + len(model_output[i])
		new_model_output.append(posterior[begin : end])
		begin = end

	return new_model_output

def get_marginal(p21, p11, p12, p22):
	alpha_T = []
	beta_T = []
	alpha_F =  []
	beta_F = []

	for i in range(len(p21)):
		for j in range(len(p21[i])):
			denominator = 1.0 - (p11[i][j] - p12[i][j]) * (p21[i][j] - p22[i][j])

			if denominator == 0:
				print 'denominator equals 0!'
				print 'execution will terminate!'
				sys.exit(1)

			at = (p12[i][j] + p22[i][j] * (p11[i][j] - p12[i][j])) / denominator
			bt = (p22[i][j] + p12[i][j] * (p21[i][j] - p22[i][j])) / denominator

			af = 1.0 - at
			bf = 1.0 - bt
			alpha_T.append(at)
			beta_T.append(bt)
			alpha_F.append(af)
			beta_F.append(bf)

	return [alpha_T, beta_T, alpha_F, beta_F]

def get_matches_per_spec(mass, param, index, title):
	spec_dict = index.spec_dict
	unique_pep = index.unique_pep[0]
#	search_index = index.search_index
	x_residue = param['x_residue']

	index_list = index.get_candidates(title)
	spec = spec_dict[title]

	matches = []
	for i in range(len(index_list)):
		index1 = index_list[i][0]
		index2 = index_list[i][1]
		
		pep1 = unique_pep[index1]
		pep2 = unique_pep[index2]

		pep_sorted = sorted([pep1, pep2], key = lambda x : x.seq)
		pep1 = pep_sorted[0]
		pep2 = pep_sorted[1]

		ch = spec_dict[title].ch

		mz = spec_dict[title].mz
		it = spec_dict[title].it

		k_pos1 = []
		k_pos2 = []
		if param['ntermxlink'] == True:
			if pep1.is_nterm == True:
				k_pos1.append(0)
			if pep2.is_nterm == True:
				k_pos2.append(0)

		pep_seq1 = pep1.seq
		k_pos1.extend(list(zip(*filter(lambda x : x[1] == x_residue, enumerate(pep_seq1[:-1])))[0]))
		pep_seq2 = pep2.seq
		k_pos2.extend(list(zip(*filter(lambda x : x[1] == x_residue, enumerate(pep_seq2[:-1])))[0]))

		for p1 in k_pos1:
			for p2 in k_pos2:
				pos = [p1, p2]
				xl = XLink(pep1, pep2, pos, ch, mass, param)
			
				match = Match(spec, xl, mass, param)
				match.match(mass)
				matches.append(match.get_match_info(index))

	return matches

def get_pep_from_pro(header, pro_seq, pattern_string, mass, param):
	missed_sites = param['missed_sites']
	min_len = param['min_length']
	max_len = param['max_length']
	mod_res = param['mod_res']
	pattern = re.compile(pattern_string)
	sites = [0]

	for i in range(len(pro_seq)):
		if i == len(pro_seq) - 1:
			sites.append(i + 1)
		elif (pro_seq[i] == 'K' or pro_seq[i] == 'R') and pro_seq[i + 1] != 'P':
			sites.append(i + 1)

	pep_seqs = []

	for i in range(len(sites)):
		if i < len(sites) - missed_sites - 1:
			for j in range(missed_sites + 1):
				pep_seq = pro_seq[sites[i] : sites[i + j + 1]]
				if len(pep_seq) >= min_len and len(pep_seq) <= max_len and pattern.match(pep_seq):
					pep_seqs.append(pep_seq)
		else:
			for j in range(i + 1, len(sites)):
				pep_seq = pro_seq[sites[i] : sites[j]]
				if len(pep_seq) >= min_len and len(pep_seq) <= max_len and pattern.match(pep_seq):
					pep_seqs.append(pep_seq)

	pep_seqs = list(set(pep_seqs))
	peps = []

	for pep_seq in pep_seqs:
		modif = dict(position=[], delta_mass=[])
		is_nterm = True if pep_seq == pro_seq[:len(pep_seq)] else False

		peps.append(Peptide(pep_seq, header, modif, mass, is_nterm))
		if len(mod_res) != 0:
			mod_mass = param['mod_mass']
			index = [i for i, aa in enumerate(pep_seq) if aa == mod_res]
			if len(index) != 0:
				max_var_mod = min(param['max_var_mod'], len(index))

				for i in range(1, max_var_mod + 1):
					mod_cfg = subset(index, i)
					for j in range(len(mod_cfg)):
						pep_seq_mod = list(pep_seq)	
						modif = dict(position=[], delta_mass=[])
						for k in range(len(mod_cfg[j])):
							pep_seq_mod[mod_cfg[j][k]] = pep_seq_mod[mod_cfg[j][k]].lower()
							modif['position'].append(mod_cfg[j][k])
							modif['delta_mass'].append(mod_mass)

						pep_seq_mod = ''.join(pep_seq_mod)
						peps.append(Peptide(pep_seq_mod, header, modif, mass, is_nterm))

	return peps

def read_param(filename):
	param = dict(
		use_a_ion=True,
		verbose=False,
		ch_pre_xlink_ions=[1, 3],
		ch_post_xlink_ions=[2, 5],
		base_peak_int = 100.0,
		dynamic_range = 0.001,
		missed_sites = 2,
		min_length = 4,
		max_length = 51,
		mod_res = '',
		mod_mass = 0.0,
		linker_mass = 136.10005,
		ms1_tol = dict(measure='ppm', val=5),
		ms2_tol = dict(measure='da', val=0.01),
		min_mz = 200,
		max_mz = 2000,
		mode = 'conservative',
		x_residue = 'K',
		aa = 'ACDEFGHIKLMNPQRSTVWY',
		rel_change_perc = 0.01,
		neutral_loss=dict(
			h2o_loss=dict(
				mass=-18.010565,
				aa=set('ACDEFGHIKLMNPQRSTVWY')),
			nh3_loss=dict(
				mass=-17.026549,
				aa=set('ACDEFGHIKLMNPQRSTVWY')),
			h2o_gain=dict(
				mass=18.010565,
				aa=set('ACDEFGHIKLMNPQRSTVWY'))),
		model_TT_TF = [0.0] * 17,
		model_TF_FF = [0.0] * 17,
		nTT = 169,
		nTF = 8568,
		nFF = 91242)
		
	mass = dict(
		A=71.037114,
		R=156.101111,
		N=114.042927,
		D=115.026943,
		C=103.009184,
		E=129.042593,
		Q=128.058578,
		G=57.021464,
		H=137.058912,
		I=113.084064,
		L=113.084064,
		K=128.094963,
		M=131.040485,
		F=147.068414,
		P=97.052764,
		S=87.032028,
		T=101.047678,
		W=186.079313,
		Y=163.063329,
		V=99.068414,
		Hatom=1.007825032,
		Oatom=15.99491462,
		neutron_mass = 1.008701,
		b_ion_res=1.0078246,
		a_ion_res=-26.9870904,
		y_ion_res=19.0183888,
		isotope_inc = [1.008701/4, 1.008701/3, 1.008701/2, 1.008701/1])

	f = open(filename)	
	lines = f.readlines()
	for l in lines:
		l = l[:-1]
		cols= l.split('\t')
		if len(l) == 0 or l[0] == '#' or len(cols) < 2:
			continue

		name = cols[0]
		val = cols[1]
		if name == 'database':
			param['database'] = val
		elif name == 'MS_data_directory':
			param['ms_data'] = val
			if val[-1] != '/':
				param['ms_data'] += '/'
		elif name == 'XLresidue':
			param['x_residue'] = val
		elif name == 'ms1tol_unit':
			param['ms1_tol']['measure'] = val
		elif name == 'ms1tol_val':
			param['ms1_tol']['val'] = float(val)	########## used to be int(val)
		elif name == 'ms2tol_unit':
			param['ms2_tol']['measure'] = val
		elif name == 'ms2tol_val':
			param['ms2_tol']['val'] = float(val)
		elif name == 'linker_mass':
			param['linker_mass'] = float(val)
		elif name == 'miss_cleave':
			param['missed_sites'] = int(val)
		elif name == 'include_a_ions':
			param['use_a_ion'] = True if val.lower() == 'true' else False 
		elif name == 'min_peplen':
			param['min_length'] = int(val)
		elif name == 'max_peplen':
			param['max_length'] = int(val)
		elif name == 'fix_mod_res':
			param['fix_mod_res'] = val
		elif name == 'fix_mod_mass':
			param['fix_mod_mass'] = float(val)
		elif name == 'var_mod_res':
			param['mod_res'] = val
		elif name == 'var_mod_mass':
			param['mod_mass'] = float(val)
		elif name == 'min_preXL_ions_ch':
			param['ch_pre_xlink_ions'][0] = int(val)
		elif name == 'max_preXL_ions_ch':
			param['ch_pre_xlink_ions'][1] = int(val)
		elif name == 'min_postXL_ions_ch':
			param['ch_post_xlink_ions'][0] = int(val)
		elif name == 'max_postXL_ions_ch':
			param['ch_post_xlink_ions'][1] = int(val)
		elif name == 'target_database':
			param['target_database'] = val
		elif name =='uniprot_database':
			param['uniprot_database'] = val
		elif name == 'max_iterations':
			param['max_iterations'] = int(val)
		elif name == 'annotate_spec':
			param['annotation'] = True if val.lower() == 'true' else False
		elif name == 'deisotope':
			param['deisotope'] = True if val.lower() == 'true' else False
		elif name == 'ndeisotope':
			param['ndeisotope'] = int(val)
		elif name == 'ntermxlink':
			param['ntermxlink'] = True if val.lower() == 'true' else False
		elif name == 'decoy_string':
			param['decoy_string'] = val
		elif name == 'cutoff':
			param['cutoff'] = float(val)
		elif name == 'is_unique':
			param['is_unique'] = val
		elif name == 'true_true_psm_file':
			param['true_true_psm_file'] = val
		elif name == 'max_var_mod':
			param['max_var_mod'] = int(val)
		elif len(name) >= 4 and name[:2] == 'CI':
			if len(name) == 4:
				s = int(name[2:])
				param['model_TT_TF'][s] = float(val)
			elif len(name) == 5:
				s = int(name[3:])
				param['model_TF_FF'][s] = float(val)
		elif name == 'nTT':
			param['nTT'] = int(val)
		elif name == 'nTF':
			param['nTF'] = int(val)
		elif name == 'nFF':
			param['nFF'] = int(val)

	param['pattern_string'] = '^[' + param['aa'] + ']*' + param['x_residue'] + '[' + param['aa'] + ']+$'
	param['prior_tr_TT_TF'] = float(param['nTT']) / (param['nTT'] + param['nTF'])
	param['prior_tr_TF_FF'] = float(param['nTF']) / (param['nTF'] + param['nFF'])

	f.close()
	if len(param['fix_mod_res']) > 0:
		mass[param['fix_mod_res']] += param['fix_mod_mass']

	return [param, mass]

def read_spec(directory, param, mass):
	files = glob.glob(directory + '*.mzXML')

	spec_dict = dict()
	total = []
	for filename in files:
		reader = MZXMLReader(filename)
		spec = reader.get_spec_list(mass, param)
		total.append(spec)

	ss = []
	for i in range(len(total)):
		ss.append(set())

		tmp = []
		for j in range(len(total[i])):
			if total[i][j].ret_time >= 0 and total[i][j].ret_time <= 110*60 and total[i][j].ch >= 2 and total[i][j].ch <= 7:
				tmp.append(total[i][j])

		tmp = sorted(tmp, key = lambda s : s.mol_weight)

		tolerance = 0.01
		lower_ratio = 0.3
		upper_ratio = 1 / float(lower_ratio)

		for j in range(len(tmp) - 1):
			MZ = []
			IT = []
			mz = tmp[j].mz
			it = tmp[j].it
			last_index = 0
			ik = 0
			jk = 0
			for ik in range(len(mz)):
				if last_index == 0:
					jk = 0
				else:
					jk = last_index
				while not (last_index > len(mz) - 1 or jk > len(mz) - 1 or ik > len(mz) - 1 or mz[jk] > mz[ik] + tolerance):
					if mz[jk] <= mz[ik] - tolerance:
						last_index = jk

					ratio = float(it[ik]) / float(it[jk])
					if abs(mz[ik] - mz[jk]) <= tolerance and ratio >= lower_ratio and ratio <= upper_ratio:
						MZ.append(mz[ik])
						IT.append(it[ik])
					jk = jk + 1	
			if len(MZ) >= 25:
				spec_dict[tmp[j].title] = tmp[j]
	return spec_dict

def get_tophits(index, result):
	model_TT_TF = index.param['model_TT_TF']
	model_TF_FF = index.param['model_TF_FF']
	prior_tr_TT_TF = index.param['prior_tr_TT_TF']
	prior_tr_TF_FF = index.param['prior_tr_TF_FF']
	max_iterations = index.param['max_iterations']
	rel_change_perc = index.param['rel_change_perc']

	p21 = []
	p11 = []
	p12 = []
	p22 = []

	for i in range(len(result)):
		print i

		p21.append([])
		p11.append([])
		p12.append([])
		p22.append([])

		for j in range(len(result[i][1])):

			feature = result[i][1][j][2]
			x = list(feature[0])
			x.extend(feature[1])
			x_flip = list(feature[1])
			x_flip.extend(feature[0])
			b = model_TT_TF

			p21[-1].append(logit(b, x))
			p11[-1].append(logit(b, x_flip))

			b = model_TF_FF

			p12[-1].append(logit(b, x))
			p22[-1].append(logit(b, x_flip))

	[alpha_T, beta_T, alpha_F, beta_F] = get_marginal(p21, p11, p12, p22)
	p21 = adjust_prior_marginal(prior_tr_TT_TF, p21, alpha_T, max_iterations, rel_change_perc)
	p11 = adjust_prior_marginal(prior_tr_TT_TF, p11, beta_T, max_iterations, rel_change_perc)
	p12 = adjust_prior_marginal(prior_tr_TF_FF, p12, beta_F, max_iterations, rel_change_perc)
	p22 = adjust_prior_marginal(prior_tr_TF_FF, p22, alpha_F, max_iterations, rel_change_perc)

	for i in range(len(result)):
		print i
#		result[i] = list(result[i])

		for j in range(len(result[i][1])):

			pep1 = index.unique_pep[0][result[i][1][j][0][0]]
			pep2 = index.unique_pep[0][result[i][1][j][0][1]]

			ap21 = p21[i][j]
			ap11 = p11[i][j]
			ap12 = p12[i][j]
			ap22 = p22[i][j]

			denominator = 1 - (ap11 - ap12) * (ap21 - ap22)
			if denominator == 0:
				print 'denominator equals 0!'
				print 'execution will terminate!'
				sys.exit(1)

			marginal_alaph_T = (ap12 + ap22 * (ap11 - ap12)) / denominator
			marginal_beta_T = (ap22 + ap12 * (ap21 - ap22)) / denominator

			prob1 = ap11 * marginal_beta_T
			prob2 = ap21 * marginal_alaph_T
			score = (prob1 + prob2) / float(2)

			info = {'alpha' : marginal_alaph_T, 'beta' : marginal_beta_T, 'prob' : [prob1, prob2], 'score' : score}

#			result[i][1][j] = list(result[i][1][j])
			result[i][1][j].append(info)

	for r in result:
		r[1] = sorted(r[1], key = lambda x : x[3]['score'], reverse = True)

	result = sorted(result, key = lambda x : x[1][0][3]['score'], reverse = True)

	tophits = []

	for r in result:
		scan = r[0]
		pep = [index.unique_pep[0][r[1][0][0][0]].seq, index.unique_pep[0][r[1][0][0][1]].seq]
		pos = [int(r[1][0][1][0]), int(r[1][0][1][1])]
		pro = [index.unique_pep[0][r[1][0][0][0]].pro_id, index.unique_pep[0][r[1][0][0][1]].pro_id]
		ch = int(scan.split('.')[-1])
		score = r[1][0][3]['score']
		alpha = r[1][0][3]['alpha']
		beta = r[1][0][3]['beta']
		tophits.append([pep, pos, pro, ch, score, alpha, beta, scan])

	return tophits

def write_results(output_file, tophits):
	f = open(output_file, 'w')
	f.write('Rank\tPep_alpha\tPep_beta\tSite_alpha\tSite_beta\tPro_alpha\tPro_beta\tCharge\tpr(alpha=T,beta=T)\tpr(alpha=T)\tpr(beta=T)\tSpectrum\n')
	for i in range(len(tophits)):
		f.write('%d\t' % (i + 1))
		f.write('%s\t%s\t' % (tophits[i][0][0], tophits[i][0][1]))
		f.write('%d\t%d\t' % (tophits[i][1][0], tophits[i][1][1]))
		f.write('%s\t%s\t' % (','.join(tophits[i][2][0]), ','.join(tophits[i][2][1])))
		f.write('%d\t' % tophits[i][3])
		f.write('%E\t' % tophits[i][4])
		f.write('%E\t' % tophits[i][5])
		f.write('%E\t' % tophits[i][6])
		f.write('%s\n' % tophits[i][7])
	f.close()

def subset(index, k):
	k = len(index) if k > len(index) else k
	if k == 0 or len(index) == 0:
		return [[]]

	sublist = subset(index[1:], k - 1)
	output = []
	for i in range(len(sublist)):
		combo = [index[0]]
		combo.extend(sublist[i])
		output.append(combo)

	if k <= len(index) - 1:
		sublist = subset(index[1:], k)
		for i in range(len(sublist)):
			output.append(sublist[i])

	return output

def get_true_true(result, index, param, mass):
	true_true = []
	for i in range(len(result)):
		title = result[i][0]
		spec = index.spec_dict[title]
		ch = spec.ch

		candidate = []
		sum_int = []
		for j in range(len(result[i][1])):
			pep1 = index.unique_pep[0][result[i][1][j][0][0]]
			pep2 = index.unique_pep[0][result[i][1][j][0][1]]
			sl = [set(), set()]
			pos = result[i][1][j][1]

			for pro in pep1.pro_id:
				cols = pro.split('|R')
				if len(cols) > 1 and len(cols[1]) > 0:
					sl[0].add(cols[1][0])

			for pro in pep2.pro_id:
				cols = pro.split('|R')
				if len(cols) > 1 and len(cols[1]) > 0:
					sl[1].add(cols[1][0])

			feature = list(result[i][1][j][2][0])
			feature.extend(result[i][1][j][2][1])
			if feature[0] / float(feature[7]) >= 0.20 and feature[1] / float(feature[7]) >= 0.20 and feature[8] / float(feature[15]) >= 0.20 and feature[9] / float(feature[15]) >= 0.20 and feature[2] >= 0.1 and feature[10] >= 0.1 and (len(sl[0]) == 0 or len(sl[1]) == 0 or len(sl[0].intersection(sl[1]))) > 0:
				xl = XLink(pep1, pep2, pos, ch, mass, param)
				match = Match(spec, xl, mass, param)
				match.match(mass)

				candidate.append(match)
				sum_int.append(feature[2] + feature[10])

		if len(candidate) == 0:
			continue
		combo = zip(candidate, sum_int)
		candidate = list(zip(*sorted(combo, key = lambda x : x[1], reverse = True))[0])
		sum_int = list(zip(*sorted(combo, key = lambda x : x[1], reverse = True))[1])
		true_true.append(candidate[0])

	for i in range(len(true_true)):
		pep1 = true_true[i].xlink.pep[0]
		pep2 = true_true[i].xlink.pep[1]
		s = pep1.seq + '\t' + pep2.seq + '\t' + ','.join(pep1.pro_id) + '\t' + ','.join(pep2.pro_id)
		print s

	if len(true_true) < 150:
		print '\nWARNING: The number of True-True PSMs(' + str(len(true_true)) + ') is too small and maybe insufficient for training an reliable model!\n'
	return true_true


def read_true_true(filename, mass, param):
	f = open(filename)
	lines = f.readlines()
	modif = dict(position = [], delta_mass = [])

	true_true = []
	mz = []
	it = []

	for l in lines:
		l = l[:-1]
		cols = l.split('\t')	
		if len(cols) == 5:
			pep = [cols[0], cols[1]]
			pos = [int(cols[2]) - 1, int(cols[3]) - 1]
			ch = int(cols[4])
		elif len(cols) == 2:
			mz.append(float(cols[0]))
			it.append(float(cols[1]))
		elif len(cols) == 1 and len(cols[0]) > 0:
			prec_mz = float(cols[0])	
		else:
			pep1 = Peptide(pep[0], 'protein_' + pep[0], modif, mass, False)
			pep2 = Peptide(pep[1], 'prptein_' + pep[1], modif, mass, False)
			xl = XLink(pep1, pep2, pos, ch, mass, param)
			spec = Spectrum('precursor_mz' + str(prec_mz), 0, prec_mz, ch, mz, it, 0, mass)	
			match = Match(spec, xl, mass, param)
			true_true.append(match)

			mz = []
			it = []
	f.close()

	if len(true_true) < 150:
			print '\nWARNING: The number of True-True PSMs(' + str(len(true_true)) + ') is too small and maybe insufficient for training an reliable model!\n'

def get_true_false(true_true, param, mass):
	true_true_seq = set()
	true_true_mass = []
	linker_mass = param['linker_mass']

	if 'uniprot_database' not in param:
		print 'No uniprot database specified!'
		print 'execution will terminate!'
		sys.exit(1)

	fasta = FastaReader(param['uniprot_database']).read_fasta()
	pattern_string = param['pattern_string']
	peps = dict()

	for match in true_true:
		true_true_seq.add(match.xlink.pep[0].seq)
		true_true_seq.add(match.xlink.pep[1].seq)
		true_true_mass.append(match.xlink.mol_weight - linker_mass)
	true_true_mass = max(true_true_mass)


	for header, seq in fasta:
		if 'MOUSE' in header:
			pep = get_pep_from_pro(header, seq, pattern_string, mass, param)

			for p in pep:
				if p.seq not in peps and p.seq not in true_true_seq and p.prec_mass < true_true_mass:
					peps[p.seq] = p


	peps = peps.values()
	alpha = []
	beta = []

	for i in range(len(true_true)):
		print i
		sys.stdout.flush()
		match = true_true[i]
		ch = match.spec.ch
		ms2tol = match.xlink.mol_weight * 5 * (10 ** (-6))
		alpha.append([])
		beta.append([])

		pep = match.xlink.pep

		for j in range(len(peps)):
			if abs(pep[0].prec_mass + peps[j].prec_mass + linker_mass - match.xlink.mol_weight) <= ms2tol:
				pepseq1 = pep[0].seq
				pepseq2 = peps[j].seq
				k_pos1 = list(zip(*filter(lambda x : x[1] == 'K', enumerate(pepseq1[:-1])))[0])
				k_pos2 = list(zip(*filter(lambda x : x[1] == 'K', enumerate(pepseq2[:-1])))[0])
				pos = [k_pos1[len(k_pos1) / 2], k_pos2[len(k_pos2) / 2]]
				xl = XLink(pep[0], peps[j], pos, ch, mass, param)

				tf = Match(match.spec, xl, mass, param)
				tf.match(mass)
				feature = tf.feature
				if (feature[1][0] + feature[1][1]) / float(feature[1][7]) >= 0.2:
					alpha[-1].append(tf)

		for j in range(len(peps)):
			if abs(pep[1].prec_mass + peps[j].prec_mass + linker_mass - match.xlink.mol_weight) <= ms2tol:
				pepseq1 = pep[1].seq
				pepseq2 = peps[j].seq
				k_pos1 = list(zip(*filter(lambda x : x[1] == 'K', enumerate(pepseq1[:-1])))[0])
				k_pos2 = list(zip(*filter(lambda x : x[1] == 'K', enumerate(pepseq2[:-1])))[0])
				pos = [k_pos1[len(k_pos1) / 2], k_pos2[len(k_pos2) / 2]]
				xl = XLink(pep[1], peps[j], pos, ch, mass, param)

				tf = Match(match.spec, xl, mass, param)
				tf.match(mass)
				feature = tf.feature
				if (feature[1][0] + feature[1][1]) / float(feature[1][7]) >= 0.2:
					beta[-1].append(tf)

	true_false = []
	for i in range(len(alpha)):
		true_false.extend(alpha[i])
	for i in range(len(beta)):
		true_false.extend(beta[i])

	return true_false

def get_false_false(true_true, param, mass):
	linker_mass = param['linker_mass']

	true_true_seq = set()
	true_true_mass = []
	for match in true_true:
		true_true_seq.add(match.xlink.pep[0].seq)
		true_true_seq.add(match.xlink.pep[1].seq)
		true_true_mass.append(match.xlink.mol_weight - linker_mass)

	min_mass = int(min(true_true_mass) - 0.2)
	max_mass = int(max(true_true_mass) + 0.2)

	if 'uniprot_database' not in param:
		print 'No uniprot database specified!'
		print 'execution will terminate!'
		sys.exit(1)

	fasta = FastaReader(param['uniprot_database']).read_fasta()
	pattern_string = param['pattern_string']
	peps = dict()

	for header, seq in fasta:
		if 'YEAST' in header:
			pep = get_pep_from_pro(header, seq, pattern_string, mass, param)
			for p in pep:
				if p.seq not in peps and p.seq not in true_true_seq:
					peps[p.seq] = p
	peps = peps.values()

	int_dict = dict()
	for pep in peps:
		num = int(pep.prec_mass)
		if num > max_mass:
			continue

		if num not in int_dict:
			int_dict[num] = [pep]
		else:
			int_dict[num].append(pep)

	false_false = []
	for k in range(len(true_true)):
		match = true_true[k]
		print k
		sys.stdout.flush()

		false_false.append([])
		prec_mass = match.xlink.mol_weight - linker_mass
		ch = match.spec.ch
		ms2tol = match.xlink.mol_weight * 3 * (10 ** (-6))

		mass_list = range(500, max_mass - 500)
		shuffle(mass_list)
		mass_list = mass_list[:25]

		for m in mass_list:
			if m not in int_dict:
				continue

			shuffle(int_dict[m])
			int_dict[m] = int_dict[m][:50]

			for i in range(len(int_dict[m])):
				num = int(prec_mass - int_dict[m][i].prec_mass)
				if num not in int_dict:
					continue
				shuffle(int_dict[num])
				int_dict[num] = int_dict[num][:50]

				for j in range(len(int_dict[num])):
					pepseq1 = int_dict[m][i].seq
					pepseq2 = int_dict[num][j].seq
					k_pos1 = list(zip(*filter(lambda x : x[1] == 'K', enumerate(pepseq1[:-1])))[0])
					k_pos2 = list(zip(*filter(lambda x : x[1] == 'K', enumerate(pepseq2[:-1])))[0])
					pos = [k_pos1[len(k_pos1) / 2], k_pos2[len(k_pos2) / 2]]
					xl = XLink(int_dict[m][i], int_dict[num][j], pos, ch, mass, param)
					if abs(match.xlink.mol_weight - xl.mol_weight) <= ms2tol:
						ff = Match(match.spec, xl, mass, param)
						ff.match(mass)

						feature = ff.feature
						if (feature[0][0] + feature[0][1]) / float(feature[0][7]) >= 0.15 and (feature[1][0] + feature[1][1]) / float(feature[1][7]) >= 0.15:
							false_false[-1].append(ff)
	l = []
	for i in range(len(false_false)):
		l.extend(false_false[i])
	false_false = l
	return false_false

def get_feature_matrix(matches):
	X = []
	for m in matches:
		x = []
		x.extend(m.feature[0])
		x.extend(m.feature[1])
		X.append(x)
	X = np.asarray(X)
	return X

def filter_by_fdr(top_hits, decoy_string, cutoff, is_unique):
	if cutoff < 0 or cutoff > 1:
		print 'fdr cutoff should be greater than 0.0 and less than 1.0!'
		sys.exit(1)

	top_hits = sorted(top_hits, key = lambda x : x[4], reverse = True)

	intra_cum_count = []
	inter_cum_count = []
	tardec_cum_count = []
	decdec_cum_count = []

	intra_count = 0
	inter_count = 0
	tardec_count = 0
	decdec_count = 0

	unique_intra_fdr = set()
	unique_inter_fdr = set()
	unique_tardec_fdr = set()
	unique_decdec_fdr = set()

	xl_type = []

	for i in range(len(top_hits)):
		pro1 = top_hits[i][2][0]
		pro2 = top_hits[i][2][1]

		is_tar = [[], []]
		is_dec = [[], []]
		pep_str = [top_hits[i][0][0], top_hits[i][0][1]]
#		pep_str = [top_hits[i][0][0], top_hits[i][0][1], str(top_hits[i][3])]
		pep_str = '_'.join(pep_str)

		for part in pro1:
			if decoy_string in part:
				is_dec[0].append(True)
				is_tar[0].append(False)
			else:
				is_dec[0].append(False)
				is_tar[0].append(True)

		for part in pro2:
			if decoy_string in part:
				is_dec[1].append(True)
				is_tar[1].append(False)
			else:
				is_dec[1].append(False)
				is_tar[1].append(True)

		if any(is_tar[0]) and any(is_tar[1]):
			if len(set(pro1).intersection(set(pro2))) > 0:

				if is_unique == False:
					intra_count += 1
				else:
					unique_intra_fdr.add(pep_str)
					intra_count = len(unique_intra_fdr)

				xl = 'intraxlink'
			else:

				if is_unique == False:
					inter_count += 1
				else:
					unique_inter_fdr.add(pep_str)
					inter_count = len(unique_inter_fdr)

				xl = 'interxlink'
		elif (any(is_tar[0]) and all(is_dec[1])) or (all(is_dec[0]) and any(is_tar[1])):

			if is_unique == False:
				tardec_count += 1
			else:
				unique_tardec_fdr.add(pep_str)
				tardec_count = len(unique_tardec_fdr)

			xl = 'target-decoy'
		elif all(is_dec[0]) and all(is_dec[1]):

			if is_unique == False:
				decdec_count += 1
			else:
				unique_decdec_fdr.add(pep_str)
				decdec_count = len(unique_decdec_fdr)

			xl = 'decoy-decoy'
		else:
			print 'execution will terminate!'
			sys.exit(1)

		intra_cum_count.append(intra_count)
		inter_cum_count.append(inter_count)
		tardec_cum_count.append(tardec_count)
		decdec_cum_count.append(decdec_count)
		xl_type.append(xl)

#	tmp = enumerate(xl_type)
#	tmp = filter(lambda x : x[1] == 'target-decoy' or x[1] == 'decoy-decoy', tmp)
#	print tmp[0][0]

	fdr_intra = []
	for i in range(len(top_hits)):
		if intra_cum_count[i] != 0:
			fdr = float(tardec_cum_count[i] - decdec_cum_count[i]) / intra_cum_count[i]
			fdr_intra.append([fdr, i])
		else:
			fdr_intra.append([float(sys.maxint), i])

	fdr_inter = []
	for i in range(len(top_hits)):
		if inter_cum_count[i] != 0:
			fdr = float(tardec_cum_count[i] - decdec_cum_count[i]) / inter_cum_count[i]
			fdr_inter.append([fdr, i])
		else:
			fdr_inter.append([float(sys.maxint), i])
#	pickle.dump([fdr_intra, fdr_inter], file('save.pickle', 'w'))
	fdr_intra = filter(lambda x : x[0] <= cutoff, fdr_intra)
	fdr_inter = filter(lambda x : x[0] <= cutoff, fdr_inter)

	if any(fdr_intra) < 0 or any(fdr_inter) < 0:
		print 'warning: negative fdr value'

	max_index_intra = fdr_intra[-1][1] if len(fdr_intra) > 0 else -1
	max_index_inter = fdr_inter[-1][1] if len(fdr_inter) > 0 else -1

	intra = []
	for i in range(len(top_hits)):
		if xl_type[i] == 'intraxlink' and i <= max_index_intra:
			intra.append(top_hits[i])

	inter = []
	for i in range(len(top_hits)):
		if xl_type[i] == 'interxlink' and i <= max_index_inter:
			inter.append(top_hits[i])

	print '#intra = %d, #TD = %d, #DD = %d' % (intra_cum_count[max_index_intra], tardec_cum_count[max_index_intra], decdec_cum_count[max_index_intra])
	print '#inter = %d, #TD = %d, #DD = %d' % (inter_cum_count[max_index_inter], tardec_cum_count[max_index_inter], decdec_cum_count[max_index_inter]) 

	unique_intra = set()
	f = open('intra' + str(cutoff), 'w')
	for i in range(len(intra)):
		pep = [intra[i][0][0], intra[i][0][1]]
		pro = [','.join(intra[i][2][0]), ','.join(intra[i][2][1])]
		pos = [intra[i][1][0], intra[i][1][1]]
		score = intra[i][4]
		ch = intra[i][3]
		scan = intra[i][-1]

		f.write('%d\t%s\t%s\t%d\t%d\t%s\t%s\t%d\t%f\t%s\n' % (i + 1, pep[0], pep[1], pos[0] + 1, pos[1] + 1, pro[0], pro[1], ch, score, scan))
		unique_intra.add('_'.join(pep))
	f.close()

	unique_inter = set()
	f = open('inter' + str(cutoff), 'w')
	for i in range(len(inter)):
		pep = [inter[i][0][0], inter[i][0][1]]
		pro = [','.join(inter[i][2][0]), ','.join(inter[i][2][1])]
		pos = [inter[i][1][0], inter[i][1][1]]
		score = inter[i][4]
		ch = inter[i][3]
		scan = inter[i][-1]

		f.write('%d\t%s\t%s\t%d\t%d\t%s\t%s\t%d\t%f\t%s\n' % (i + 1, pep[0], pep[1], pos[0] + 1, pos[1] + 1, pro[0], pro[1], ch, score, scan))
		unique_inter.add('_'.join(pep))
	f.close()

	return [intra, unique_intra, inter, unique_inter]

def get_matches_from_tophits(tophit, param, mass, spec_dict):
	pep1 = tophit[0][0]
	pep2 = tophit[0][1]
	pro1 = tophit[2][0]
	pro2 = tophit[2][1]
	pos1 = tophit[1][0]
	pos2 = tophit[1][1]
	ch = tophit[3]
	title = tophit[7]
	filename = pep1 + '_' + pep2 + '_' + str(pos1) + '_' + str(pos2) + '_' + str(ch) + '_' + title + '.annotation'

	modif = dict(position=[], delta_mass=[])
	for j in range(len(pep1)):
		if pep1[j].islower():
			modif['position'].append(j)
			modif['delta_mass'].append(param['mod_mass'])
	pep1 = Peptide(pep1, ', '.join(pro1), modif, mass, pos1 == 0 and pep1[0] != param['x_residue'])

	modif = dict(position=[], delta_mass=[])
	for j in range(len(pep2)):
		if pep2[j].islower():
			modif['position'].append(j)
			modif['delta_mass'].append(param['mod_mass'])
	pep2 = Peptide(pep2, ', '.join(pro2), modif, mass, pos2 == 0 and pep2[0] != param['x_residue'])
	xl = XLink(pep1, pep2, [pos1, pos2], ch, mass, param)
	match = Match(spec_dict[title], xl, mass, param)
	return [match, filename]
