def read_kojak_result(filename):
	f = open(filename)
	lines = f.readlines()
	basename = filename[:-10]

	unique_psm = dict()
	for line in lines[2:]:
		cols = line.split('\t')
		if cols[9] == '-' or cols[12] == '-':
			continue

		scan = basename + '.' + cols[0] + '.' + cols[3]
		score = float(cols[6])
		ch = int(cols[3])
		pep = [cols[9], cols[12]]
		pos = [int(cols[10]) - 1, int(cols[13]) - 1]
		pro = [cols[11], cols[14]]
		pro[0] = pro[0].split(';')[:-1]
		pro[1] = pro[1].split(';')[:-1]

		for i in range(len(pro[0])):
			pro[0][i] = pro[0][i].split(' ')[0][1:]
		for i in range(len(pro[1])):
			pro[1][i] = pro[1][i].split(' ')[0][1:]
	
		tmp = sorted([[pep[0], pos[0], pro[0]], [pep[1], pos[1], pro[1]]], key = lambda x : x[0])
		pep[0] = tmp[0][0]
		pep[1] = tmp[1][0]
		pos[0] = tmp[0][1]
		pos[1] = tmp[1][1]
		pro[0] = tmp[0][2]
		pro[1] = tmp[1][2]
	
		psm = [pep, pos, pro, ch, score, scan]

		is_tar = [[], []]
		is_dec = [[], []]

		for part in pro[0]:
			if 'reverse' in part:
				is_dec[0].append(True)
				is_tar[0].append(False)
			else:
				is_dec[0].append(False)
				is_tar[0].append(True)

		for part in pro[1]:
			if 'reverse' in part:
				is_dec[1].append(True)
				is_tar[1].append(False)
			else:
				is_dec[1].append(False)
				is_tar[1].append(True)

		if any(is_tar[0]) and any(is_tar[1]):
			psm.append('TT')
		elif (any(is_tar[0]) and all(is_dec[1])) or (all(is_dec[0]) and any(is_tar[1])):
			psm.append('TD')
		elif all(is_dec[0]) and all(is_dec[1]):
			psm.append('DD')

		if scan not in unique_psm:
			unique_psm[scan] = psm
		elif unique_psm[scan][4] <= psm[4]:
			if psm[-1] == 'TT' and (unique_psm[scan][-1] == 'TD' or unique_psm[scan][-1] == 'DD'):
				unique_psm[scan] = psm
			elif psm[-1] == 'DD' and unique_psm[scan][-1] == 'TD':
				unique_psm[scan] = psm

	tophits = unique_psm.values()
	for i in range(len(tophits)):
		tophits[i].pop()

	return tophits
