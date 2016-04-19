def readKojakResult(fileName):
	f = open(fileName)
	lines = f.readlines()
	baseName = fileName[:-10]

	uniquePPSM = dict()
	for line in lines[2:]:
		cols = line.split('\t')
		if cols[9] == '-' or cols[12] == '-':
			continue

		scan = baseName + '.' + cols[0] + '.' + cols[3]
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
	
		tmp = sorted(((pep[0], pos[0], pro[0]), (pep[1], pos[1], pro[1])), key = lambda x : x[0])
		pep[0] = tmp[0][0]
		pep[1] = tmp[1][0]
		pos[0] = tmp[0][1]
		pos[1] = tmp[1][1]
		pro[0] = tmp[0][2]
		pro[1] = tmp[1][2]
	
		ppsm = [pep, pos, pro, ch, score, scan]

		if scan not in uniquePPSM:
			uniquePPSM[scan] = ppsm
		elif uniquePPSM[scan][4] < ppsm[4]:
			uniquePPSM[scan] = ppsm

	tophits = uniquePPSM.values()
	return tophits
