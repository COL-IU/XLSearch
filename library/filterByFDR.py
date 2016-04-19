
def filterByFDR(tophits, cutoff):
	tophits = sorted(tophits, key = lambda x : x[4], reverse = True)

	intraCumCount = []
	interCumCount = []
	tardecCumCount = []
	decdecCumCount = []

	intraCount = 0
	interCount = 0
	tardecCount = 0
	decdecCount = 0

	xlType = []

	for i in range(len(tophits)):
		pro1 = tophits[i][2][0]
		pro2 = tophits[i][2][1]

		isTar = [[], []]
		isDec = [[], []]

		for part in pro1:
			if 'reverse' in part:
				isDec[0].append(True)
				isTar[0].append(False)
			else:
				isDec[0].append(False)
				isTar[0].append(True)

		for part in pro2:
			if 'reverse' in part:
				isDec[1].append(True)
				isTar[1].append(False)
			else:
				isDec[1].append(False)
				isTar[1].append(True)

		if any(isTar[0]) and any(isTar[1]):
			if len(set(pro1).intersection(set(pro2))) > 0:
				intraCount += 1
				xl = 'intraxlink'
			else:
				interCount += 1
				xl = 'interxlink'
		elif (any(isTar[0]) and all(isDec[1])) or (all(isDec[0]) and any(isTar[1])):
			tardecCount += 1
			xl = 'target-decoy'
		elif all(isDec[0]) and all(isDec[1]):
			decdecCount += 1
			xl = 'decoy-decoy'
		else:
			print '???????????'

		intraCumCount.append(intraCount)
		interCumCount.append(interCount)
		tardecCumCount.append(tardecCount)
		decdecCumCount.append(decdecCount)
		xlType.append(xl)

	fdrIntra = []
	for i in range(len(tophits)):
		if intraCumCount[i] != 0:
			fdr = float(tardecCumCount[i] - decdecCumCount[i]) / intraCumCount[i]
			fdrIntra.append([fdr, i])

	fdrInter = []
	for i in range(len(tophits)):
		if interCumCount[i] != 0:
			fdr = float(tardecCumCount[i] - decdecCumCount[i]) / interCumCount[i]
			fdrInter.append([fdr, i])

	fdrIntra = filter(lambda x : x[0] <= cutoff, fdrIntra)
	fdrInter = filter(lambda x : x[0] <= cutoff, fdrInter)

	maxIndexIntra = fdrIntra[-1][1] if len(fdrIntra) > 0 else -1
	maxIndexInter = fdrInter[-1][1] if len(fdrInter) > 0 else -1

	INTRA = []
	for i in range(len(tophits)):
		if xlType[i] == 'intraxlink' and i <= maxIndexIntra:
			INTRA.append(tophits[i])
	INTER = []
	for i in range(len(tophits)):
		if xlType[i] == 'interxlink' and i <= maxIndexInter:
			INTER.append(tophits[i])

	uniqueIntra = set()
	f = open('intra' + str(cutoff), 'w')
	for i in range(len(INTRA)):
		pep = [INTRA[i][0][0], INTRA[i][0][1]]
		pro = [','.join(INTRA[i][2][0]), ','.join(INTRA[i][2][1])]
		pos = [INTRA[i][1][0], INTRA[i][1][1]]
		score = INTRA[i][4]
		ch = INTRA[i][3]
		scan = INTRA[i][-1]

		f.write('%d\t%s\t%s\t%d\t%d\t%s\t%s\t%d\t%f\t%s\n' % (i + 1, pep[0], pep[1], pos[0] + 1, pos[1] + 1, pro[0], pro[1], ch, score, scan))
		uniqueIntra.add('_'.join(pep))
	f.close()

	uniqueInter = set()
	f = open('inter' + str(cutoff), 'w')
	for i in range(len(INTER)):
		pep = [INTER[i][0][0], INTER[i][0][1]]
		pro = [','.join(INTER[i][2][0]), ','.join(INTER[i][2][1])]
		pos = [INTER[i][1][0], INTER[i][1][1]]
		score = INTER[i][4]
		ch = INTER[i][3]
		scan = INTER[i][-1]

		f.write('%d\t%s\t%s\t%d\t%d\t%s\t%s\t%d\t%f\t%s\n' % (i + 1, pep[0], pep[1], pos[0] + 1, pos[1] + 1, pro[0], pro[1], ch, score, scan))
		uniqueInter.add('_'.join(pep))
	f.close()

	return [INTRA, uniqueIntra, INTER, uniqueInter]
