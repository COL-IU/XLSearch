import xml.etree.ElementTree as et

def readXQuestResult(fileName):

	xmlObj = et.parse(fileName)
	root = xmlObj.getroot()
	spectrum_search = root.getchildren()	

	tophits = []

	for spectrum in spectrum_search:
		hit = spectrum.getchildren()
		if len(hit) > 0:
			info = hit[0].attrib
			type = info['type']
			if type != 'xlink':
				continue

			pep = [info['seq1'], info['seq2']]
			pos = info['xlinkposition'].split(',')
			pos[0] = int(pos[0]) - 1
			pos[1] = int(pos[1]) - 1
			pro = [info['prot1'], info['prot2']]
			pro[0] = sorted(pro[0].split(','))
			pro[1] = sorted(pro[1].split(','))

			tmp = sorted(((pep[0], pos[0], pro[0]), (pep[1], pos[1], pro[1])), key = lambda x : x[0])
			pep[0] = tmp[0][0]
			pep[1] = tmp[1][0]
			pos[0] = tmp[0][1]
			pos[1] = tmp[1][1]
			pro[0] = tmp[0][2]
			pro[1] = tmp[1][2]

			ch = int(info['charge'])
			score = float(info['score'])
			scan = spectrum.attrib['spectrum']

			scan = scan[:(len(scan) - 1)/2]
			scan = scan.split('.')
			scan = scan[0] + '.' + str(int(scan[1])) + '.' + scan[3]

			one = [pep, pos, pro, ch, score, scan]

			tophits.append(one)

	return tophits
