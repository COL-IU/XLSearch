

def printPPSM(hit, specdict, fileName):
	pep = hit[0]
	pos = hit[1]
	pro = hit[2]
	ch = hit[3]
	scan = hit[-1]
	mz = specdict[scan].mz
	it = specdict[scan].it

	f = file(fileName, 'w')
	f.write('%s\n' % scan)
	f.write('%s\t%s\n' % (pep[0], pep[1]))
	f.write('%s\t%s\n' % (pro[0], pro[1]))
	f.write('%d\t%d\n' % (pos[0], pos[1]))
	f.write('%d\n' % ch)

	for i in range(len(mz)):
		f.write('%f\t%f\n' % (mz[i], it[i]))
	f.close()

