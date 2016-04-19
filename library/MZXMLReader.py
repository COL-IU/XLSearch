import xml.etree.ElementTree as et
import base64
import struct
from Spectrum import *

class MZXMLReader:
	def __init__(self, fileName):
		self.fileName = fileName
		self.baseName = fileName[:fileName.index('.')].split('/')[-1]
	def getSpectraList(self, mass, param):
		fileName = self.fileName
		baseName = self.baseName

		basepeakInt = param['basepeakint']
		dynamicRange = param['dynamicrange']

		xmlObj = et.parse(fileName)
		root = xmlObj.getroot()
		children = root.getchildren()
		children = children[0].getchildren()

		spectra = []

		for i in range(0, len(children)):
			if children[i].tag[-4:] != 'scan':
				continue

			scanNum = children[i].attrib['num']
			retentionTime = int(float(children[i].attrib['retentionTime'][2:-1]))

			info = children[i].getchildren()
			for j in range(0, len(info)):
				if info[j].tag[-11:] == 'precursorMz':
					ch = int(info[j].attrib['precursorCharge'])
					precursorMZ = float(info[j].text)
				elif info[j].tag[-5:] == 'peaks':
					base64Peaklist = info[j].text
					data = base64.b64decode(base64Peaklist)
					if len(data) % 8 != 0:
						print 'MZXMLReader: incorrect format of peak content'
					numPeaks = len(data) / 8

					mz = []
					it = []
					for k in range(0, numPeaks):
						val = data[(k * 8 + 0) : (k * 8 + 4)]
						val = val[::-1]
						mz.append(struct.unpack('f', val)[0])
						val = data[(k * 8 + 4) : (k * 8 + 8)]
						val = val[::-1]
						it.append(struct.unpack('f', val)[0])

					maxInt = max(it)

					peaks = zip(mz, it)
					peaks = filter(lambda x:x[1] >= dynamicRange * maxInt, peaks)
					peaks = zip(*peaks)
					mz = list(peaks[0]);
					it = list(peaks[1]);
					it = map(lambda x : x * basepeakInt / (maxInt), it)

			title = baseName + '.' + scanNum + '.' + str(ch)
			spectra.append(Spectrum(title, scanNum, precursorMZ, ch, mz, it, retentionTime, mass))
		return spectra
