import xml.etree.ElementTree as et
import base64
import struct
from spectrum import Spectrum

class MZXMLReader:
	def __init__(self, filename):
		self.filename = filename
		self.basename = filename[:filename.index('.')].split('/')[-1]
	def get_spec_list(self, mass, param):
		filename = self.filename
		basename = self.basename

		base_peak_int = param['base_peak_int']
		dynamic_range = param['dynamic_range']

		xml_obj = et.parse(filename)
		root = xml_obj.getroot()
		children = root.getchildren()
		children = children[0].getchildren()

		spec = []

		for i in range(0, len(children)):
			if children[i].tag[-4:] != 'scan':
				continue

			scan_num = children[i].attrib['num']
			ret_time = int(float(children[i].attrib['retentionTime'][2:-1]))

			info = children[i].getchildren()
			for j in range(0, len(info)):
				if info[j].tag[-11:] == 'precursorMz':
					if 'precursorCharge' in info[j].attrib :
						ch = int(info[j].attrib['precursorCharge'])
					else:
						ch = 4
					prec_mz = float(info[j].text)
				elif info[j].tag[-5:] == 'peaks':
					base64Peaklist = info[j].text
					data = base64.b64decode(base64Peaklist)
					if len(data) % 8 != 0:
						print 'MZXMLReader: incorrect format of peak content'
					num_peaks = len(data) / 8

					mz = []
					it = []
					for k in range(0, num_peaks):
						val = data[(k * 8 + 0) : (k * 8 + 4)]
						val = val[::-1]
						mz.append(struct.unpack('f', val)[0])
						val = data[(k * 8 + 4) : (k * 8 + 8)]
						val = val[::-1]
						it.append(struct.unpack('f', val)[0])

					max_int = max(it)

					peaks = zip(mz, it)
					peaks = filter(lambda x:x[1] >= dynamic_range * max_int, peaks)
					peaks = zip(*peaks)
					mz = list(peaks[0]);
					it = list(peaks[1]);
					it = map(lambda x : x * base_peak_int / (max_int), it)

			title = basename + '.' + scan_num + '.' + str(ch)
			spec.append(Spectrum(title, scan_num, prec_mz, ch, mz, it, ret_time, mass))
		return spec
