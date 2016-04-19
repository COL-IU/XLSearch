class FastaReader:
	def __init__(self, fileName):
		self.fileName = fileName
	def readFasta(self):
		f = open(self.fileName)
		lines = f.readlines()
		
		fasta = []
		header = str()
		sequence = str()
		for l in lines:
			if l[0] == '>':
				if len(header) != 0:
					fasta.append((header, sequence))
					header = l.strip()[1:]
					header = header[0:header.index(' ')]
					sequence = str()
				else:
					header = l.strip()[1:]
					header = header[0:header.index(' ')]
			elif len(header) != 0:
				sequence += l.strip()
		fasta.append((header, sequence))	

		return	fasta
