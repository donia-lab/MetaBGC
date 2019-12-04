

class HMMRecord:
	def __init__(self, acc_numb, sampleType, sampleID, protType, bitscore, window, interval):
		self.acc_numb = acc_numb
		self.sampleType = sampleType
		self.sampleID = sampleID
		self.protType = protType
		self.bitscore = bitscore
		self.window = window
		self.interval = interval

class HMMFile:
	def __init__(self,intervalStart, intervalEnd, hmmFile):
		self.intervalStart = intervalStart
		self.intervalEnd = intervalEnd
		self.hmmFile = hmmFile