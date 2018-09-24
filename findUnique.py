#!/usr/bin/env python3
# Name: Jasrajveer Malhi (jmalhi)


"""
Input a FASTA file in command line, an example, 
"python3 findUnique.py <bos-tRNA.fa ". 
Program will read thefasta file and output unique 
subsequences for a single tRNA, output will be sorted by tRNA 
header.

"""

import sys
class FastAreader :
	"""FastAreader will open the fasta file."""
	def __init__ (self, fname=''):
		'''contructor: saves attribute fname '''
		self.fname = fname
	def doOpen (self):
		if self.fname is '':
			return open(self.fname)
		else:
			return sys.stdin
	def readFasta (self):
		header = ''
		sequence = ''
		with self.doOpen() as fileH:
			header = ''
			sequence = ''
			
			line = fileH.readline()
			while not line.startswith('>') :
				line = fileH.readline()
			header = line[1:].rstrip()
			for line in fileH:
				if line.startswith ('>'):
					yield header,sequence
					header = line[1:].rstrip()
					sequence = ''
				else :
					sequence += ''.join(line.rstrip().split()).upper()
		yield header,sequence

class TRNA:
	"""
	TRNA class will find unique subsets of the tRNA sequence and hold each one it finds. 

	"""
	def __init__(self,head,tails):
		"""Initialize all of the sets being used. uniqueSet will store sequences 
		that will be printed later in program."""
		self.front = head
		self.sequence = tails
		self.powerSet = set()
		self.uniqueSet = set()
		self.deleteSet = set()
		self.removedSet = set()
		start = 0
		
		
		while start != len(self.sequence): 
		# Will make a powerset, starts at begining of sequence; while not equal to the length.	
		
			for end in range(len(self.sequence), start, - 1):
				kMer = self.sequence[start : end]
				self.powerSet.add(kMer)
				# Goes from big to small and adds to the powerset
			start += 1  
	
	def smallestLength(self):
		"""The class will find the shortest value for the tRNA sequence. Each sequence found
		will be added to the deleted set."""
		small = sorted(list(self.uniqueSet), key = len)
		count = 0
		# While goes until the count is equal to the length
		while (count != len(small)): 
			index = self.sequence.find(small[count])
			
			if (len(small[count]) + index ) < len(self.sequence):
				right = self.sequence[index:len(small[count]) + index + 1]
				self.deleteSet.add(right) 
			# Both add to deleteSet
			if index > 0:
				left = self.sequence[index - 1: len(small[count]) + index]
				self.deleteSet.add(left) 
				
			count += 1 #Will add 1 to the count
		self.uniqueSet -= self.deleteSet # Sets the values for the uniqueSet

def main ():
	"""Main will load the stdin and print the unique subsets in the correct format. """
	tRNA = {} # Set dictionary called tRNA
	myReader = FastAreader(sys.stdin)
	for front, sequence in myReader.readFasta():
		name = TRNA(front, sequence)
		tRNA [front] = name
	for key, value in tRNA.items():
		for that, thatValue in tRNA.items():
			if key != that:
				value.removedSet = value.removedSet | (value.powerSet & thatValue.powerSet)
	for sets in tRNA.values():
		sets.uniqueSet = sets.powerSet - sets.removedSet

		sets.smallestLength()
		filteredSequence = sets.sequence.replace('-', '') # Filters out the dashes
		Filterunder = filteredSequence.replace('_', '') # Filters out the underscore
		
		print(sets.front, '\n', Filterunder) # Prints for both the front and sequence
		total = [] # Set a list called total
		for i in sorted(list(sets.uniqueSet)):
			index = sets.sequence.find( i )
			total.append('.' * index + i) #Adds to the end of total
			
			# Makes the program print through all the tRNA sequences
		for t_RNA in sorted(total, key = len):
			print(t_RNA.replace('-','').replace('_','').replace('?', ''))
main()