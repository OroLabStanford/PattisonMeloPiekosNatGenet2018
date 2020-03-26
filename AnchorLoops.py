"""
Author:		Samantha Piekos
Date: 		3/28/17
Title: 		AnchorLoops.py
Version: 	python/3.3.2

Filters HiChIP data for loops that contain a feature of interest in at least one of their contact bins. HiChIP_file and 
feature_of_interest_file need to be tab deliminated. The first  6 columns for the HiChIP_file must contain [chr1, start1,
stop1, chr2, start2, stop2]. The first 4 columns of the feature_of_interest file must contain [chrom, start, stop, ID]. 
Header not expected in either file. Outputs entire line of the HiChIP file containing loop anchored in the feature of 
interest with the chr start stop ID of the corresponding feature of interest added to the end.

python3 anchorLoops.py HiChIP_file feature_of_interest_file output_file

@param 	HiChIP_file		path to text file containing HiChIP loop coordinates (chr1 start1 stop1 chr2 start2 stop2)
				plus additional information such as loop count and fdr in the remaining columns
@param 	feature_of_interest_file		path to bed file of genomic element of interest
@param 	output_file		path to output file to write the HiChIP line and the coordinates + ID of the anchored
						element of interest
"""

import sys
import subprocess
import re
from collections import OrderedDict

def write2dict(key, value, dictionary):
	'''
	Input is a key, value, and dictionary. Writes unique entries to a list that is the value of the 
	given key. Returns the updated dictionary

	@key	key that the entry is to be associated with
	@value	variable to be appended to the value list associated with the given key
	@dictionary	dictionary to which the entry is to be written
	@return 	returns the updated dictionary
	'''

	if key in dictionary:
		if value not in dictionary[key]:
			dictionary[key].append(value)
	else:
		dictionary[key] = [value]
	return(dictionary)

def unpackHiChIPfile(file):
	'''
	Input is a HiChIP file path and an empty dictionary. Writes the input file to the dictionary with
	the chromosome as the key and the value as a list of lists containing the line in list format. 
	Returns the updated dictionary.

	@file	file path to text file to be written to the ditctionary
	@return 	dictionary containing the data in the input file
	'''
	dictionary = {}
	with open(file, 'r') as file:  # write input file line by line to dictionary with chrom as key
		for line in file:
			line = line.rstrip('\r\n').split('\t')
			chrom, line[1], line[2], line[4], line[5] = line[0], int(line[1]), int(line[2]), int(line[4]), int(line[5])
			if chrom[0:3] != 'chr':  # add chr in front of the chrom # if it's not already there
				line[0] = 'chr' + line[0]
				line[3] = 'chr' + line[3]
			dictionary = write2dict(line[0], line, dictionary)
	file.close()
	return(dictionary)

def checkBin(val, Bin):
	'''
	Returns a boolean concerning if a given integer falls within a given range.

	@val	integer to be checked if its in the bin
	@Bin	a list or tuple of two integer or floats that are the [start, stop] of a bin
	@return	boolean regarding if the val is between the two numbers defining the bin
	'''
	start, stop = int(Bin[0]), int(Bin[1])
	val = int(val)
	return(val <= stop and val >= start)

def middleBin(start, stop, Bin):
	'''
	Returns boolean regarding if a given range of integers is fully contained
	within a second range of integers.

	@start	integer of the lower bound parameter of the element to be checked
	@stop	integer of the upper bound parameter of the element to be checked
	@Bin	a list or tuple of two integers that are the [start, stop] of a bin
	@return	boolean regarding if the element is fully contained within a bin
	'''
	bin_start, bin_stop = Bin[0], Bin[1]
	return(start <= bin_start and stop >= bin_stop)

def loopChecker(start_feature, stop_feature, HiChIP_coordinates):
	'''
	Takes in the lower and upper bounds of an element to be checked as well as a list of the lower and upper
	bounds of 2 chromatin loop contact bins [chr1, start1, stop1, chr2, start2, stop2]. Assumes that all regions 
	are on the same chromosome for both element and chromatin contact bins. Returns boolean if any part of the
	element is contained within either of the loop contact bins (anchored in the chromatin contact).

	@start_feature	integer of the lower bound parameter of the element to be checked
	@stop_feature	integer of the upper bound parameter of the element to be checked
	@HiChIP_coordinates		a list or tuple of 4 integer or floats that are the [start1, stop1, start2, stop2]
							of the two contact bins
	@return	boolean of if the element is in either of the two loop contact bins
	'''
	HiChIP1 = [HiChIP_coordinates[1], HiChIP_coordinates[2]]
	HiChIP2 = [HiChIP_coordinates[4], HiChIP_coordinates[5]]
	# return HiChIP bin coordinates if start or stop of the peak coordinates falls within one of the bins
	if checkBin(start_feature, HiChIP1) or checkBin(stop_feature, HiChIP1) or checkBin(start_feature, HiChIP2)\
	or checkBin(stop_feature, HiChIP2):
		return(True)
	# return HiCHIP coordinates if the peak overlaps with the entirity of either of the loop bins
	elif middleBin(start_feature, stop_feature, HiChIP1) or middleBin(start_feature, stop_feature, HiChIP2):
		return(True)
	return(False)  # return False if there is no overlap between the peak and either of the loop bins

def fileLineCounter(file):
	"""
	Counts and returns the number of lines in a file provided a file path.

	@file 	file path for file to count the number of lines in
	@return 	the number of lines in the provided file
	"""
	p = subprocess.Popen(["wc", "-l", file], stdout=subprocess.PIPE)  # unix command to count number of lines in file
	out, err = p.communicate()
	decoded = out.decode('utf-8')  # decode output from unix shell
	count = re.search(r"[0-9]+", str(decoded))  # non-greedy search for the line count
	return(count)

def checkAllLoops(line, HiChIP_dict, output_dict, anchored_features, anchored_loops):
	"""
	For a given feature identifies all loops in which it's anchored in one of the contact bins. Writes these loop and
	features to the output_dict and adds the feature and loop to the anchored_features and anchored_loops set respectively.
	Reutrns updated output_dict, anchored_features, and anchored_loops.
	
	@line 	a line from the features_file contained in a list
	@HiChIP_dict 	dictionary containing the HiChIP_file data with chromosome as key and value a list of list of the data
	@output_dict 	dictionary to which write info about the loops + feature pairs
	@anchored_features 	set containing genomic coordinates + ID of features that are contained in a chromatin loop
	@anchored_loops 	set containing the loop genomic coordinates of loops that contain at least one feature
	@return 	updated output_dict
	@return 	updated anchored_features
	@return 	updated anchored_loops
	"""
	chrom, start, stop = line[0], line[1], line[2]
	if chrom in HiChIP_dict.keys():  # check if feature_of_interest chr is represented in the HiChIP data
		for value in HiChIP_dict[chrom]:  # loop through the HiChIP values associated with the chr of interest
			loop = value[:6]
			check = loopChecker(start, stop, loop)  # check if feature peak overlaps with either loop bin
			if check:  # write lines of HiChIP files to dictionary if they are anchored at one end with the feature of interest
				write = value[:]
				write.extend(line[:4])
				anchored_features.add(tuple(line[3]))
				anchored_loops.add(tuple(loop))
				output_dict = write2dict(chrom, write, output_dict)
	return(output_dict, anchored_features, anchored_loops)

def identifyAnchoredLoops(feature_file, HiChIP_dict):
	"""
	Takes in feature of interest filepath and HiChIP data contained in a dictionary. Identifies loops that contain features
	in one of the bins and writes these loop + feature pairs. Writes these to an output dictionary in which the chromosome 
	is the key and the values are a list of lists with the inner list containing the loop line plus the paired feature genomic 
	coordinates + ID. Returns this dictionary. Prints to console the number of unique features anchored in loops and number 
	of unique loops that have  at least one feature anchored in one of it's contact bins.

	@file 	file path to the feature_file
	@HiChIP_dict 	dictionary containing the HiChIP_file data with chromosome as key and value a list of list of the data
	@return 	dictionary containing chromosome as the key and the value as a list of list of loop line and paired feature 
				genomic coordinates + ID
	"""
	output_dict, anchored_features, anchored_loops = {}, set(), set()
	with open(feature_file, 'r') as feature_file:  # loop through the feature_of_interest file
		for line in feature_file:
			line = line.strip().split('\t')
			keep = line[:4]
			line[1], line[2] = int(line[1]), int(line[2])
			checkAllLoops(line, HiChIP_dict, output_dict, anchored_features, anchored_loops)
	feature_file.close()
	print('# Number of anchored features =', len(anchored_features))
	print('# Number of anchored loops =', len(anchored_loops))
	return(output_dict)

def orderChrDict(chr_dict):
	'''
	Orders the list of values in a dictionary by the second and third items in the inner lists. Also orders the keys.
	Returns sorted dictionary.

	@chr_dict	dictionary whose values are a list of lists
	@return 	dictionary value lists sorted by the second and third items in the lists contained within the value lists
	'''
	for key, value in chr_dict.items():  # sort values for each key
		chr_dict[key] = sorted(value, key=lambda element: (element[1], element[2]))
	final_dict = OrderedDict(sorted(chr_dict.items(), key=lambda t: t[0]))  # sort keys and write sorted info to an Ordered Dict

	return(final_dict)

def outputFile(output_dict, file):
	'''
	Writes a given dictionary containing lists of string lists to a given output file.

	@output_dict	dictionary containing as the values lists of lists in which the inner lists are lines to be written to the output
	@file 	string of the file path to which to write the output
	'''
	with open(file, 'w') as file:
		# parse through dictionary writing string lists in all the value lists to the output file
		for chrom, value in output_dict.items():
			for entry in value:
				line = []
				for item in entry:
					item = str(item)
					line.append(item)
				file.write('\t'.join(line) + '\n')
	file.close()
	return('done')

def main():
	HiChIP_file = sys.argv[1]  # assign terminal input to variables
	feature_file = sys.argv[2]
	output_file = sys.argv[3]
	HiChIP_dict = unpackHiChIPfile(HiChIP_file)  # write HiChIP file to a dictionary
	output_dict = identifyAnchoredLoops(feature_file, HiChIP_dict)  # identify loop and feature pairs
	sort_output_dict = orderChrDict(output_dict)
	outputFile(sort_output_dict, output_file)  # write loop and feature pairs to an output file
	count = fileLineCounter(output_file)
	if count:  # print the number of loops anchored at one end by the feature of interest
		print('# Number of lines = {}'.format(count.group(0)))

if __name__ == '__main__':
	main()
