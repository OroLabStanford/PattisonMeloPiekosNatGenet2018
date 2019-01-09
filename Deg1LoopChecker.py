'''
Author:		Samantha Piekos
Date:		6/14/2017
Title:		Deg1LoopCheckerFinal.py
Version:	python/3.3.2

This program identifies two distal genomic elements that are connected to each other via chromatin looping ensuring that
each element is in it's own individual contact bin. One genomic element is already associated with a loop (in one of the 
contact bins) in the HiChIP input file. The second element is contained in a bed file. The output is the coordinates and 
ID of the second element, the loop read count, fdr, & ID, and the anchored first element coordinates plus ID. Prints to
console the number of unique second elements attached to at least one of the anchored elements.

python EQTLfinder.py elementFile GTExDirectory outputDirectory tag=output

@param 	HiChIP_file		path to text file containing HiChIP loop coordinates (chr1 start1 stop1 chr2 start2 stop2)
						loop_count label, loop count, fdr, contact ID, and anchored element coordinates + ID (14 columns)
@param 	anchor_name 	text label of data type of anchored element
@param 	target_file		path to bed file of second genomic element of interest
@param 	target_name		text label of data type of the second element
@param 	output_file		path to output file to write the second element coordinates + ID, loop count + fdr + ID, and the
						anchored elements coordinates + ID
'''

import sys
from collections import OrderedDict

'''
Input is a key, value, and dictionary. Writes entry to a list that is the value of the given key
only if it is not already in the list. Returns the updated dictionary

@key	key that the entry is to be associated with
@entry	variable to be appended to the value list associated with the given key
@dictionary	dictionary to which the entry is to be written
@return 	returns the updated dictionary
'''
def write2dict(key, value, dictionary):
# writes only unique entries to a dictionary in which the value is a set to which the entry is appended
	if key in dictionary:
		if value not in dictionary[key]:
			dictionary[key].append(value)
	else:
		dictionary[key] = [value]
	return(dictionary)

'''
Input is a file path and an empty dictionary. Writes the input file to the dictionary with the
chromosome as the key and the value as a list of lists containing the line in list format. 
Returns the updated dictionary.

@file	file path to text file to be written to the ditctionary
@dictionary	dictionary to which the entry is to be written
@return 	updated dictionary
'''
def unpackFile2ChrDict(file, dictionary):
	with open(file, 'r') as file:  # write input file line by line to dictionary with chrom as key
		for line in file:
			line = line.rstrip('\r\n').split('\t')
			chrom, line[1], line[2] = line[0], int(line[1]), int(line[2])
			dictionary = write2dict(chrom, line, dictionary)
	return(dictionary)

'''
Returns a boolean concerning if a given integer falls within a given range.

@val	integer to be checked if its in the bin
@bin	a list or tuple of two integer or floats that are the [start, stop] of a bin
@return	boolean regarding if the val is between the two numbers defining the bin
'''
def checkBin(val, Bin):
	start, stop = int(Bin[0]), int(Bin[1])
	val = int(val)
	return(val < stop and val > start)

'''
Returns boolean regarding if a given range of integers is fully contained
within a second range of integers.

@start	integer of the lower bound parameter of the element to be checked
@stop	integer of the upper bound parameter of the element to be checked
@Bin	a list or tuple of two integer or floats that are the [start, stop] of a bin
@return	boolean regarding if the element is fully contained within a bin
'''
def middleBin(start, stop, Bin):
	bin_start, bin_stop = int(Bin[0]), int(Bin[1])
	start, stop = int(start), int(stop)
	return(start <= bin_start and stop >= bin_stop)

'''
Returns boolean regarding if a given range of integers partially or completely
overlaps with a second range of integers.

@start	integer of the lower bound parameter of the element to be checked
@stop	integer of the upper bound parameter of the element to be checked
@Bin	a list or tuple of two integer or floats that are the [start, stop] of a bin
@return	boolean regarding if the element partially or completely overlaps with the given bin
'''
def binChecker(start, stop, Bin):
# checks if an element is contained at all within a given bin - returns boolean
	start, stop = int(start), int(stop)
	Bin[0], Bin[1] = int(Bin[0]), int(Bin[1])
	# return HiChIP bin coordinates if start or stop of the peak coordinates falls within one of the bins
	if checkBin(start, Bin) or checkBin(stop, Bin):
		return(True)
	# return HiCHIP coordinates if the peak overlaps with the entirity of the loop bins
	elif middleBin(start, stop, Bin):
		return(True)
	else:  # return False if there is no overlap between the peak and either of the loop bins
		return(False)

'''
Input is a list containing the genomic start stop coordinates of an element, another list containing the 
genomic start stop coordinates of a second element, and the last containing the genomic coordinates of a loop
(chr1 start1 stop1 chr2 start2 stop2). Assumes that the loop is intra-chromasomal and that both elements are
known to be contained on the same chromosome as the loop. Returns True if the first element is in one contact 
bin of the loop and the second element is on the other contact bin of the loop. Otherwise it returns False.

@feature 	list containing genomic start and stop coordinates (but not chr) of the first element of interest
@anchor 	list containing genomic start and stop coordinates (but not chr) of the second element of interest
@loop 		genomic coordinates of the loop (chr1 start1 stop1 chr2 start2 stop2)
'''
def DistalConnectCheck(feature, anchor, loop):
# identifies which end of the loop the anchor is associated and checks if the feature is on the other end of the loop
	feature_start, feature_stop = int(feature[0]), int(feature[1])
	anchor_start, anchor_stop = int(anchor[0]), int(anchor[1])
	loop1 = [int(loop[1]), int(loop[2])]
	loop2 = [int(loop[4]), int(loop[5])]

	check = binChecker(anchor_start, anchor_stop, loop1)
	if check:
		isin0 = binChecker(feature_start, feature_stop, loop2)
		if isin0:
			return('yes')
	else:  # anchor is in bin 2 and check if feature is in bin 1
		isin1 = binChecker(feature_start, feature_stop, loop1)
		if isin1:  # write lines of HiChIP files to dictionary if they are anchored at one end with the feature of interest
			return('yes')
	return(False)

'''
Input is a dictionary containing as values a list of lists, in which the inner list contains an ID in the
index[3] position. The number of unique IDs in this index[3] position is counted and returned.
@dictionary 	dictionary containing list of lists whose inner list has an ID to be counted in index[3]
@return 	the number of unique IDs in the dictionary
'''
def countUniqueID(dictionary):
	UniqueID = set()
	for key, value in dictionary.items():  # loop over dictionary adding new IDs to the UniqueID set
		for item in value:
			UniqueID.add(item[3])  # ID has to be in index[3] (fourth column) position in the loop list
	return(len(UniqueID))

'''
Input is the dictionaries containing the info from the HiChIP file and second element file as well as the identities of the 
elements in each of these files. Identify contacts in which the anchor element is in one contact bin and the target is in the
other contact bin of a chromosome loop. If it is write the second element coordinates + ID, loop count + fdr + ID, and anchor
element coordinates + ID to a dictionary with the chromosome as the key and list containing that info in a list as the value.
Return the dictionary.

@HiChIP_dict 	dictionary with HiChIP input file info with the chromosome as the key and a list of lists as the value
@anchor_name 	the type of data the anchor element is in the HiChIP file
@target_dict	dictionary with second element input file info with the chromosome as the key and a list of lists as the value
@target_name 	the type of data the element in the second element file
@return 	dictionary containing desired output info
'''
def deg1Analysis(HiChIP_dict, anchor_name, target_dict, target_name):
	output_dict = {}
	for chrom, value in HiChIP_dict.items():  # loop over HiChIP_dict and target_dicts
		if chrom in target_dict:
			for item in value:  # define elements start stop and loop coordinates
				loop = item[:6]
				anchor = [int(item[11]), int(item[12])]
				for line in target_dict[chrom]:  # check if one element is in one contact bin and the other is the other bin
					feature = [line[1], line[2]]
					check = DistalConnectCheck(feature, anchor, loop)
					if check:  # write info of elements and loop to dictionary if it's in the correct configuration
						keep = line[0:4]
						keep.extend(item[6:14])
						output_dict = write2dict(chrom, keep, output_dict)
	return(output_dict)

'''
Orders the list of values in a dictionary by the second and third items in the inner lists. Returns sorted dictionary.

@chr_dict	dictionary whose values are a list of lists
@return 	dictionary value lists sorted by the second and third items in the lists contained within the value lists
'''
def orderChrDict(chr_dict):  # orders keys in dictionary and sorts values by second component of values
	for key, value in chr_dict.items():  # sort values for each key
		chr_dict[key] = sorted(value, key=lambda element: (element[1], element[2]))
	final_dict = OrderedDict(sorted(chr_dict.items(), key=lambda t: t[0]))  # sort keys and write sorted info to an Ordered Dict

	return(final_dict)

'''
Writes a given dictionary containing lists of string lists to a given output file.

@output_dict	dictionary containing as the values lists of lists in which the inner lists are lines to be written to the output
@file 	string of the file path to which to write the output
'''
def outputFile(output_dict, file):
	with open(file, 'w') as file:
		# parse through dictionary writing string lists in all the value lists to the output file
		for chrom, value in output_dict.items():
			for entry in value:
				line = []
				for item in entry:
					item = str(item)
					line.append(item)
				file.write('\t'.join(line) + '\n')
	return('done')

def main():
	HiChIP_file = sys.argv[1]  # write input parameters to variables
	anchor_name = sys.argv[2]
	target_file = sys.argv[3]
	target_name = sys.argv[4]
	output_file = sys.argv[5]
	HiChIP_dict = {}
	target_dict = {}
	unpack = [(HiChIP_file, HiChIP_dict), (target_file, target_dict)]  # pair together file path and associated empty data structure
	for file, dictionary in unpack:  # unpack input files to dictionaries with chrom as key
		dictionary = unpackFile2ChrDict(file, dictionary)
	output_dict = deg1Analysis(HiChIP_dict, anchor_name, target_dict, target_name)  # identify contacts with target in opposite bin as the anchor
	count = countUniqueID(output_dict)  # 
	print('# Number of contacts in which', target_name, 'is directly looped to', anchor_name, '(i.e. 1°', anchor_name, ') =', count)
	sort_output_dict = orderChrDict(output_dict)  # order the output data by genomic location in each on

	print(outputFile(sort_output_dict, output_file))  # write contacts of interest to new output file (subset of original input file)

if __name__ == '__main__':
	main()
