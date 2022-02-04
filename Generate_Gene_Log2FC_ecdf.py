""" 
title: Generate_Gene_Log2FC_ecdf.py
author: Samantha Piekos
date: 2/3/22
version: Python 3.7.9
description: Makes ecdf comparing log2 FC in gene expression level of a subset of genes
			 compared to all protein-coding genes. It then evaluates via student's two-tail
			 t-test if there is a significant difference in the expression level changes of
			 the subset of genes compared to all protein-level genes.

python3 Generate_Gene_Log2FC_ecdf.py file_input_gene_fc file_input_subset_gene_list file_output_pdf subset_gene_name color=color

@file_input_gene_fc				csv file of the log2 fold change in gene expression level
@file_input_subset_gene_list	file with one gene symbol per line of subset of genes of interest
@file_output_pdf				filepath to save output ecdf graph as a pdf
@subset_gene_name				name of subset of genes (for print out statements)
@color							specify color for subset of colors in graph; string of the appropriate 
								color denoter for plt.plot; optional arguement; default value is green ('g')
"""

import sys
import numpy as np  # version=1.20.2
import matplotlib.pyplot as plt  # version=3.4.1
from scipy import stats as stats  # version=1.6.2
plt.rcParams['pdf.fonttype'] = 42


# define universal variables
X_LIM_MIN = -2
X_LIM_MAX = 2


def get_optional_arguments(optional_arguments):
	for item in optional_arguments:
		command, value = item.split('=')
		if command == 'color':
			c = value
	return c


def add_fc_to_dict(d, k, v):
	if k in d:
		d[k].append(v)
	else:
		d[k] = [v]
	return d


def average_values(d):
	# average values for multiple log2 fold changes reported for the same gene symbol
	d_final = {}
	for k, v in d.items():
		v_final = sum(v)/len(v)
		d_final[k] = v_final
	return d_final


def write_fc_to_dict(file_input_gene_fc):
	dict_fc = {}
	f = open(file_input_gene_fc, 'r')
	f.readline()
	for l in f:
		l = l.rstrip('\r\n').split(',')
		fc, gene_type, gene_symbol = l[2], l[10], l[11]
		gene_type = gene_type.strip('"')
		gene_symbol = gene_symbol.strip('"')
		if fc == 'NA':  # replace fc missing values with 0
			fc = 0
		fc = float(fc)
		if gene_type == 'protein_coding':
			dict_fc = add_fc_to_dict(dict_fc, gene_symbol, fc)
	return average_values(dict_fc)


def read_subset_file_to_list(file_input_subset_gene_list):
	list_subset_genes = []
	f = open(file_input_subset_gene_list, 'r')
	for l in f:
		l = l.rstrip('\r\n').split(',')
		list_subset_genes.append(l[0])
	return list_subset_genes


def make_fc_lists(file_input_gene_fc, file_input_subset_gene_list):
	# write protein-coding genes fold change to dict
	dict_fc = write_fc_to_dict(file_input_gene_fc)

	# generate fc lists for subset of genes and all protein-coding genes
	list_all_fc = list(dict_fc.values())
	list_subset_genes = read_subset_file_to_list(file_input_subset_gene_list)
	list_subset_fc = list({k: v for k, v in dict_fc.items() if k in list_subset_genes}.values())

	return list_subset_fc, list_all_fc


def make_ecdf(output_file, *data):
	'''
	Makes a ecdf. It takes as many datasets as input as desired.
	Input for each dataset in data list is  (list of fac values, color).
	Color specification for each dataset must be a string of the appropriate 
	color denoter for plt.plot.
	'''
	num_bins = 20  # number of bins for log2FC ecdf

	for item in data:
		entry, color = item
		sort_entry = np.sort(entry)
		p = 1. * np.arange(len(entry))/(len(entry) - 1)
		plt.plot(sort_entry, p, color = color)
	plt.xlim(X_LIM_MIN, X_LIM_MAX)
	plt.ylim(0.0, 1.0)
	plt.xlabel('Log2 Fold Change', fontsize = 16)
	plt.ylabel('Cumulative Fraction of Genes', fontsize = 16)
	plt.yticks(fontsize = 16)
	plt.xticks(fontsize = 16)
	plt.savefig(output_file)
	plt.close()


def perform_t_test(subset_gene_name, list_subset_fc, list_all_fc):
	t, p = stats.ttest_ind(list_subset_fc, list_all_fc, equal_var=True)
	print('#', subset_gene_name, "vs all protein-coding genes by Student's Two-Tail t-test p-value =", p)
	print('# Number of', subset_gene_name, '=', len(list_subset_fc))
	print('# Number of all protein-coding genes =', len(list_all_fc))


def main():
	c = 'g'

	# read in arguments
	file_input_gene_fc = sys.argv[1]
	file_input_subset_gene_list = sys.argv[2]
	file_output_pdf = sys.argv[3]
	subset_gene_name = sys.argv[4]
	if len(sys.argv) > 5:
		c = get_optional_arguments(sys.argv[5:])

	# generate list of fc for subset of genes and all protein-coding genes
	list_subset_fc, list_all_fc = make_fc_lists(file_input_gene_fc, file_input_subset_gene_list)

	# make ecdf and calculate p-value for the mean log2 FC in gene expression levels of 
	# subset of genes compared to all genes
	output = make_ecdf(file_output_pdf, (list_subset_fc, c), (list_all_fc, 'k'))
	perform_t_test(subset_gene_name, list_subset_fc, list_all_fc)


if __name__ == '__main__':
	main()
