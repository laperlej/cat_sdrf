# -*- coding: utf-8 -*-
"""
analyser la colonne 5)antibody et 6)taget pour ajouter la cible ou le tag à la colonne clean_target
analyser la colonne 8)strain, 9)genotype, 11)description pour déterminer la protéine-cible avec les flag (regex +match avec dictionnaire de genes)
"""

import re
from collections import OrderedDict

#Joins the content of columns (input_cols) with '|' as separator
def merge_cols(row, input_cols):
	return "|".join(row[input_col].lower() for input_col in input_cols if row[input_col] is not None)

"""
def raw_files_filter_rows(rows, output_col1, output_col2 ):
	for row in rows:
		row = raw_files_filter_row(row, output_col1, output_col2)
	return rows
#probably not necessary with chip-chip
def raw_files_filter_row(row, output_col1, output_col2):
	#Fills the 'raw_files column with .CEL, .pair, .GPR or .txt files
	raw_file_dict = {'CEL':'(\S+\.CEL\.gz)', 'PAIR':'(\S+\.pair\.gz)', 'GPR':'(\S+\.gpr\.gz)', 'TXT':'(\S+\.txt\.gz)'}
	new_value = ''
	for raw_file in raw_file_dict:
		#important to keep the column's upper case for the file's URL (not using 'merge_col' here)
		searchtarget = row['19)all_supp_files']
		match =  re.findall(raw_file_dict[raw_file], searchtarget)
		for pattern in match:
			new_value = " | ".join(match)
		
	row[output_col1] = new_value
	return row """
#iterates on each line of the dictionnary that is rows (metadatas from xml files) 
def condition_rows(rows, species):
	for row in rows:
		row = condition_row(row, species)
	return rows	
	
def condition_row(row, species):
	"""explanation for the selection in the «clean» or «discard» file 
	returns «yes» if all conditions are met in row; returns «wrong org», «no raw data» or «wrong assay» otherwise
	input_cols = 3)organism, 4)clean_assay, 18)raw_files
	output_col = «Selection»
	"""
	species_dict={
		"saccer": "Saccharomyces cerevisiae",
		"pombe": "Schizosaccharomyces pombe",
		"celegans":"Caenorhabditis elegans"}
	discard_assays=["rip-seq", "unwanted", 'non-genomic', 'wgs', 'faire', 'rna-seq', 'chip-chip', 'chip-seq', 'brdu', 'mnase-chip', 'dnase', 'other']
	
	if species_dict[species] not in row['3)organism']:
		row['Selection'] = 'Wrong org' 
	elif not row['18)raw_files']:
		row['Selection'] = 'No raw data'
	elif any(discard_assay in row['4)clean_assay'].lower() for discard_assay in discard_assays):
		row['Selection'] = 'Wrong assay'
	else:
		row['Selection'] = 'Yes'
	return row	
	
# Iterate on each line of the dictionnary that is rows (contains info from the xml files)
def filter_rows(rows, target_dico, histones_dico, input_cols, output_col):
	""" Concatenate the content of 2 dictionnaries (necessary when filtering for the antibody target; otherwise use an empty dictionnary plus the one needed"""
	all_targets = OrderedDict ([])
	all_targets.update(histones_dico)
	#this dict comes last since it ends with a regex catching anything (when not using an empty dict)
	all_targets.update(target_dico)
	for row in rows:
		row = filter_row(row, all_targets, input_cols, output_col) 
	return rows

def filter_row(row, all_targets, input_cols, output_col):
	"""
	iterate on the regex of target-dico and compares it to the information in input_col until a match is found.
	Multi-task function, can be used to filter the 'antibody' column and the 'assaytype' column (each with their own arguments) 
	input: 
		row: dictionnary where the key is the column's title and the value is the content of said column
		all_targets: target-regex dictionnary (concatenation of 2 dict)
		input_cols: list of the concatenated columns in which we search.
		output_col: column changed if there is a match
	output:
		row: column output_col = key of the dictionnary (info) if there is a match
	"""
	
	#Iterates on a target-regex dictionnary	
	new_value = ""
	for info in all_targets:
		# Defines the search target as the concatenation of some columns (by merge_cols; also converts to lowercase)
		searchtarget = merge_cols(row, input_cols)
		#If there is a match between the regex of the target and the concatenated columns	
		if re.search(all_targets[info], searchtarget):
			new_value = info
			break

	#Either a target (if one was found) or nothing will assigned to the output column
	row[output_col] = new_value
	return row


def assign_tag_multiple(rows, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico, antibody_dico):
	"""
	input: 
		tag_dico: regex dictionnary matching the protein (taget) linked with the tag (..... - tag)
		histone_dico: histone marks dictionnary (with regex)
		gene_dico: gene dictionnary (with regex)
		gene_descrip_dico: gene-alias dictionnary (not redundant with the gene dictionnary)
		chip_dico: regex dictionnary matching the target of the ChIP (.... chip or chip of ...)
		antibody_dico: antibodies' catalog number dictionnary
	"""
	#Function permiting the use of more CPU (accelerates the process)
	import multiprocessing
	multiprocessing.cpu_count()
	num_cpu = multiprocessing.cpu_count()
	try:	
		p = multiprocessing.Pool(num_cpu)
		arguments = ((row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico, antibody_dico) for row in rows)
		results = p.map(assign_tag_parallel, arguments)
		#zip returns a list of tuples
		for row, result in zip(rows, results):		
			row["5)clean_target"], row["6)reliability"] = result
		return rows	 		
	finally:
		p.close()
		p.join()
	exit(1)

def assign_tag_parallel(data):
	return assign_tag(*data)


def assign_tag(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico, antibody_dico):
	""" Function filtering inputs, mocks, tags and also according to the assay type;
		Makes use of other functions, such as compare_tag, compare_chip and compare_directly (different levels or comparison);
		This function will overwrite the content of 'clean_assay'
	"""
	#Assign 'N/A' to 5)clean_target column if the assay type is mnase, dnase, with ssDNA, bisulfite, with BrdU or FAIRE-Seq
	assays_list = ['mnase', 'dnase', 'faire', 'ssdna', 'bisulfite-seq', 'atac-seq', 'brdu']
	if any(assay in merge_cols(row,["4)clean_assay"]) for assay in assays_list):
		return "N/A", "assay type (1)"
	# With ChIP-chip, seems that BrdU assays have ideed no target
	#elif 'brdu' in row["4)clean_assay"].lower() and 'brdu' in row["8)antibody"].lower():
	#	return "N/A", "assay type (1)"	
	if 'whole cell extract' in merge_cols(row,["17)Sample_description"]) or 'whole-cell extract' in merge_cols(row,["17)Sample_description"]):
		return 'WCE', 'keyword (1)'
	elif 'not specified' in merge_cols(row,["8)antibody"])  and 'input' in merge_cols(row,["13)cell_type","17)Sample_description", "1,1)Sample_title"]):
		return 'input', 'keyword (1)'
	#Assign 'input' to column 'clean_target' if one of the following keyword is found in specific lines
	elif "input" in merge_cols(row, ["7)assaytype", "8)antibody"]) or "reference dna" in merge_cols(row, ["7)assaytype", "8)antibody"]):
		return "input", "keyword (1)"

	input_word_list = ['chromatin input', 'input sonicated dna', 'wce fraction used for the nomalization', 'wce fraction used for normalization', 'input lane', 'input dataset', 'channel ch1 is input dna', 'noip' ]
	if any(input_word in merge_cols(row, ["17)Sample_description", "1,1)Sample_title", '13)cell_type', '11)Material_type']) for input_word in input_word_list):
		return "input", "keyword (1)"

	elif 'input' in merge_cols(row, ["13)cell_type"]) and 'input control' not in merge_cols(row, ["17)Sample_description", "1,1)Sample_title"]) and 'input dna' in merge_cols(row, ["17)Sample_description", "1,1)Sample_title"]):
		return 'input', 'keyword (2)'	
	elif 'input dna' not in merge_cols(row, ["17)Sample_description", "1,1)Sample_title"]) and 'input control' not in merge_cols(row, ["17)Sample_description", "1,1)Sample_title"]) and 'input' in merge_cols(row, ["17)Sample_description", "1,1)Sample_title"]) and 'compared to input' not in merge_cols(row, ["17)Sample_description"]):
		return 'input', 'keyword (3)'
#	elif '[input dna]' in merge_cols(row, ["17)Sample_description"]):
	elif 'input dna' in merge_cols(row, ["17)Sample_description", "1,1)Sample_title"]):
		return 'input', 'keyword (3)'
	#Assign 'mock' to column 'clean_target' if one of the following keyword is found
	mock_list = ['non antibody control', 'no epitope tag', 'no-epitope', 'untagged', 'un-tagged', 'no tag', 'notag', 'no tap tag', 'null-tap', 'no-tag']
	#added the column 'all_supp_files', sometimes a keyword is found in the file name
	if any(mock in merge_cols(row, [ '9)target' ,"7)assaytype", "17)Sample_description", "1,1)Sample_title", '15)genotype', '14)strain', '19)all_supp_files']) for mock in mock_list):
		return "Mock", "keyword (1)"
	# special case were 'IP-vs-mock' is a sample name
	elif 'mock' in merge_cols(row, [ '9)target' ,"7)assaytype", "17)Sample_description", "1,1)Sample_title", '15)genotype', '14)strain', '19)all_supp_files']) and 'ip-vs-mock' not in merge_cols(row, [ '9)target' ,"7)assaytype", "17)Sample_description", "1,1)Sample_title", '15)genotype', '14)strain', '19)all_supp_files']):
		return "Mock", "keyword (1)"
	#Assign 'control' to column 'clean_target' if one of the following keyword is found
	control_list = ['control for', 'control_for', 'control replicate', 'degron', 'wild type control']
	if any(control in merge_cols(row, ["7)assaytype", "17)Sample_description", "1,1)Sample_title","Other"]) for control in control_list):
		return "control", "keyword (1)"
	elif 'control' in  merge_cols(row, ["11)Material_type", "17)Sample_description", "1,1)Sample_title"]) and 'input control' not in merge_cols(row, ["17)Sample_description", "1,1)Sample_title"]) and 'mock' not in merge_cols(row, ["5)clean_target"]):
		return "control", "keyword (2)"

	elif "empty" in row["5)clean_target"]:
		#Calls to different functions in order to find a target
		var_search_target = search_target(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico)
		if var_search_target is not None:
			return var_search_target
		var_search_antibody = search_antibody(row, antibody_dico)
		if var_search_antibody is not None:
			return var_search_antibody	
		var_compare_tag2 = 	compare_tag2(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico)
		if var_compare_tag2 is not None:
			return var_compare_tag2
		var_compare_chip1 = compare_chip1(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico)
		if var_compare_chip1 is not None:
			return var_compare_chip1
		var_compare_tag_larger2 = compare_tag_larger2(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico)
		if var_compare_tag_larger2 is not None:
			return var_compare_tag_larger2	
		var_compare_chip2 = compare_chip2(row, histones_dico, gene_dico, gene_descrip_dico, chip_dico)
		if var_compare_chip2 is not None:
			return var_compare_chip2
		var_compare_directly =	compare_directly(row, histones_dico, gene_dico, gene_descrip_dico)
		if var_compare_directly is not None:
			return var_compare_directly	
				
	elif "tag" in row["5)clean_target"]:
		#Calls to different functions in order to find a target
		var_search_target = search_target(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico)
		if var_search_target is not None:
			return var_search_target
		var_search_antibody = search_antibody(row, antibody_dico)
		if var_search_antibody is not None:
			return var_search_antibody	
		var_compare_tag1 = 	compare_tag1(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico)
		if var_compare_tag1 is not None:
			return var_compare_tag1
		var_compare_chip1 = compare_chip1(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico)
		if var_compare_chip1 is not None:
			return var_compare_chip1
		var_compare_tag_larger1 = compare_tag_larger1(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico)
		if var_compare_tag_larger1 is not None:
			return var_compare_tag_larger1
		var_compare_chip2 = compare_chip2(row, histones_dico, gene_dico, gene_descrip_dico, chip_dico)
		if var_compare_chip2 is not None:
			return var_compare_chip2
		var_compare_directly =	compare_directly(row, histones_dico, gene_dico, gene_descrip_dico)
		if var_compare_directly is not None:
			return var_compare_directly

	else:
		return row["5)clean_target"], "target_dico (1)"

def search_target(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico):
	""" This function searches for a target that is a gene by comparing the content of the gene dictionnary with the columns 'assaytype', 'antibody' and 'target'"""
	for gene in gene_dico:
		if re.search(gene_dico[gene],merge_cols(row,["7)assaytype", "8)antibody", "9)target"])):
			return gene, "target (2)"
	return None	

def search_antibody(row, antibody_dico):
	""" This function searches for antibodies' catalog number in columns 'cell_type', 'description' and 'Sample_title"""
	for antibody in antibody_dico:
		if re.search(antibody_dico[antibody],merge_cols(row,['13)cell_type', "17)Sample_description"])):
			return antibody, "antibody no (2)"	
		elif re.search(antibody_dico[antibody],merge_cols(row,[ "1,1)Sample_title"])):	
			return antibody, "antibody no (3)"

def compare_tag1(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico): 
	""" This function searches the target of the tag found in the columns 'antibody' or 'target'; very specific """
	#the tag found in clean_target determines which regex is used to find the target (gain in specificity)
	tagged = row["5)clean_target"]
	#compares the regex to the content of columns 'assaytype', 'cell_type' and 'descriptiom'
	match = re.search(tag_dico[tagged],merge_cols(row,["7)assaytype",'13)cell_type', "17)Sample_description"]))
	match2 = re.search(tag_dico[tagged],merge_cols(row,["1,1)Sample_title"]))
	if match:
		#compares tag's match with the histone dict
		for hist in histones_dico:
			if re.search(histones_dico[hist], match.group(1)):	
				return hist, "tag to histone (2)"
		#compares tag's match with the gene dict
		for gene in gene_dico:
			if re.search(gene_dico[gene], match.group(1)):
				return gene, "tag to gene (3)"
		
		#compares tag's match with the alias dict
		for gene in gene_descrip_dico:
			if re.search(gene_descrip_dico[gene], match.group(1)):	
				return gene, "tag to gene descr (4)"			
	elif match2:
		#compares tag's match with the histone dict
		for hist in histones_dico:
			if re.search(histones_dico[hist], match2.group(1)):	
				return hist, "tag to histone (2)"
		#compares tag's match with the gene dict
		for gene in gene_dico:
			if re.search(gene_dico[gene], match2.group(1)):
				return gene, "tag to gene (3)"
		
		#compares tag's match with the alias dict
		for gene in gene_descrip_dico:
			if re.search(gene_descrip_dico[gene], match2.group(1)):	
				return gene, "tag to gene descr (4)"
	return None 

def compare_chip1(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico): 
	""" This function searches for the keyword 'ChIP' and compares what comes before 'ChIP' with the histone, gene and alias dictionnaries"""
	#Iterates on the chip_dico regex (in order to find the target of the ChIP)
	for regex in chip_dico:
		match = re.search(chip_dico[regex], merge_cols(row,['13)cell_type', "17)Sample_description"]))
		match2 = re.search(chip_dico[regex], merge_cols(row,[ "1,1)Sample_title"]))
		if match:
			for hist in histones_dico:
				#compares regex's match with the histone dict
				if re.search(histones_dico[hist], match.group(1)):
					return hist, "chip to histone (3)"
				
			for gene in gene_dico:
				#compares regex's match with the gene dict
				if re.search(gene_dico[gene], match.group(1)):
					return gene, "chip to gene (3)"	
				
			for gene in gene_descrip_dico:
				#compares regex's match with the alias dict
				if re.search(gene_descrip_dico[gene], match.group(1)):
					return gene, "chip to gene descr (4)"
	if match2:
			for hist in histones_dico:
				#compares regex's match with the histone dict
				if re.search(histones_dico[hist], match2.group(1)):
					return hist, "chip to histone (3)"
				
			for gene in gene_dico:
				#compares regex's match with the gene dict
				if re.search(gene_dico[gene], match2.group(1)):
					return gene, "chip to gene (3)"	
				
			for gene in gene_descrip_dico:
				#compares regex's match with the alias dict
				if re.search(gene_descrip_dico[gene], match2.group(1)):
					return gene, "chip to gene descr (4)"
	return None 

def compare_tag_larger1(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico): 
	""" This function searches the target of the tag found in the columns 'antibody' or 'target'; columns 'strain' and 'genotype' are sometimes very specific and sometimes very unspecific (when it lists many genes) """
	#the tag found in clean_target determines which regex is used to find the target (gain in specificity)
	tagged = row["5)clean_target"]
	#compares the regex to the content of columns 'strain' and 'genotype' 
	match = re.search(tag_dico[tagged],merge_cols(row,["14)strain", "15)genotype"]))
	if match:
		#compares tag's match with the histone dict
		for hist in histones_dico:
			if re.search(histones_dico[hist], match.group(1)):	
				return hist, "tag to histone (3)"
		#compares tag's match with the gene dict
		for gene in gene_dico:
			if re.search(gene_dico[gene], match.group(1)):
				return gene, "tag to gene (4)"
		
		#compares tag's match with the alias dict
		for gene in gene_descrip_dico:
			if re.search(gene_descrip_dico[gene], match.group(1)):	
				return gene, "tag to gene descr (4)"
	return None


def compare_tag2(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico):
	""" This function searches for any tags that were not found with columns 'antibody' and 'target'=> for 'empty' in 5)clean_target; quite specific search"""
	for tag in tag_dico:
		#compares the regex (value as (something)::flag) to the content of the columns 'cell_type' and 'Sample_title'
		match = re.search(tag_dico[tag],merge_cols(row,[ '13)cell_type', "17)Sample_description"]))
		match2 = re.search(tag_dico[tag],merge_cols(row,[ "1,1)Sample_title"]))
		if match:
			#compares tag's match with the histone dict
			for hist in histones_dico:
				if re.search(histones_dico[hist], match.group(1)):	
					return hist, "search tag to histone (2)"
			#compares tag's match with the gene dict
			for gene in gene_dico:
				if re.search(gene_dico[gene], match.group(1)):
					return gene, "search tag to gene (4)"
			#compares tag's match with the alias dict
			for gene in gene_descrip_dico:
				if re.search(gene_descrip_dico[gene], match.group(1)):	
					return gene, "search tag to gene descr (4)"
		if match2:
			#compares tag's match with the histone dict
			for hist in histones_dico:
				if re.search(histones_dico[hist], match2.group(1)):	
					return hist, "search tag to histone (2)"
			#compares tag's match with the gene dict
			for gene in gene_dico:
				if re.search(gene_dico[gene], match2.group(1)):
					return gene, "search tag to gene (4)"
			#compares tag's match with the alias dict
			for gene in gene_descrip_dico:
				if re.search(gene_descrip_dico[gene], match2.group(1)):	
					return gene, "search tag to gene descr (4)"			
	return None

def compare_tag_larger2(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico):
	""" This function searches for any tags that were not found with columns 'antibody' and 'target'=> for 'empty' in 5)clean_target; larger search"""
	for tag in tag_dico:
		#compares the regex (value as (something)::flag) to the content of the columns 'strain' and 'genotype'
		match = re.search(tag_dico[tag],merge_cols(row,[ "14)strain", "15)genotype"]))
		if match:
			#compares tag's match with the histone dict
			for hist in histones_dico:
				if re.search(histones_dico[hist], match.group(1)):	
					return hist, "search tag to histone (2)"
			#compares tag's match with the gene dict
			for gene in gene_dico:
				if re.search(gene_dico[gene], match.group(1)):
					return gene, "search tag to gene (4)"
			#compares tag's match with the alias dict
			for gene in gene_descrip_dico:
				if re.search(gene_descrip_dico[gene], match.group(1)):	
					return gene, "search tag to gene descr (4)"				
	return None		

def compare_chip2(row, histones_dico, gene_dico, gene_descrip_dico, chip_dico):
	""" This function searches for the keyword 'ChIP' and compares what comes before 'ChIP' with the histone, gene and alias dictionnaries"""
	#Iterates on the chip_dico regex (in order to find the target of the ChIP)
	for regex in chip_dico:
		#match = re.search(chip_dico[regex], merge_cols(row,["1)identifier", "1,1)Sample_title", "14)strain", "15)genotype"]))
		match = re.search(chip_dico[regex], merge_cols(row,["1,1)Sample_title", "14)strain", "15)genotype"]))
		if match:
			for hist in histones_dico:
				#compares the regex's match with the histone dictionnary
				if re.search(histones_dico[hist], match.group(1)):
					return hist, "chip to histone (4)"
				
			for gene in gene_dico:
				#compares the regex's match with the gene dictionnary
				if re.search(gene_dico[gene], match.group(1)):
					return gene, "chip to gene (4)"	
				
			for gene in gene_descrip_dico:
				#compares the regex's match with the alias dictionnary
				if re.search(gene_descrip_dico[gene], match.group(1)):
					return gene, "chip to gene descr (4)"
	return None			


def compare_directly(row, histones_dico, gene_dico, gene_descrip_dico):
	""" Compares the content of certain columns with dictionnaries in sequence (histone, gene, alias)"""
	#1rst sweep, most specific (normally)
	#iterate on the histone dictionnary and compares with specific column
	for hist in histones_dico:
		#if re.search(histones_dico[hist], merge_cols(row, ["1)identifier", "8)antibody", "9)target", "17)Sample_description"])):
		if re.search(histones_dico[hist], merge_cols(row, ["8)antibody", "9)target", "17)Sample_description"])):
			return hist, "histone mark (3)"
		elif re.search(histones_dico[hist], row["15)genotype"].lower()):
			return hist, "histone mark (4)"
	#iterates on the gene dictionnary and compares with specific column
	for gene in gene_dico:
		if re.search(gene_dico[gene], merge_cols(row, ["8)antibody","9)target", "17)Sample_description"])):
			return gene, "gene (4)"
	#iterates on the alias dictionnary and compares with specific columns 
	for gene in gene_descrip_dico:
		if re.search(gene_descrip_dico[gene], merge_cols(row, ["8)antibody","9)target", "17)Sample_description"])):
			return gene, "gene descr (5)"
	#2nd sweep, less specific
	#iterate on the gene dictionnary and compares with less specific columns
	for gene in gene_dico:
		if re.search(gene_dico[gene], row['1,1)Sample_title'].lower()):
			return gene, "gene (5)"	
		elif re.search(gene_dico[gene], merge_cols(row,[ "14)strain"])):
			return gene, "gene (5)"
		elif re.search(gene_dico[gene], row["15)genotype"].lower()):
			return gene, "gene (5)"		
	#iterate on the alias dictionnary and compares with less specific columns
	for gene in gene_descrip_dico:
		if re.search(gene_descrip_dico[gene], row['1,1)Sample_title'].lower()):
			return gene, "gene descr (5)"
		elif re.search(gene_descrip_dico[gene], merge_cols(row,[ "14)strain"])):
			return gene, "gene descr (5)"
		elif re.search(gene_descrip_dico[gene], row["15)genotype"].lower()):
			return gene, "gene descr (5)"

	#3rd sweep, least specific
	#iterate on the histone dictionnary
	for hist in histones_dico:
		if re.search(histones_dico[hist], row["Other"].lower()):
			return hist, "histone mark (5)"
	#iterate on the gene dictionnary and compares with least specific columns
	for gene in gene_dico:
		#if re.search(gene_dico[gene], merge_cols(row, ["1)identifier", "Other"])):
		if re.search(gene_dico[gene], merge_cols(row, ["Other"])):
			return gene, "gene (6)"
	#iterate on the alias dictionnary and compares with least specific columns
	for gene in gene_descrip_dico:
		#if re.search(gene_descrip_dico[gene], merge_cols(row, ["1)identifier", "Other"])):
		if re.search(gene_descrip_dico[gene], merge_cols(row, ["Other"])):
			return gene, "gene descr (6)"		

	return row["5)clean_target"], "target_dico (1)"


	
