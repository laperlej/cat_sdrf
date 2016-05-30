# -*- coding: utf-8 -*-
"""
analyser la colonne 5)antibody et 6)taget pour ajouter la cible ou le tag à la colonne clean_target
analyser la colonne 8)strain, 9)genotype, 11)description pour déterminer la protéine-cible avec les flag (regex +match avec dictionnaire de genes)
"""

import re
from collections import OrderedDict

def filter_rows(rows, target_dico, input_cols, output_col):
	for row in rows:
		row = filter_row(row, target_dico, input_cols, output_col) 
	return rows

def merge_cols(row, input_cols):
	return "".join(row[input_col].lower() for input_col in input_cols)

def filter_row(row, target_dico, input_cols, output_col):
	"""
	itère sur les regex de gene_dico (puis de target-dico) et les compare avec l'information dans input_col jusqu'à ce qu'il y ait un match.
	input: 
		row: dictionnaire, la clé est le titre de la colonne et la valeur est le contenu de la colonne
		target_dico: le dictionnaire de regex pour les cibles
		input_cols: liste des colonnes dans lesquelles on cherche (concaténées).
		output_col: colonne changée lorsqu'il y a un match
	output:
		row: colonne output_col = clé du dictionnaire (info) s'il y a eu un match
	"""
	
	#section qui itère sur target_dico	
	new_value = ""
	for info in target_dico.keys():
		searchtarget = merge_cols(row, input_cols)
		if re.search(target_dico[info], searchtarget):
			new_value = info
			break
	 
	row[output_col] = new_value
	return row


def assign_tag_multiple(rows, tag_dico, gene_dico, chip_dico):
	for row in rows:
		row["clean_target"] = assign_tag(row, tag_dico, gene_dico, chip_dico)
	return rows

def assign_tag(row, tag_dico, gene_dico, chip_dico):
	#s'il n'y a pas d'anticorps
	if "MNase" in row["clean_assay"] or "DNase" in row["clean_assay"]:
		return "N/A"

	elif "none" in row["clean_target"]:
	#	et qu'il y a input dans les col 4, 11 et 13 concaténées
		if "input" in merge_cols(row, ["4)assaytype", "11)description","13)other"]).lower():
			return "input"
		else:
			return row["clean_target"]		
	
	elif "tag" in row["clean_target"]:
	#	si untagged dans les col 4, 11 et 13
		if "untagged" in merge_cols(row, ["4)assaytype", "11)description","13)other"]).lower():
	#		ajouter Mock à clean target:
			return "Mock"
			
		elif "mock" in 	merge_cols(row, ["4)assaytype", "11)description","13)other"]).lower():
	#		ajouter Mock à clean target:
			return "Mock"

		elif "input" in merge_cols(row, ["4)assaytype", "11)description","13)other"]).lower():
			return "input" 	
		else:	
	#		si tagged, retourne fonction compare_tag	
			return compare_tag(row,tag_dico, gene_dico, chip_dico)
	else:
		return row["clean_target"]

def compare_tag(row, tag_dico, gene_dico, chip_dico): 
	tagged = row["clean_target"]
	#compare le regex du tag de clean_target à ce qui a dans les colonnes 8-9-11
	match = re.search(tag_dico[tagged],merge_cols(row,["6)target","8)strain", "9)genotype", "11)description"]).lower())
	if match:	
		#compare match avec gene_dico
		for gene in gene_dico:
			match2 = re.search(gene_dico[gene], match.group(1))		
			if match2:
				return gene
	else: 
		for regex in chip_dico.keys():
			match3 = re.search(chip_dico[regex], merge_cols(row,["11)description"]).lower())
			if match3:
				for gene in gene_dico:
					match4 = re.search(gene_dico[gene], match3.group(1))		
					if match4:
						return gene
	return row["clean_target"]

		
