# -*- coding: utf-8 -*-
"""
analyser la colonne 5)antibody et 6)taget pour ajouter la cible ou le tag à la colonne clean_target
analyser la colonne 8)strain et 9)genotype pour déterminer la protéine-cible avec les flag (regex +match avec dictionnaire de genes)
"""

import re
def filter_rows(rows, target_dico, input_cols, output_col):
	for row in rows:
		row = filter_row(row, target_dico, input_cols, output_col) 
	return rows

def merge_cols(row, input_cols):
	return "".join(row[input_col].lower() for input_col in input_cols)

def filter_row(row, target_dico, input_cols, output_col):
	"""
	itère sur les regex de target-dico et les compare avec l'information dans input_col jusqu'à ce qu'il y ait un match.
	input: 
		row: dictionnaire, la clé est le titre de la colonne et la valeur est le contenu de la colonne
		target_dico: le dictionnaire de regex
		input_cols: liste des colonnes dans lesquelles on cherche (concaténées).
		output_col: colonne changée lorsqu'il y a un match
	output:
		row: colonne output_col = clé du target_dico (info) s'il y a eu un match
	"""
	new_value = ""
	for info in target_dico.keys():
		searchtarget = merge_cols(row, input_cols)
		if re.search(target_dico[info], searchtarget):
			new_value = info
			break
	 
	row[output_col] = new_value
	return row

def assign_tag(row, tag_dico):
	#	if "tag" in row[clean_target]:
	if "tag" in row[clean_target]:

	#		si <untagged< dans les col 4, 11 et 13
		if "Untagged" in merge_cols(row, ["4)assaytype", "11)description","13)other"]):
	#			ajouter Mock a clean target:
			row["clean_target"]= "Mock"
	#		sinon 
		else:	
	#			faitdequoi()	
			row["clean_target"] = compare_tag(row,tag_dico)


def compare_tag(row, tag_dico):
	if "tag" in row[clean_target]:
	#compare tag de clean_taget avec regex dans tag_dico
		for tagged in tag_dico.keys():
			match = re.search(tag_dico[tagged],"tag")
	#si match
			if match:		
	#compare match avec gene_dico
				for gene in gene_dico:
					match2 = re.search(gene_dico[gene], match)		
	#si match
					if match2:
	#ecrit gene dans clean_target
						return gene
		return ""		
