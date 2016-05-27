""" finds the protein tagged """

import re
def tagfilter_rows(rows, tags_dico, input_cols, output_col):
	for row in rows:
		row = filter_row(row, tags_dico, input_cols, output_col) 
	return rows

def tagfilter_row(row, tags_dico, input_cols, output_col):
	"""
	itère sur les regex de target-dico et les compare avec l'information dans input_col jusqu'à ce qu'il y ait un match.
	input: 
		row: dictionnaire, la clé est le titre de la colonne et la valeur est le contenu de la colonne
		untagged_dico: dictionnaire de regex pour tagged ou untagged
		tags_dico: le dictionnaire de regex pour les tags
		genes_dico: dictionnaire de regex pour tous les gènes de S. cerevisiae
		input_cols1: liste des colonnes dans lesquelles on cherche (tagged ou untagged)
		input_cols2: liste des colonnes dans lesquelles on cherche (tag-gène)
		input_cols3: liste des colonnes dans lesquelles on cherche (les gènes)
		output_col: colonne changée lorsqu'il y a un match
	output:
		row: colonne output_col = match
	"""
	new_value = ""
	for info in tags_dico.keys():
		searchtarget = ''
		for input_col2 in input_cols2:
			searchtarget += row[input_col2].lower()
			match = re.search(tags_dico[info], searchtarget)
		if match:
			new_value = match
			print new_value
			for gene in genes_dico.keys():
				searchtarget2 = ''
				for input_col3 in input_cols3:
					searchtarget2 += row[input_col3].lower()
					match = re.search(genes_dico[])

			break	
	 
	row[output_col] = new_value
	return row



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
		searchtarget = ''
		for input_col in input_cols:
			searchtarget += row[input_col].lower()
		if re.search(target_dico[info], searchtarget):
			new_value = info
			break	
	 
	row[output_col] = new_value
	return row
