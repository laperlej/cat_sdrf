# -*- coding: utf-8 -*-
"""
analyser la colonne 5)antibody et 6)taget pour ajouter la cible ou le tag à la colonne clean_target
analyser la colonne 8)strain, 9)genotype, 11)description pour déterminer la protéine-cible avec les flag (regex +match avec dictionnaire de genes)
"""

import re
from collections import OrderedDict

# Itère sur chaque ligne du dictionnaire rows (contient les info des fichiers sdrf)
def filter_rows(rows, target_dico, input_cols, output_col):
	for row in rows:
		row = filter_row(row, target_dico, input_cols, output_col) 
	return rows

#Joint le contenu des colonnes (input_cols) sans séparateur
def merge_cols(row, input_cols):
	return "".join(row[input_col].lower() for input_col in input_cols)

def filter_row(row, target_dico, input_cols, output_col):
	"""
	itère sur les regex de target-dico et les compare avec l'information dans input_col jusqu'à ce qu'il y ait un match.
	Fonction multi-task, peut être utilisée aussi pour filtrer les essais 
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


def assign_tag_multiple(rows, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico):
	"""
	input: 
		tag_dico: dictionnaire de regex pour matcher la cible des tag (..... - tag)
		histone_dico: dictionnaire des marques d'histones avec regex
		gene_dico: dictionnaire de gènes avec regex
		gene_descrip_dico: dictionnaire de noms standards et des alias pour les gènes (non redondant)
		chip_dico: dictionaire de regex pour matcher la cible des chip (.... chip)
	"""

	for row in rows:
		row["clean_target"], row["reliability"] = assign_tag(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico)
	return rows
#section qui trie entre les input, mock, tag et types d'essai
#appelle compare_tag pour tag; 
def assign_tag(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico):
	#Assigne 'N/A' à la colonne clean_target si l'essai est mnase ou dnase
	if "MNase" in row["clean_assay"] or "DNase" in row["clean_assay"] or "FAIRE-seq" in row["13)other"] or "FAIRE-seq" in row["11)description"]:
		return "N/A", "assay type"
	#Assigne 'input' à la colonne clean_target si le mot-clé input est trouvé
	elif "input" in merge_cols(row, ["4)assaytype", "11)description","13)other"]).lower():
		return "input", "keyword in concat"
	#Assigne 'mock' à la colonne clean_target si le mot-clé 'mock' est trouvé
	elif "mock" in merge_cols(row, ["4)assaytype", "11)description","13)other"]).lower():
		return "Mock", "keyword in concat"
	
	elif "empty" in row["clean_target"]:
		#section qui itère sur histones_dico pour trouver un match sur deux niveaux, moins au plus large
		for hist in histones_dico.keys():
			if re.search(histones_dico[hist], merge_cols(row, ["1)identifier", "5)antibody", "6)target", "9)genotype", "11)description", "13)other"]).lower()):
				return hist, "histone mark"	
		#section qui itère sur le dictionnaire de gènes, à deux niveaux (moins au plus large)
		for gen in gene_dico.keys():
			if re.search(gene_dico[gen], row["11)description"].lower()):
				return gen, "gene"
			elif re.search(gene_dico[gen], merge_cols(row, ["1)identifier", "5)antibody", "6)target", "9)genotype", "13)other"]).lower()):
				return gen, "gene at large"
		#section qui itère sur le dictionnaire de nom standard et description de gènes, à deux niveaux (moins au plus large)
		for gen_descr in gene_descrip_dico.keys():
			if re.search(gene_descrip_dico[gen_descr], row["11)description"].lower()):
				return gen_descr, "gene descr"
			elif re.search(gene_descrip_dico[gen_descr], merge_cols(row, ["1)identifier", "5)antibody", "6)target", "9)genotype", "13)other"]).lower()):
				return gen_descr, "gene descr at large"	
		return row["clean_target"], "target_dico"	

	# section qui appelle compare_tag pour les lignes qui contiennent 'tag' et qui ne sont pas des 'Mock'
	elif "tag" in row["clean_target"]:
	#	si untagged dans les col 4, 11 et 13
		if "untagged" in merge_cols(row, ["4)assaytype", "11)description","13)other"]).lower() or "no tag" in merge_cols(row, ["4)assaytype", "11)description","13)other"]).lower():
	#		ajouter Mock à clean target:
			return "Mock","keyword in concat"
		else:	
	#		si tagged, retourne fonction compare_tag	
			return compare_tag(row,tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico)	
	else:
		return row["clean_target"], "target_dico"

def compare_tag(row, tag_dico, histones_dico, gene_dico, gene_descrip_dico, chip_dico): 
	#ajouter filtre avec histones_dico ici
	tagged = row["clean_target"]
	#compare le regex du tag de clean_target à ce qu'il y a dans les colonnes 4-6-8-9-11
	#idéalement, enlever ici la colonne strain, amène des faux positifs
	match = re.search(tag_dico[tagged],merge_cols(row,["4)assaytype", "6)target","8)strain", "9)genotype", "11)description"]).lower())
	if match:
		#compare match avec histones_dico
		for hist in histones_dico.keys():
			match2 = re.search(histones_dico[hist], match.group(1))		
			if match2:
				return hist, "tag to histone"
		#compare match avec gene_dico
		for gene in gene_dico.keys():
			match2 = re.search(gene_dico[gene], match.group(1))		
			if match2:
				return gene, "tag to gene"
		
		#compare match avec gene_descrip_dico
		for gene in gene_descrip_dico.keys():
			match2 = re.search(gene_descrip_dico[gene], match.group(1))		
			if match2:
				return gene, "tag to gene descr"
	else: 
		# itère sur les regex du chip_dico (pour trouver la cible des ChIP)
		for regex in chip_dico.keys():
			match3 = re.search(chip_dico[regex], merge_cols(row,["11)description"]).lower())
			if match3:
				#dictionnaire d'histones?
				for gene in gene_dico.keys():
					#compare match du regex avec gene_dico
					match4 = re.search(gene_dico[gene], match3.group(1))		
					if match4:
						return gene, "chip to gene"
				
				for gene in gene_descrip_dico.keys():
					#compare match du regex avec gene_descrip_dico
					match4 = re.search(gene_descrip_dico[gene], match3.group(1))		
					if match4:
						return gene, "chip to gene descr"	
									
	return row["clean_target"], "target_dico"

		
